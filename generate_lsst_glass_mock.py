#!/usr/bin/env python
"""generate_lsst_glass_mock.py — generate an LSST Y1/Y10 lens+source GLASS mock catalogue.

Self-contained CLI script: cosmology (CAMB) -> matter Cls (disk-cached) -> lens/source
n(z) tomography via binny's LSST survey preset -> optional theory 3x2pt Cls (gg/gG/GG,
saved as a SACC file) -> a lognormal GLASS simulation streamed straight to two HDF5
catalogues (source, lens), with a survey mask applied both to galaxy sampling and to a
per-galaxy SYST_WEIGHT (1/visibility). Diagnostic plots are opt-in PNGs (no interactive
display).

Examples:
    # full, physically valid Y1 mock (default cosmology: Planck 2018, default density: real survey density)
    python generate_lsst_glass_mock.py --nside 1024 --year 1 --outdir out

    # fast dev-loop smoke test: shrink the galaxy count directly (--nside does NOT do this)
    python generate_lsst_glass_mock.py --nside 64 --year 1 --outdir out --debug-density-scale 0.01

    # Y10 mock, with theory Cls (SACC) and diagnostic plots, tagged with a fixed mock id
    python generate_lsst_glass_mock.py --nside 1024 --year 10 --outdir out \\
        --save-theory-cls 1 --save-diagnostics 1 --mock-id run042

    # custom cosmology via a small YAML file (only the keys you set are overridden from Planck18),
    # with an individual CLI override layered on top of that
    python generate_lsst_glass_mock.py --nside 1024 --year 1 --outdir out \\
        --cosmo-config my_cosmology.yaml --H0 70.0

    # full CAMB ini for total control of the parameter space
    python generate_lsst_glass_mock.py --nside 1024 --year 1 --outdir out --camb-ini my_params.ini
"""
from __future__ import annotations

import argparse
import hashlib
import logging
import os
import pickle
import secrets
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import healpy as hp
from tqdm import tqdm

import camb
from cosmology.compat.camb import Cosmology

import glass
import glass.ext.camb

import binny

log = logging.getLogger("glass_lsst_mock")

BIAS = 1.5           # constant linear galaxy bias applied to all lens tomobins
N_SHELLS = 25         # number of comoving-distance shells (matches the source notebook)

# Bundled survey mask: an LSST Y10 visibility map (0-1, fraction of the pixel observed),
# already inverted to visibility convention — do NOT re-invert it. Resolved relative to
# this script's own location, not the caller's cwd, so it works from anywhere.
DEFAULT_MASK_PATH = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "data", "lsst_y10_final_map4096.fits.gz"
)

# Planck 2018 (TT,TE,EE+lowE+lensing), Planck 2018 VI, Table 2.
PLANCK18_H0 = 67.36
PLANCK18_OMCH2 = 0.1200
PLANCK18_OMBH2 = 0.02237


# ────────────────────────────── logging / caching ──────────────────────────────

def setup_logging(outdir: str, level: str) -> None:
    os.makedirs(outdir, exist_ok=True)
    logging.basicConfig(
        level=getattr(logging, level.upper()),
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        handlers=[
            logging.StreamHandler(),
            logging.FileHandler(os.path.join(outdir, "run.log")),
        ],
    )


def cached(cache_dir: str, name: str, key_parts, compute_fn):
    os.makedirs(cache_dir, exist_ok=True)
    key = hashlib.sha256(repr(key_parts).encode()).hexdigest()[:16]
    path = os.path.join(cache_dir, f"{name}_{key}.pkl")
    if os.path.exists(path):
        log.info("cache hit for %s -> loading %s", name, path)
        with open(path, "rb") as f:
            return pickle.load(f)
    log.info("cache miss for %s -> computing", name)
    result = compute_fn()
    with open(path, "wb") as f:
        pickle.dump(result, f)
    return result


def _file_hash(path):
    """sha256 of a file's bytes, or None if path is falsy — used to invalidate disk
    caches when a --cosmo-config/--camb-ini file's content changes, without needing to
    enumerate every possible CAMB field it might set."""
    if not path:
        return None
    with open(path, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


# ────────────────────────────── cosmology / shells ──────────────────────────────

def load_cosmo_config(path):
    import yaml
    with open(path) as f:
        return yaml.safe_load(f) or {}


def resolve_cosmology(args):
    """Build (pars, cosmo) from a CAMB ini, a YAML config layered over Planck 2018
    defaults, or plain Planck 2018 defaults — then apply any explicit --H0/--omch2/
    --ombh2 CLI overrides last, on top of whichever path built `pars`."""
    if args.camb_ini:
        pars = camb.read_ini(args.camb_ini)
        log.info("cosmology loaded from CAMB ini: %s", args.camb_ini)
    else:
        params = {"H0": PLANCK18_H0, "omch2": PLANCK18_OMCH2, "ombh2": PLANCK18_OMBH2}
        if args.cosmo_config:
            params.update(load_cosmo_config(args.cosmo_config))
            log.info("cosmology defaults overridden from %s", args.cosmo_config)
        pars = camb.set_params(NonLinear=camb.model.NonLinear_both, **params)

    for name in ("H0", "omch2", "ombh2"):
        val = getattr(args, name)
        if val is not None:
            setattr(pars, name, val)
            log.info("cosmology override: %s=%s", name, val)

    results = camb.get_background(pars)
    cosmo = Cosmology(results)
    return pars, cosmo


def cosmology_cache_key(pars, args):
    """Cache-key component covering the resolved cosmology: H0/omch2/ombh2 read back
    off `pars` (covers the plain-default and CLI-override cases), plus a hash of the
    --cosmo-config/--camb-ini file's content (covers everything else such a file might
    set, e.g. ns/As, without enumerating every CAMB field)."""
    return (
        pars.H0, pars.omch2, pars.ombh2,
        _file_hash(args.camb_ini), _file_hash(args.cosmo_config),
    )


def build_shells(cosmo, zmax: float, num: int = N_SHELLS):
    zb = glass.distance_grid(cosmo, 0.0, zmax, num=num)
    shells = glass.linear_windows(zb)
    return zb, shells


def get_matter_cls(pars, lmax, shells, cache_dir, cosmo_key, zb):
    return cached(
        cache_dir,
        "matter_cls",
        (*cosmo_key, lmax, zb.tolist()),
        lambda: glass.ext.camb.matter_cls(pars, lmax, shells),
    )


# ────────────────────────────── tomography (binny) ──────────────────────────────

class Tomography:
    """Per-role (lens/source) n(z) tomography, sourced from binny's LSST preset."""

    def __init__(self, role, year, density_scale=1.0):
        tomo = binny.NZTomography()
        res = tomo.build_survey_bins("lsst", role=role, year=str(year), include_tomo_metadata=True)
        fractions = tomo.population_stats()["fractions"]  # int-keyed by bin index

        n_gal_arcmin2 = res.spec["sample_properties"]["number_density"]["n_gal_arcmin2"] * density_scale
        sigma_z0 = res.spec["uncertainties"]["scatter_scale"]
        shape_noise = res.spec["sample_properties"].get("shape_noise") or {}
        sigma_e = shape_noise.get("sigma_e")

        self.role = role
        self.year = str(year)
        self.z = np.asarray(res.z)
        self.nz = np.asarray(res.nz) * n_gal_arcmin2
        self.n_bins = len(res.bin_keys)
        self.tomo_nz = np.array(
            [res.bins[k] * n_gal_arcmin2 * fractions[k] for k in res.bin_keys]
        )
        self.bin_edges = np.asarray(res.tomo_meta["bins"]["bin_edges"])
        self.n_gal_arcmin2 = n_gal_arcmin2
        self.sigma_z0 = sigma_z0
        self.sigma_e = sigma_e

        log.info(
            "%s Y%s tomography: n_bins=%d n_gal_arcmin2=%.4f sigma_z0=%.3f%s%s",
            role, self.year, self.n_bins, n_gal_arcmin2, sigma_z0,
            f" sigma_e={sigma_e:.3f}" if sigma_e is not None else "",
            f" [DENSITY SCALED x{density_scale:g} — not a physical mock]" if density_scale != 1.0 else "",
        )


# ────────────────────────────── theory 3x2pt Cls (SACC) ──────────────────────────────

def compute_theory(pars, lmax, cls, lens: Tomography, src: Tomography, bias, cache_dir,
                    outdir, year, mock_id, save_diagnostics, cosmo_key):
    """Compute the GLASS/CAMB 3x2pt theory Cls and save them as a single SACC file.

    Tracer names (lens{i}/src{j}) and Cl data types (galaxy_density_cl / galaxy_shear_cl_ee /
    galaxy_shearDensity_cl_e, shear-tracer-first for the mixed pair) follow the Firecrown
    convention already used by cosmotron's own SACC tooling, so this file is directly
    consumable by that ecosystem without renaming anything.
    """
    import sacc

    lens_tomo_ngal = glass.partition(lens.z, lens.tomo_nz, shells_global)
    src_tomo_ngal = glass.partition(src.z, src.tomo_nz, shells_global)
    src_tomo_lens = glass.multi_plane_weights(src_tomo_ngal, shells_global, cosmo_global)

    cl_pos = glass.effective_cls(cls, lens_tomo_ngal)
    cl_she = glass.effective_cls(cls, src_tomo_lens)
    cl_cross = glass.effective_cls(cls, lens_tomo_ngal, src_tomo_lens)

    pars.Accuracy.TimeStepBoost = 5
    pars.SourceTerms.counts_density = True
    pars.SourceTerms.counts_redshift = False
    pars.SourceTerms.counts_lensing = False
    pars.SourceTerms.counts_velocity = False
    pars.SourceTerms.counts_radial = False
    pars.SourceTerms.counts_timedelay = False
    pars.SourceTerms.counts_ISW = False
    pars.SourceTerms.counts_potential = False
    pars.SourceTerms.counts_evolve = False
    pars.min_l = 1
    pars.set_for_lmax(lmax)
    pars.SourceWindows = (
        [camb.sources.SplinedSourceWindow(z=lens.z, W=lens.tomo_nz[i], source_type="counts")
         for i in range(lens.n_bins)]
        + [camb.sources.SplinedSourceWindow(z=src.z, W=src.tomo_nz[i], source_type="lensing")
           for i in range(src.n_bins)]
    )

    theory = cached(
        cache_dir,
        "camb_theory_dict",
        (
            *cosmo_key, lmax, lens.n_bins, src.n_bins,
            lens.tomo_nz.tobytes(), src.tomo_nz.tobytes(),
            tuple(lens.bin_edges.tolist()), tuple(src.bin_edges.tolist()),
        ),
        lambda: camb.get_results(pars).get_source_cls_dict(lmax=lmax, raw_cl=True),
    )

    def key_pos(i):
        return f"W{i + 1}"

    def key_she(i):
        return f"W{lens.n_bins + i + 1}"

    ell = np.arange(lmax + 1)
    lens_nbar = np.array([np.trapezoid(lens.tomo_nz[i], lens.z) for i in range(lens.n_bins)])
    glass_pos = cl_pos / lens_nbar[:, None, None] / lens_nbar[None, :, None]
    glass_cross = cl_cross / lens_nbar[:, None, None]
    glass_she = cl_she

    gG_pairs = [(il, js) for il in range(lens.n_bins) for js in range(src.n_bins)]

    # the catalogue is sampled with linear galaxy bias `bias`, but the GLASS/CAMB theory
    # spectra above are the UNBIASED matter spectra, so apply the galaxy bias here:
    #   gg -> bias**2, gG -> bias**1, GG -> bias**0
    save_gg = bias**2 * glass_pos
    save_gG = bias * glass_cross
    save_GG = glass_she

    ell_sacc = ell[1:].astype(float)

    s = sacc.Sacc()
    for i in range(lens.n_bins):
        s.add_tracer("NZ", f"lens{i}", lens.z, lens.tomo_nz[i])
    for j in range(src.n_bins):
        s.add_tracer("NZ", f"src{j}", src.z, src.tomo_nz[j])

    for i in range(lens.n_bins):
        for j in range(i, lens.n_bins):
            s.add_ell_cl("galaxy_density_cl", f"lens{i}", f"lens{j}", ell_sacc, save_gg[i, j, 1:])
    for i in range(src.n_bins):
        for j in range(i, src.n_bins):
            s.add_ell_cl("galaxy_shear_cl_ee", f"src{i}", f"src{j}", ell_sacc, save_GG[i, j, 1:])
    for i in range(lens.n_bins):
        for j in range(src.n_bins):
            # shear tracer first, per the Firecrown/cosmotron convention for a mixed pair
            s.add_ell_cl("galaxy_shearDensity_cl_e", f"src{j}", f"lens{i}", ell_sacc, save_gG[i, j, 1:])

    sacc_path = os.path.join(outdir, f"lsst_y{year}_theory_glass_mock_{mock_id}.sacc.h5")
    s.save_hdf5(sacc_path, overwrite=True)
    log.info("wrote theory 3x2pt Cls (SACC) to %s", sacc_path)

    mask = (ell[1:] >= 10) & (ell[1:] <= 100)
    ratios = {}
    for i in range(lens.n_bins):
        r_gg = np.median((glass_pos[i, i, 1:] / theory[f"{key_pos(i)}x{key_pos(i)}"][1:])[mask])
        ratios[f"lens_gg_bin{i}"] = float(r_gg)
    for i in range(src.n_bins):
        r_GG = np.median((glass_she[i, i, 1:] / theory[f"{key_she(i)}x{key_she(i)}"][1:])[mask])
        ratios[f"src_GG_bin{i}"] = float(r_GG)
    for (il, js) in gG_pairs:
        r_gG = np.median((glass_cross[il, js, 1:] / theory[f"{key_pos(il)}x{key_she(js)}"][1:])[mask])
        ratios[f"gG_lens{il}_src{js}"] = float(r_gG)
    log.info("GLASS/CAMB median ratios (ell 10-100): %s",
             ", ".join(f"{k}={v:.3f}" for k, v in ratios.items()))

    if save_diagnostics:
        diag_dir = os.path.join(outdir, "diagnostics")
        os.makedirs(diag_dir, exist_ok=True)
        fig, axes = plt.subplots(1, 3, figsize=(15, 4.2), sharey=True)
        for ax, t in zip(axes, ["clustering  gg", "galaxy-galaxy lensing  gG", "cosmic shear  GG"]):
            ax.axhline(1.0, c="k", ls="--", lw=0.8)
            ax.set_title(t)
            ax.set_xscale("log")
            ax.set_xlabel(r"$\ell$")
            ax.set_ylim(0.7, 1.3)
        axes[0].set_ylabel("GLASS / CAMB")
        for i in range(lens.n_bins):
            axes[0].plot(ell[1:], glass_pos[i, i, 1:] / theory[f"{key_pos(i)}x{key_pos(i)}"][1:],
                         label=f"bin {i}")
        for i in range(src.n_bins):
            axes[2].plot(ell[1:], glass_she[i, i, 1:] / theory[f"{key_she(i)}x{key_she(i)}"][1:],
                         label=f"bin {i}")
        for (il, js) in gG_pairs:
            axes[1].plot(ell[1:], glass_cross[il, js, 1:] / theory[f"{key_pos(il)}x{key_she(js)}"][1:],
                         label=f"lens {il} x src {js}")
        axes[0].legend()
        axes[1].legend(ncol=2, fontsize=7)
        plt.tight_layout()
        fig.savefig(os.path.join(diag_dir, "theory_ratio.png"), dpi=120)
        plt.close(fig)

    return ratios


# ────────────────────────────── mask ──────────────────────────────

def load_mask(mask_path: str, nside: int, outdir: str, save_diagnostics: bool):
    # the bundled mask is already a visibility map (0-1) -- do not invert it
    y_mask = hp.read_map(mask_path)
    y_ud = hp.ud_grade(y_mask, nside)

    mean_vis = float(np.mean(1.0 / y_mask[y_mask != 0]))
    log.info("mask loaded from %s: mean(1/visibility) over unmasked pixels = %.4f", mask_path, mean_vis)

    if save_diagnostics:
        diag_dir = os.path.join(outdir, "diagnostics")
        os.makedirs(diag_dir, exist_ok=True)
        hp.mollview(y_ud, title="LSST Forecast Mask", unit="Visibility")
        plt.savefig(os.path.join(diag_dir, "mask.png"), dpi=120)
        plt.close("all")

    hp.write_map(os.path.join(outdir, f"survey_footprint_N{nside}.fits"), y_ud, overwrite=True)
    return y_mask, y_ud


# ────────────────────────────── HDF5 streaming catalogue ──────────────────────────────

def make_h5_catalogue(path, fields):
    import h5py
    f = h5py.File(path, "w")
    for name, dtype in fields:
        f.create_dataset(name, shape=(0,), maxshape=(None,), dtype=dtype, chunks=True)
    return f


def append_h5(f, **columns):
    n_new_rows = len(next(iter(columns.values())))
    n_old = f[next(iter(columns))].shape[0]
    n_new = n_old + n_new_rows
    for name, arr in columns.items():
        f[name].resize((n_new,))
        f[name][n_old:n_new] = arr


def run_simulation(cls, shells, cosmo, src: Tomography, lens: Tomography, bias, mask_ud,
                    nside, rng, outdir, year, mock_id):
    fields = glass.lognormal_fields(shells)
    gls = glass.solve_gaussian_spectra(fields, cls)
    matter = glass.generate(fields, gls, nside, ncorr=3, rng=rng)
    convergence = glass.MultiPlaneConvergence(cosmo)

    src_ngal = glass.partition(src.z, src.nz, shells)
    lens_ngal = glass.partition(lens.z, lens.nz, shells)

    source_fields = [
        ("RA", "f8"), ("DEC", "f8"), ("PHOTO_Z", "f8"), ("ZBIN", "i8"),
        ("G1", "f8"), ("G2", "f8"), ("SYST_WEIGHT", "f8"),
    ]
    lens_fields = [
        ("RA", "f8"), ("DEC", "f8"), ("PHOTO_Z", "f8"), ("ZBIN", "i8"),
        ("SYST_WEIGHT", "f8"),
    ]
    source_path = os.path.join(outdir, f"lsst_y{year}_source_glass_mock_{mock_id}.h5")
    lens_path = os.path.join(outdir, f"lsst_y{year}_lens_glass_mock_{mock_id}.h5")
    h5_source = make_h5_catalogue(source_path, source_fields)
    h5_lens = make_h5_catalogue(lens_path, lens_fields)

    n_source_total = 0
    n_lens_total = 0
    n_source_dropped = 0
    n_lens_dropped = 0

    with tqdm(total=len(shells), desc="Simulating shells", unit="shell") as pbar:
        for i, delta_i in enumerate(matter):
            convergence.add_window(delta_i, shells[i])
            kappa_i = convergence.kappa
            gamm1_i, gamm2_i = glass.shear_from_convergence(kappa_i)

            for gal_lon, gal_lat, gal_count in glass.positions_from_delta(
                src_ngal[i], delta_i, bias, mask_ud, rng=rng,
            ):
                gal_z = glass.redshifts(gal_count, shells[i], rng=rng)
                gal_phz = glass.gaussian_phz(gal_z, src.sigma_z0, rng=rng)
                gal_zbin = np.digitize(gal_phz, src.bin_edges) - 1

                # drop galaxies whose (photo-)z scattered outside the tomographic range —
                # they don't belong to any defined bin, same as a real survey would exclude them
                in_range = (gal_zbin >= 0) & (gal_zbin < src.n_bins)
                n_source_dropped += gal_count - int(in_range.sum())
                if not np.any(in_range):
                    continue
                gal_lon, gal_lat = gal_lon[in_range], gal_lat[in_range]
                gal_z, gal_zbin = gal_z[in_range], gal_zbin[in_range]
                gal_count = int(in_range.sum())

                gal_eps = glass.ellipticity_intnorm(gal_count, src.sigma_e, rng=rng, xp=np)
                gal_she = glass.galaxy_shear(gal_lon, gal_lat, gal_eps, kappa_i, gamm1_i, gamm2_i)
                gal_pix = hp.ang2pix(nside, gal_lon, gal_lat, lonlat=True)
                gal_weight = 1.0 / mask_ud[gal_pix]

                append_h5(
                    h5_source,
                    RA=gal_lon, DEC=gal_lat, PHOTO_Z=gal_z, ZBIN=gal_zbin,
                    G1=gal_she.real, G2=gal_she.imag, SYST_WEIGHT=gal_weight,
                )
                n_source_total += gal_count

            for gal_lon, gal_lat, gal_count in glass.positions_from_delta(
                lens_ngal[i], delta_i, bias, mask_ud, rng=rng,
            ):
                gal_z = glass.redshifts(gal_count, shells[i], rng=rng)
                gal_phz = glass.gaussian_phz(gal_z, lens.sigma_z0, rng=rng)
                gal_zbin = np.digitize(gal_phz, lens.bin_edges) - 1

                in_range = (gal_zbin >= 0) & (gal_zbin < lens.n_bins)
                n_lens_dropped += gal_count - int(in_range.sum())
                if not np.any(in_range):
                    continue
                gal_lon, gal_lat = gal_lon[in_range], gal_lat[in_range]
                gal_z, gal_zbin = gal_z[in_range], gal_zbin[in_range]
                gal_count = int(in_range.sum())

                gal_pix = hp.ang2pix(nside, gal_lon, gal_lat, lonlat=True)
                gal_weight = 1.0 / mask_ud[gal_pix]

                append_h5(
                    h5_lens,
                    RA=gal_lon, DEC=gal_lat, PHOTO_Z=gal_z, ZBIN=gal_zbin,
                    SYST_WEIGHT=gal_weight,
                )
                n_lens_total += gal_count

            pbar.set_postfix(source=n_source_total, lens=n_lens_total)
            pbar.update(1)

    h5_source.close()
    h5_lens.close()
    log.info("wrote %s source galaxies to %s (%s dropped: outside tomographic z-range)",
             f"{n_source_total:,}", source_path, f"{n_source_dropped:,}")
    log.info("wrote %s lens galaxies to %s (%s dropped: outside tomographic z-range)",
             f"{n_lens_total:,}", lens_path, f"{n_lens_dropped:,}")
    return source_path, lens_path, n_source_total, n_lens_total


def save_catalogue_diagnostics(source_path, lens_path, src: Tomography, lens: Tomography, outdir):
    import h5py
    diag_dir = os.path.join(outdir, "diagnostics")
    os.makedirs(diag_dir, exist_ok=True)

    with h5py.File(source_path, "r") as f:
        src_zbin = f["ZBIN"][:]
        src_photoz = f["PHOTO_Z"][:]
    with h5py.File(lens_path, "r") as f:
        lens_zbin = f["ZBIN"][:]
        lens_photoz = f["PHOTO_Z"][:]

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    axes[0].set_title("source catalogue - PHOTO_Z by ZBIN")
    for i in range(src.n_bins):
        axes[0].hist(src_photoz[src_zbin == i], density=True, histtype="step", label=f"bin {i}")
    axes[0].set_xlabel("PHOTO_Z")
    axes[0].legend()

    axes[1].set_title("lens catalogue - PHOTO_Z by ZBIN")
    for i in range(lens.n_bins):
        axes[1].hist(lens_photoz[lens_zbin == i], density=True, histtype="step", label=f"bin {i}")
    axes[1].set_xlabel("PHOTO_Z")
    axes[1].legend()
    plt.tight_layout()
    fig.savefig(os.path.join(diag_dir, "catalogue_zbin_histograms.png"), dpi=120)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    axes[0].set_title("source redshifts in catalogue")
    for i in range(int(src_zbin.max()) + 1):
        axes[0].hist(src_photoz[src_zbin == i], histtype="stepfilled", edgecolor="none",
                     alpha=0.5, bins=50, density=1, label=f"cat. bin {i}")
    for i in range(src.n_bins):
        axes[0].plot(src.z, src.tomo_nz[i] / src.n_gal_arcmin2 * src.n_bins, alpha=0.5,
                     label=f"inp. bin {i}")
    axes[0].plot(src.z, src.nz / src.n_gal_arcmin2 * src.n_bins, ls="--", c="k")
    axes[0].legend(ncol=2, fontsize=7)

    axes[1].set_title("lens redshifts in catalogue")
    for i in range(int(lens_zbin.max()) + 1):
        axes[1].hist(lens_photoz[lens_zbin == i], histtype="stepfilled", edgecolor="none",
                     alpha=0.5, bins=50, density=1, label=f"cat. bin {i}")
    for i in range(lens.n_bins):
        axes[1].plot(lens.z, lens.tomo_nz[i] / lens.n_gal_arcmin2 * lens.n_bins, alpha=0.5,
                     label=f"inp. bin {i}")
    axes[1].plot(lens.z, lens.nz / lens.n_gal_arcmin2 * lens.n_bins, ls="--", c="k")
    axes[1].legend(ncol=2, fontsize=7)
    plt.tight_layout()
    fig.savefig(os.path.join(diag_dir, "catalogue_nz_comparison.png"), dpi=120)
    plt.close(fig)

    log.info("source ZBIN max: %d, lens ZBIN max: %d", int(src_zbin.max()), int(lens_zbin.max()))


# ────────────────────────────── CLI ──────────────────────────────

def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    sim = p.add_argument_group("simulation")
    sim.add_argument(
        "--nside", type=int, required=True,
        help="HEALPix resolution for the lognormal density/shear maps. Does NOT change the "
             "total simulated galaxy count (that's density x footprint area) — only how "
             "finely clustering is resolved on the sky, and therefore runtime.",
    )
    sim.add_argument(
        "--lmax", type=int, default=None,
        help="Maximum angular multipole for the matter power spectra / lognormal fields. "
             "Default: 2*nside.",
    )
    sim.add_argument(
        "--year", choices=["1", "10"], required=True,
        help="LSST survey year preset (Y1 or Y10), resolved via binny's 'lsst' survey "
             "config — sets n_gal_arcmin2, number of tomographic bins, bin edges, and "
             "photo-z scatter for both the lens and source samples.",
    )
    sim.add_argument(
        "--seed", type=int, default=42,
        help="Seed for the numpy random generator driving the whole simulation "
             "(matter fields, positions, redshifts, ellipticities).",
    )

    cosmo_grp = p.add_argument_group(
        "cosmology",
        "Default (no flags at all): Planck 2018 TT,TE,EE+lowE+lensing "
        f"(H0={PLANCK18_H0}, omch2={PLANCK18_OMCH2}, ombh2={PLANCK18_OMBH2}). "
        "--cosmo-config/--camb-ini replace the whole parameter set (or layer over the "
        "Planck18 default, for --cosmo-config); --H0/--omch2/--ombh2 then override "
        "individual values on top of whichever of those was used, last.",
    )
    cosmo_src = cosmo_grp.add_mutually_exclusive_group()
    cosmo_src.add_argument(
        "--cosmo-config", default=None,
        help="YAML file of camb.set_params kwargs (e.g. ns, As, mnu, nnu, w, omk, ...). "
             "Keys present in the file override the Planck 2018 default; keys not "
             "mentioned keep their Planck18 value. Mutually exclusive with --camb-ini.",
    )
    cosmo_src.add_argument(
        "--camb-ini", default=None,
        help="A full CAMB .ini parameter file, passed straight to camb.read_ini() — full "
             "control of the CAMB parameter space for power users, bypassing the Planck18 "
             "defaults entirely. Mutually exclusive with --cosmo-config.",
    )
    cosmo_grp.add_argument(
        "--H0", type=float, default=None,
        help="Hubble constant [km/s/Mpc]. Overrides whatever --cosmo-config/--camb-ini/"
             "Planck18-default would otherwise set. Default: no override.",
    )
    cosmo_grp.add_argument(
        "--omch2", type=float, default=None,
        help="Physical cold dark matter density, Omega_c * h^2. Overrides whatever "
             "--cosmo-config/--camb-ini/Planck18-default would otherwise set. "
             "Default: no override.",
    )
    cosmo_grp.add_argument(
        "--ombh2", type=float, default=None,
        help="Physical baryon density, Omega_b * h^2. Overrides whatever "
             "--cosmo-config/--camb-ini/Planck18-default would otherwise set. "
             "Default: no override.",
    )

    out_grp = p.add_argument_group("mock identity / output")
    out_grp.add_argument(
        "--outdir", required=True,
        help="Output directory (created if missing). Holds the catalogues, run.log, the "
             "cache/ dir (disk-cached CAMB results), and (if enabled) theory Cls / "
             "diagnostics/.",
    )
    out_grp.add_argument(
        "--mock-id", default=None,
        help="Identifier embedded in output filenames: "
             "lsst_y{year}_{lens,source}_glass_mock_{MOCK_ID}.h5 (and the theory SACC file, "
             "if enabled). Default: an auto-generated random 8-character hex id (logged at "
             "the start of the run, so an auto-generated run is still traceable).",
    )
    out_grp.add_argument(
        "--mask-path", default=DEFAULT_MASK_PATH,
        help="Path to a HEALPix visibility map (values 0-1, fraction of each pixel "
             "observed) used both to restrict galaxy sampling to the survey footprint and "
             "to set each galaxy's SYST_WEIGHT = 1/visibility. Default: the bundled LSST "
             "Y10 mask (data/lsst_y10_final_map4096.fits.gz, gzip-compressed but readable "
             "directly — healpy/astropy decompress transparently), reused for both Y1 and "
             "Y10 (no separate Y1 mask currently exists).",
    )

    toggles = p.add_argument_group("optional outputs / dev-loop controls")
    toggles.add_argument(
        "--save-theory-cls", type=int, default=0, choices=[0, 1],
        help="If 1, additionally compute the theory 3x2pt Cls (gg/gG/GG) and save them as "
             "a single SACC file (HDF5-backed), including the tomographic n(z)s actually "
             "used in this run. This costs one extra CAMB reference call. Default: 0 (off).",
    )
    toggles.add_argument(
        "--save-diagnostics", type=int, default=0, choices=[0, 1],
        help="If 1, save diagnostic PNGs under <outdir>/diagnostics/ (mask mollview, "
             "GLASS-vs-CAMB theory ratio [only if --save-theory-cls 1], catalogue "
             "PHOTO_Z/ZBIN histograms and n(z) sanity-check plots). These look the same "
             "run to run, so they're off by default.",
    )
    toggles.add_argument(
        "--debug-density-scale", type=float, default=1.0,
        help="Scale factor applied to both lens/source n_gal_arcmin2. Total galaxy count "
             "is density x footprint area (NOT affected by --nside), so this is how you "
             "shrink a run. Default is 1.0 (the real survey density, a physically valid "
             "mock — a Y1 run is then O(10^8-10^9) rows per catalogue). Pass a smaller "
             "value (e.g. 0.01) explicitly for a fast, NON-physical dev-loop smoke test; "
             "anything other than 1.0 logs a warning that the output is not a physical "
             "mock.",
    )
    toggles.add_argument(
        "--log-level", default="INFO",
        help="Python logging level (DEBUG/INFO/WARNING/ERROR). Logged to both the console "
             "and <outdir>/run.log.",
    )

    args = p.parse_args(argv)
    if args.lmax is None:
        args.lmax = 2 * args.nside
    return args


# module-level handles used by compute_theory (set in main); kept simple since this
# is a one-shot CLI script, not a library.
shells_global = None
cosmo_global = None


def main(argv=None):
    args = parse_args(argv)
    setup_logging(args.outdir, args.log_level)
    cache_dir = os.path.join(args.outdir, "cache")

    mock_id = args.mock_id or secrets.token_hex(4)

    t0 = time.time()
    log.info("starting glass_lsst_mock: mock_id=%s nside=%d lmax=%d year=Y%s outdir=%s",
              mock_id, args.nside, args.lmax, args.year, args.outdir)

    rng = np.random.default_rng(seed=args.seed)

    pars, cosmo = resolve_cosmology(args)
    log.info("cosmology resolved: H0=%.4f omch2=%.5f ombh2=%.5f", pars.H0, pars.omch2, pars.ombh2)
    cosmo_key = cosmology_cache_key(pars, args)

    lsst_cfg = binny.load_survey_config("lsst")
    zmax = float(lsst_cfg["z_grid"]["stop"])
    zb, shells = build_shells(cosmo, zmax)
    global shells_global, cosmo_global
    shells_global = shells
    cosmo_global = cosmo

    cls = get_matter_cls(pars, args.lmax, shells, cache_dir, cosmo_key, zb)

    if args.debug_density_scale != 1.0:
        log.warning(
            "--debug-density-scale=%g set: output is NOT a physical mock, dev-loop testing only",
            args.debug_density_scale,
        )
    log.info("building lens/source tomography via binny (lsst Y%s)", args.year)
    lens = Tomography("lens", args.year, density_scale=args.debug_density_scale)
    src = Tomography("source", args.year, density_scale=args.debug_density_scale)

    if args.save_theory_cls:
        log.info("computing theory 3x2pt Cls")
        compute_theory(pars, args.lmax, cls, lens, src, BIAS, cache_dir, args.outdir,
                        args.year, mock_id, args.save_diagnostics, cosmo_key)

    cls_discretized = glass.discretized_cls(cls, nside=args.nside, lmax=args.lmax, ncorr=3)

    mask, mask_ud = load_mask(args.mask_path, args.nside, args.outdir, args.save_diagnostics)

    log.info("running simulation")
    source_path, lens_path, n_src, n_lens = run_simulation(
        cls_discretized, shells, cosmo, src, lens, BIAS, mask_ud, args.nside, rng,
        args.outdir, args.year, mock_id,
    )

    if args.save_diagnostics:
        save_catalogue_diagnostics(source_path, lens_path, src, lens, args.outdir)

    log.info("done in %.1f s. source=%s lens=%s", time.time() - t0, source_path, lens_path)


if __name__ == "__main__":
    main()
