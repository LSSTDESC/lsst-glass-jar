#!/usr/bin/env python
"""glass_lsst_mock.py — generate an LSST Y1/Y10 lens+source GLASS mock catalogue.

Self-contained CLI script: cosmology (CAMB) -> matter Cls (disk-cached) -> lens/source
n(z) tomography via binny's LSST survey preset -> theory 3x2pt Cls (gg/gG/GG) -> a
lognormal GLASS simulation streamed straight to two HDF5 catalogues (source, lens),
with a survey mask applied both to galaxy sampling and to a per-galaxy SYST_WEIGHT
(1/visibility). Diagnostic plots are written as PNGs (no interactive display).

Example:
    python glass_lsst_mock.py --nside 1024 --year 1 --outdir /path/to/output
"""
from __future__ import annotations

import argparse
import hashlib
import logging
import os
import pickle
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
DEFAULT_MASK_PATH = "/home/loureiro/Projects/lsst_y10_final_map4096.fits"

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


# ────────────────────────────── cosmology / shells ──────────────────────────────

def build_cosmology(H0: float, omch2: float, ombh2: float):
    pars = camb.set_params(
        H0=H0,
        omch2=omch2,
        ombh2=ombh2,
        NonLinear=camb.model.NonLinear_both,
    )
    results = camb.get_background(pars)
    cosmo = Cosmology(results)
    return pars, cosmo


def build_shells(cosmo, zmax: float, num: int = N_SHELLS):
    zb = glass.distance_grid(cosmo, 0.0, zmax, num=num)
    shells = glass.linear_windows(zb)
    return zb, shells


def get_matter_cls(pars, lmax, shells, cache_dir, H0, omch2, ombh2, zb):
    return cached(
        cache_dir,
        "matter_cls",
        (H0, omch2, ombh2, lmax, zb.tolist()),
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


# ────────────────────────────── theory 3x2pt Cls ──────────────────────────────

def compute_theory(pars, lmax, cls, lens: Tomography, src: Tomography, bias, cache_dir, outdir):
    diag_dir = os.path.join(outdir, "diagnostics")
    os.makedirs(diag_dir, exist_ok=True)

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
            pars.H0, pars.omch2, pars.ombh2, lmax, lens.n_bins, src.n_bins,
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

    save_gg = bias**2 * glass_pos
    save_gG = bias * glass_cross
    save_GG = glass_she

    for i in range(lens.n_bins):
        for j in range(i, lens.n_bins):
            np.savetxt(
                os.path.join(outdir, f"glass_theory_cls_gg_bin{i}_{j}.txt"),
                np.c_[ell[1:], save_gg[i, j, 1:]],
                header=f"GLASS clustering (gg) Cls, lens bins {i}x{j}, bias={bias}\nell, C_ell",
            )
    for i in range(src.n_bins):
        for j in range(i, src.n_bins):
            np.savetxt(
                os.path.join(outdir, f"glass_theory_cls_GG_bin{i}_{j}.txt"),
                np.c_[ell[1:], save_GG[i, j, 1:]],
                header=f"GLASS cosmic shear (GG) Cls, source bins {i}x{j}\nell, C_ell",
            )
    for i in range(lens.n_bins):
        for j in range(src.n_bins):
            np.savetxt(
                os.path.join(outdir, f"glass_theory_cls_gG_lens{i}_src{j}.txt"),
                np.c_[ell[1:], save_gG[i, j, 1:]],
                header=f"GLASS galaxy-galaxy lensing (gG) Cls, lens bin {i} x source bin {j}, "
                       f"bias={bias}\nell, C_ell",
            )

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

def load_mask(mask_path: str, nside: int, outdir: str):
    diag_dir = os.path.join(outdir, "diagnostics")
    os.makedirs(diag_dir, exist_ok=True)

    y_mask = 1.0 - hp.read_map(mask_path)
    y_ud = hp.ud_grade(y_mask, nside)

    mean_vis = float(np.mean(1.0 / y_mask[y_mask != 0]))
    log.info("mask loaded from %s: mean(1/visibility) over unmasked pixels = %.4f", mask_path, mean_vis)

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
                    nside, rng, outdir, year):
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
    source_path = os.path.join(outdir, f"glass_lsst_y{year}_source_catalogue.h5")
    lens_path = os.path.join(outdir, f"glass_lsst_y{year}_lens_catalogue.h5")
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
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--nside", type=int, required=True)
    p.add_argument("--lmax", type=int, default=None, help="default: 2*nside")
    p.add_argument("--year", choices=["1", "10"], required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--mask-path", default=DEFAULT_MASK_PATH)
    p.add_argument("--H0", type=float, default=PLANCK18_H0)
    p.add_argument("--omch2", type=float, default=PLANCK18_OMCH2)
    p.add_argument("--ombh2", type=float, default=PLANCK18_OMBH2)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--log-level", default="INFO")
    p.add_argument(
        "--debug-density-scale", type=float, default=1.0,
        help=(
            "scale factor on lens/source n_gal_arcmin2, for a fast dev-loop smoke test "
            "(e.g. 0.01 cuts total galaxy count ~100x). Total galaxy count is density x "
            "footprint area, NOT affected by --nside, so this is how you shrink a test run. "
            "Default 1.0 = full survey density (a physically valid mock). Anything else is "
            "NOT a physical mock."
        ),
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

    t0 = time.time()
    log.info("starting glass_lsst_mock: nside=%d lmax=%d year=Y%s outdir=%s",
              args.nside, args.lmax, args.year, args.outdir)

    rng = np.random.default_rng(seed=args.seed)

    log.info("building cosmology: H0=%.4f omch2=%.5f ombh2=%.5f", args.H0, args.omch2, args.ombh2)
    pars, cosmo = build_cosmology(args.H0, args.omch2, args.ombh2)

    lsst_cfg = binny.load_survey_config("lsst")
    zmax = float(lsst_cfg["z_grid"]["stop"])
    zb, shells = build_shells(cosmo, zmax)
    global shells_global, cosmo_global
    shells_global = shells
    cosmo_global = cosmo

    cls = get_matter_cls(pars, args.lmax, shells, cache_dir, args.H0, args.omch2, args.ombh2, zb)

    if args.debug_density_scale != 1.0:
        log.warning(
            "--debug-density-scale=%g set: output is NOT a physical mock, dev-loop testing only",
            args.debug_density_scale,
        )
    log.info("building lens/source tomography via binny (lsst Y%s)", args.year)
    lens = Tomography("lens", args.year, density_scale=args.debug_density_scale)
    src = Tomography("source", args.year, density_scale=args.debug_density_scale)

    log.info("computing theory 3x2pt Cls")
    compute_theory(pars, args.lmax, cls, lens, src, BIAS, cache_dir, args.outdir)

    cls_discretized = glass.discretized_cls(cls, nside=args.nside, lmax=args.lmax, ncorr=3)

    mask, mask_ud = load_mask(args.mask_path, args.nside, args.outdir)

    log.info("running simulation")
    source_path, lens_path, n_src, n_lens = run_simulation(
        cls_discretized, shells, cosmo, src, lens, BIAS, mask_ud, args.nside, rng,
        args.outdir, args.year,
    )

    save_catalogue_diagnostics(source_path, lens_path, src, lens, args.outdir)

    log.info("done in %.1f s. source=%s lens=%s", time.time() - t0, source_path, lens_path)


if __name__ == "__main__":
    main()
