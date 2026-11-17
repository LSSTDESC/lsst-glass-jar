<p align="center">
  <img src="lsst_in_a_jar.png" width="409.6" height="338"
 alt="LSST Glass Jar"/>
</p>
<div align="center">

# LSST GLASS JAR (Just Another Realisation)

</div>

This is a self contained repository with all you need to generate LSST Y1/Y10 lens+source simulations using [GLASS](https://github.com/glass-dev/glass)

## Setup
After cloning this repository, create a new conda environment using: 

    conda env create -f environment.yaml 

This will download and install everything you need to create these simulations. 
Don't forget to activate the conda env with `conda activate lsst-glass-jar`.

## Controllable variables:
You can check which variables can be controlled via command line by 
    python generate_lsst_glass_mock.py -h 
    
In any case:
* `--nside`: HEALPix resolution for the density/shear maps. Does NOT change the total galaxy count (that's density x footprint area), only clustering resolution and runtime [default: `2048`; the 4096 mask is downgraded. Higher nside lowers the lens clustering, see [Known limitations](#known-limitations)]
* `--shell-dx`: comoving width of the matter shells in Mpc [default: `100`, ~69 shells to z=3.5]
* `--exact-cls`: flag, compute the shell matter Cls exactly at every ell [default: off, Limber above ell=100]
* `--lmax`: maximum angular multipole [default: `nside`]
* `--year`: LSST survey year preset, `1` or `10` — resolves n(z)/bins/density/photo-z scatter via [binny](https://binny-org.github.io/binny/)'s `lsst` config
* `--seed`: random seed
* `--outdir`: path to save the mocks
* `--mock-id`: identifier embedded in output filenames [default: auto-generated random hash]
* `--mask-path`: HEALPix visibility map to use [default: the bundled LSST Y10 mask in `masks/`]
* `--mask-cut`: binarize the mask at this visibility threshold instead of using it as a continuous weight — see [Visibility Cut](#visibility-cut) below [default: `None`]
* `--cosmo-config`: YAML file of `camb.set_params` kwargs, layered over the Planck 2018 default
* `--camb-ini`: a full CAMB `.ini` parameter file (mutually exclusive with `--cosmo-config`)
* `--H0` / `--omch2` / `--ombh2`: individual cosmology overrides, applied last on top of whichever of the above was used
* `--save-theory-cls`: flag, compute+save the theory 3x2pt Cls as a SACC file [default: off]
* `--save-diagnostics`: flag, save mask/theory/catalogue diagnostic PNGs [default: off]
* `--debug-density-scale`: scale factor on the survey density, for shrinking a run [default: `1.0`, the real physical density]
* `--log-level`: Python logging level [default: `INFO`]

## Testing
To smoke-test the pipeline without waiting on a full physical mock (a real run is O(10^8-10^9) rows per catalogue, set by survey density x footprint area — not by `--nside`), shrink the density directly:

    python generate_lsst_glass_mock.py --nside 64 --year 1 --outdir out --debug-density-scale 0.01

This runs the full pipeline (cosmology, tomography, mask, lognormal sim, HDF5 catalogues) at 1% of the real LSST Y1 density, in seconds rather than minutes.

## Simulation Parameters

> [!IMPORTANT]
> **No intrinsic alignments.** These catalogues have NO intrinsic alignment (IA) model applied. Source
> galaxy shapes are pure shape-noise + gravitational shear (`G1`/`G2`) — if your analysis needs an IA
> signal, it isn't in here.

The following is independent of `--nside`, `--year`, or the survey mask:

| Parameter | Value |
|---|---|
| Cosmology | Planck 2018 TT,TE,EE+lowE+lensing: H0=67.36, omch2=0.1200, ombh2=0.02237 |
| redshift range | (0.0, 3.5) — pulled from [binny's `lsst` z-grid](https://binny-org.github.io/binny/main/examples/surveys/survey_presets.html), not hardcoded |
| matter shells | linear windows every 100 Mpc comoving (`--shell-dx`), ~69 to z=3.5 |
| lmax | = nside (default) |
| galaxy bias | constant = 1.5 |
| redshift-distr | Smail et al. distr. (via binny's `lsst` preset, per role/year) |

The following depend on `--year` and role (lens = clustering, source = shear), resolved via [binny's `lsst` survey config](https://binny-org.github.io/binny/main/examples/surveys/survey_presets.html) rather than hardcoded per-population like a single-sample pipeline:

| | lens Y1 | source Y1 | lens Y10 | source Y10 |
|---|---|---|---|---|
| number density [gal/arcmin²] | 18 | 10 | 48 | 27 |
| tomographic bins | 5 | 5 | 10 | 5 |
| photo-z error (z=0), sigma_z0 | 0.03 | 0.05 | 0.03 | 0.05 |
| shape noise, sigma_e | — | 0.26 | — | 0.26 |

## Known limitations

**Lens clustering sits below the theory Cls.** Two effects of the lognormal + linear-bias model, both in
the catalogue itself, not in any analysis code:

* **Shell width.** A galaxy clusters with the whole matter shell it is drawn in, so a shell wider than a
  lens tomographic bin dilutes that bin's clustering. The old fixed 25-shell grid (~250 Mpc) lost 13-20%
  of the lens auto-spectra; 100 Mpc shells lose 2-4%.
* **Clipping of the linear bias.** Galaxies are sampled from `1 + b*delta` with `b = 1.5`, set to zero
  where negative. The clipped fraction grows with the pixel variance of a shell: thinner shells and
  higher `--nside` clip more, mostly at low redshift.

Expected lens C_gg relative to continuous theory (CAMB/CCL), ell 30-300, Y1:

| `--shell-dx`, `--nside` | lens0 | lens1 | lens2 | lens3 | lens4 |
|---|---|---|---|---|---|
| old default (25 shells, 4096) | 0.82 | 0.84 | 0.83 | 0.81 | 0.79 |
| **100 Mpc, 2048 (default)** | **0.88** | **0.93** | **0.96** | **0.96** | **0.96** |
| 100 Mpc, 4096 | 0.82 | 0.88 | 0.93 | 0.95 | 0.96 |
| measured: one 100 Mpc, 2048 mock (±0.01-0.03) | 0.85 | 0.96 | 0.97 | 0.96 | 0.91 |

The rows without "measured" are predictions (shell weights x clipping model). The model matched full-sky
GLASS runs at nside 512 to 0.01 and the old 4096 mock to within 0.03.

Galaxy-galaxy lensing is 0.95-1.00 of theory at the defaults, and cosmic shear 0.98-1.01. Fitting these
catalogues with a linear-bias model will therefore pull the lens biases low (most for the lowest bins).
The theory SACC written by `--save-theory-cls` is GLASS's own prediction for the shells, without the
clipping; it does not include these losses either.

Very thin shells (well below 100 Mpc) at high `lmax` can make GLASS stop with "covariance matrix is not
positive definite" for shells at z < 0.25.

**Shear convention.** `G1`/`G2` are the components in the local (east, north) basis: positive `G1` is
elongation along RA. HEALPix (and NaMaster) use the (theta, phi) basis, which is rotated by 90 degrees:
there, `Q + iU = -(G1 + i G2)`. Feeding `G1`/`G2` as `Q`/`U` unchanged flips the sign of
galaxy-galaxy lensing (cosmic shear is unaffected).

## Visibility Cut

The bundled mask has values in `[0, 1]` — the probability of a galaxy being observed in that pixel. `--mask-cut` controls how that's used:

* **No cut given (`--mask-cut` unset, default)**: the continuous `[0, 1]` map is used as-is, both to weight galaxy sampling and as `1/visibility` in the catalogue's `SYST_WEIGHT` column.
* **Cut given (e.g. `--mask-cut 0.5`)**: after downgrading the mask to the simulation's `--nside`, pixels with value `>= cut` are set to `1` and pixels `< cut` are set to `0` — a binary mask. `SYST_WEIGHT` is then always `1.0` for every galaxy, since a binary mask needs no incompleteness correction.

The mask actually used for sampling (continuous or binarized) is saved to `survey_footprint_N{nside}.fits` either way.

## Outputs

Written to `--outdir`:

| File | Contents |
|---|---|
| `lsst_y{year}_source_glass_mock_{MOCK_ID}.h5` | Source catalogue. Columns: `RA`, `DEC`, `PHOTO_Z`, `ZBIN`, `G1`, `G2`, `SYST_WEIGHT` |
| `lsst_y{year}_lens_glass_mock_{MOCK_ID}.h5` | Lens catalogue. Columns: `RA`, `DEC`, `PHOTO_Z`, `ZBIN`, `SYST_WEIGHT` |
| `lsst_y{year}_theory_glass_mock_{MOCK_ID}.sacc.h5` | Theory 3x2pt Cls (SACC, HDF5-backed) — only if `--save-theory-cls`. NZ tracers `lens0..N`/`src0..N` plus `galaxy_density_cl`/`galaxy_shear_cl_ee`/`galaxy_shearDensity_cl_e` data points |
| `survey_footprint_N{nside}.fits` | The mask actually used for sampling, at simulation `nside` (continuous or binarized per `--mask-cut`) |
| `run.log` | Full log of the run |
| `cache/` | Disk-cached CAMB results (not a data product — safe to delete, just costs a recompute) |
| `diagnostics/*.png` | Mask/theory-ratio/catalogue sanity-check plots — only if `--save-diagnostics` |
