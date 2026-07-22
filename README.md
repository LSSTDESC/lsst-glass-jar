# LSST GLASS JAR (Just Another Realisation)

This is a self contained repository with all you need to generate LSST Y1/Y10 lens+source simulations using [Glass](https://github.com/glass-dev/glass)

## Setup
After cloning this repository, create a new conda environment using: `conda env create -f environment.yaml`. This will download and install everything you need to create these simulations. Don't forget to activate the conda env by `source activate lsst-glass-jar` or `conda activate lsst-glass-jar`

## Controllable variables:
You can check which variables can be controlled via command line by `python generate_lsst_glass_mock.py -h`. In any case:
* `--nside`: HEALPix resolution for the density/shear maps. Does NOT change the total galaxy count (that's density x footprint area), only clustering resolution and runtime
* `--lmax`: maximum angular multipole [default: 2*nside]
* `--year`: LSST survey year preset, `1` or `10` — resolves n(z)/bins/density/photo-z scatter via [binny](https://binny-org.github.io/binny/)'s `lsst` config
* `--seed`: random seed
* `--outdir`: path to save the mocks
* `--mock-id`: identifier embedded in output filenames [default: auto-generated random hash]
* `--mask-path`: HEALPix visibility map to use [default: the bundled LSST Y10 mask in `masks/`]
* `--cosmo-config`: YAML file of `camb.set_params` kwargs, layered over the Planck 2018 default
* `--camb-ini`: a full CAMB `.ini` parameter file (mutually exclusive with `--cosmo-config`)
* `--H0` / `--omch2` / `--ombh2`: individual cosmology overrides, applied last on top of whichever of the above was used
* `--save-theory-cls`: `0`/`1`, compute+save the theory 3x2pt Cls as a SACC file [default: `0`]
* `--save-diagnostics`: `0`/`1`, save mask/theory/catalogue diagnostic PNGs [default: `0`]
* `--debug-density-scale`: scale factor on the survey density, for shrinking a run [default: `1.0`, the real physical density]
* `--log-level`: Python logging level [default: `INFO`]

## Testing
To smoke-test the pipeline without waiting on a full physical mock (a real run is O(10^8-10^9) rows per catalogue, set by survey density x footprint area — not by `--nside`), shrink the density directly:

    python generate_lsst_glass_mock.py --nside 64 --year 1 --outdir out --debug-density-scale 0.01

This runs the full pipeline (cosmology, tomography, mask, lognormal sim, HDF5 catalogues) at 1% of the real LSST Y1 density, in seconds rather than minutes.

## Simulation Parameters
The following is independent of `--nside`, `--year`, or the survey mask:

| Parameter | Value |
|---|---|
| Cosmology | Planck 2018 TT,TE,EE+lowE+lensing: H0=67.36, omch2=0.1200, ombh2=0.02237 |
| redshift range | (0.0, 3.5) — pulled from [binny's `lsst` z-grid](https://binny-org.github.io/binny/main/examples/surveys/survey_presets.html), not hardcoded |
| number of shells | 25 (comoving-distance, linear windows) |
| lmax | = 2 * nside (default) |
| galaxy bias | constant = 1.5 |
| redshift-distr | Smail et al. distr. (via binny's `lsst` preset, per role/year) |

The following depend on `--year` and role (lens = clustering, source = shear), resolved via [binny's `lsst` survey config](https://binny-org.github.io/binny/main/examples/surveys/survey_presets.html) rather than hardcoded per-population like a single-sample pipeline:

| | lens Y1 | source Y1 | lens Y10 | source Y10 |
|---|---|---|---|---|
| number density [gal/arcmin²] | 18 | 10 | 48 | 27 |
| tomographic bins | 5 | 5 | 10 | 5 |
| photo-z error (z=0), sigma_z0 | 0.03 | 0.05 | 0.03 | 0.05 |
| shape noise, sigma_e | — | 0.26 | — | 0.26 |
