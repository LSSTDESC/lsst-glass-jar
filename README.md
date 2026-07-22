# LSST GLASS JAR (Just Another Realisation)
A self-contained repository for generating LSST log-normal mock galaxy catalogues with [GLASS](https://glass.readthedocs.io/stable/)

## Install

    conda env create -f environment.yaml
    conda activate lsst-glass-jar

Note: `binny` is imported as `binny` but its PyPI/conda package name is `pybinny` — already
handled by `environment.yaml`, just flagging it in case you install anything by hand.

## Run

    # full, physically valid Y1 mock (default cosmology: Planck 2018)
    python generate_lsst_glass_mock.py --nside 1024 --year 1 --outdir out

    # fast dev-loop smoke test (shrink the galaxy count directly; --nside does not do this)
    python generate_lsst_glass_mock.py --nside 64 --year 1 --outdir out --debug-density-scale 0.01

See `python generate_lsst_glass_mock.py --help` for the full option reference (cosmology
overrides/config files, theory-Cls/diagnostics toggles, mock-id, mask path, ...).
