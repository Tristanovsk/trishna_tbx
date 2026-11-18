# Changelog

All notable changes to trishna_tbx are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and the project follows [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Read the Docs documentation (`docs/source`, `.readthedocs.yaml`): user guide (installation, usage,
  methods, pre-launch assessment, processing chain), the VRE gallery notebooks as tutorials, API
  reference, roadmap and changelog. `docs` extra in `pyproject.toml`.
- Docstrings of the public API in reST fields (`:param:`, `:return:`), with the equations.
- Module `validation`: HYPERNETS reader and quality control, matchups with Sentinel-2 L1C subsets,
  GRS runs at 60 m with the TRISHNA-like and all S2 bands (cached), look-up of the archived GRS L2A,
  validation statistics.
- `emulate.l1c_dataset_as_trishna`: L1C subsets (NetCDF of GRSdriver) to TRISHNA-like bands at 60 m.
- VRE notebook *Validation with HYPERNETS: Sentinel-2 vs TRISHNA bands* (`notebooks/5_validation`).

### Changed
- `docs/roadmap.md` moved to `ROADMAP.md`.
- `notebooks/build_notebooks.py`: no timeout of the notebook cells (first GRS runs of the
  validation notebook).

## [0.0.1] - 2026-11

### Added
- Sensor models of TRISHNA VSWIR (super-Gaussian responses, SNR) and Sentinel-2 MSI (ESA responses),
  noise model (`sensor`).
- Emulation of TRISHNA from Sentinel-2: band restriction, 60 m, spectral band adjustment factors,
  OWT-dependent adjustment, noise, GRS band set (`emulate`).
- Optical water types with SRF-convolved references (Bi et al. 2024, Spyrakos et al. 2018), spectral
  angle, merged TRISHNA class set and Monte-Carlo benchmark (`owt`).
- Water-quality algorithm registry with the wavelengths each algorithm needs, coverage table and
  strict runs (`algorithms`).
- Two-band ratio water vapour retrieval, TRISHNA 860/910 nm vs S2 865/945 nm (`water_vapor`).
- radcalnet_oc wrapper for the TOA simulations, AERONET-OC inputs, sensitivity grids (`simulate`).
- invRrs inversion on the TRISHNA bands and information content (`inversion`).
- Readers of GRS L2A and CAMS (`io`).
- Scripts `prelaunch_assessment.py` and `l2a_s2_vs_trishna.py`; 9 VRE notebooks.
