# Roadmap and notebook plan

Mapped to MAG-25-PT-121 Task 1 (2.2.1.x) and to the VRE gallery (Task 2).

## A. Integration work in the sibling packages (small patches)

| Package | Issue found | Fix |
|---|---|---|
| grs | `high_swir` flag needs an exact 2190 nm band (mask.py:561) | make it optional, or fall back to 1610 nm for TRISHNA; until then `grs_band_set(strict=False)` keeps B12 for masking only |
| grs | aerosol/glint use the last 2 processed bands | OK with the TRISHNA set (865 + 1610); document it and add a test |
| grs | TCWV comes only from CAMS | option to read an image-based TCWV map (`trishna_tbx.water_vapor`) |
| grs | no Python way to pass a band subset except a pre-built Dataset | expose `band_idx` / `wl_to_process` in `Process.execute` |
| GRSl2bgen | `sel(method='nearest')` silently swaps 443→485 and 705→670 | add a tolerance check (`algorithms.Algo.covered`) |
| GRSl2bgen | OWT references interpolated at band centres; default `wl_range` drops 860 nm | convolve with the SRF (`owt.convolve_reference`); pass `wl_range` through |
| GRSl2bgen | `Chl.process` docstring is stale (M09B/NIRB called without Rrs) | clean up |
| invRrs | sensors are Py6S presets only | accept a custom RSR (`inversion.SensorRSR` works now via injection) |
| pyOWT | needs a sensor entry plus a 5th-order AVW regression | fit the TRISHNA AVW regression from hyperspectral libraries convolved with the TRISHNA SRF |
| radcalnet_oc | output requires `aot865`/`ang_exp` in the input | done in `simulate.build_input`; consider defaults upstream |

## B. Science tasks

1. **Atmospheric correction (2.2.1.1)**
   - GRS on S2 L1C restricted to TRISHNA bands (B2, B3, B4, B8A, B11, plus B10 for cirrus).
   - Compare with full-band GRS on AERONET-OC/HYPERNETS matchups (satmatchup).
   - Sensitivity to the aerosol model, AOT, TCWV and O3 through radcalnet_oc simulations
     (`simulate.sensitivity_grid`).
   - POLYMER is the CNES baseline: compare it too.
2. **TOA simulations (2.2.1.2)**
   - AERONET-OC and HYPERNETS spectra through radcalnet_oc, convolved to TRISHNA and S2.
   - Use the paired values for SBAF (by OWT), noise studies and re-tuning of
     OC2/KD2/Nechad (`algorithms.retune_ocx`).
   - To do: a HYPERNETS reader (WATERHYPERNET portal) and a 910 nm in situ check.
3. **OWT (2.2.1.3)**
   - Replace the synthetic members with the Bi et al. 2024 training spectra and the in situ databases.
   - Define a *TRISHNA OWT set* by merging the collinear sub-types (1/2, 3a/3b, 4a/4b).
   - Add the pyOWT AVW/Area/NDI approach with a TRISHNA regression.
   - Evaluate OWT-conditioned algorithm blending restricted to TRISHNA-compatible algorithms.
4. **Water vapour**
   - Retrieve over land, then interpolate over water. Port the hgrs coupled path-length model
     for water pixels.
   - Assess the continuum error over vegetation and soils (910 nm sits closer to 860 nm than
     945 nm does to 865 nm).
   - Get the actual WV-band SNR.
5. **Thermal synergy**
   - Use TIR (SST) for cloud/water masking, and as a prior to break the chl/CDOM degeneracy
     in time series (bloom detection, upwelling).
6. **Case studies (2.2.1.4)**: ODATIS/THEIA sites, AERONET-OC/HYPERNETS sites, REPHY/COAST-HF.
   Indian Ocean Chl-a internship (contact via ODATIS).

## C. VRE notebook gallery (proposal Task 2)

Written. See `notebooks/README.md` for the gallery cards and data. Generate them with
`notebooks/build_notebooks.py`.

| Gallery section | Notebook |
|---|---|
| 0_basics | `trishna_vs_sentinel2.ipynb` |
| 5_validation | `validation_hypernets_trishna.ipynb` (HYPERNETS BEFR matchups, both band sets) |
| 1_importdata | `emulate_trishna_l1c.ipynb` (complements *Simulate TRISHNA data*) |
| 3_atmosphericcorrection | `grs_trishna_bands.ipynb`, `water_vapor_trishna.ipynb` |
| 4_visualisation | `optical_water_types_trishna.ipynb`, `wq_products_trishna.ipynb`, `datacube_trishna.ipynb` |
| 6_trishna_simulations (new) | `toa_simulation_radcalnet_oc.ipynb` (the RadCalNet-OC tutorial requested at the ODATIS meeting), `inversion_information_content.ipynb` |

Next:
- extend the HYPERNETS validation to other sites (MAFR, VEIT, GAIT...) and to the 2023-2025 data
- replace the Magellium data paths with TREX datalake paths
- run on the ODATIS/THEIA case-study sites
