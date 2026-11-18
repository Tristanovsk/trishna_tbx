# trishna_tbx

Pre-launch toolbox to assess what **TRISHNA VSWIR** can monitor over inland and coastal
waters compared with Sentinel-2. It covers water-quality parameters, Optical Water Types (OWT)
and atmospheric water vapour. Context: CNES study *Support aux développements d'un produit
TRISHNA qualité des eaux* (MAG-25-PT-121), Task 1 (scientific developments) feeding the VRE (Task 2).

TRISHNA has been planned for launch in Q1 2028, with first data about 9 months later
(ODATIS/CNES meeting notes). Until then, the sensor is emulated from S2/Landsat or simulated
with radcalnet_oc.

## TRISHNA VSWIR vs Sentinel-2 (water-relevant bands)

| TRISHNA | centre/FWHM (nm) | SNR @ Lmin | S2 proxy | S2 centre/FWHM | note |
|---|---|---|---|---|---|
| Blue | 485 / 70 | 50 | B2 | 490 (real centroid 496.5) / 66 | TRISHNA is about 11 nm bluer |
| Green | 555 / 70 | 100 | B3 | 560 / 36 | 2x wider |
| Red | 670 / 60 | 100 | B4 | 665 / 31 | 2x wider |
| NIR | 860 / 40 | 100 | B8A | 865 / 21 | |
| WV | 910 / 30 | 50 | (B9) | 945 / 20 | **no S2 equivalent**: different absorption depth |
| Cirrus | 1380 / 30 | 25 | B10 | 1375 / 31 | |
| SWIR | 1610 / 100 | 100 | B11 | 1610 / 91 | aerosol/glint (no 2190 nm) |

TRISHNA has no band at 443 nm, none in the red edge (705/740/783 nm) and none at 2190 nm.
Its pixel size is 60 m (57–91 m across track). Source: Table 6 of the TRISHNA mission paper
(EGUsphere-2026-4492). SRFs are modelled as super-Gaussian (`expon=4`) until ISRO publishes
the measured ones.

## Package layout

```
trishna_tbx/
  sensor.py       TRISHNA/S2 band definitions, SRF (S2 from GRSdriver xlsx), SNR noise model
  emulate.py      S2 L1C/L2A -> TRISHNA-like (band restriction, 60 m, SBAF, noise); GRS band sets
  simulate.py     radcalnet_oc wrapper: Rrs (in situ / OWT / invRrs) -> TOA, convolved to both sensors
  water_vapor.py  two-band ratio TCWV (TRISHNA 860/910 vs S2 865/945), REPTRAN H2O from radcalnet_oc
  owt.py          OWT with SRF-convolved references (Bi2024, Spyrakos2018), SAM, Monte-Carlo benchmark
  algorithms.py   WQ algorithm registry with explicit band needs -> coverage table, strict runs
  inversion.py    invRrs on TRISHNA bands (RSR injection) + information-content analysis
  validation.py   HYPERNETS matchups: S2 vs TRISHNA band sets through GRS, statistics
  io.py           GRS L2A and CAMS readers
scripts/
  prelaunch_assessment.py   reproduces all numbers below (no image needed)
  l2a_s2_vs_trishna.py      real GRS L2A image: full S2 vs S2-as-TRISHNA products and OWT
tests/            pytest (8 tests)
docs/             Read the Docs documentation (Sphinx)
```

What it reuses:
- **GRSdriver**: L1C reading, `band_idx` restriction, S2 SRF.
- **grs**: atmospheric correction, run on the restricted band set (`emulate.grs_band_set`).
- **GRSl2bgen**: OWT reference data and algorithm formulas.
- **radcalnet_oc**: TOA simulation and the gas LUT.
- **invRrs**: forward model and inversion.
- **hgrs**: coupled water-vapour model (to be ported for water pixels).

## First results (pre-launch, `scripts/prelaunch_assessment.py`)

**1. Algorithm coverage.**
- **Run on TRISHNA:**
  - OC2-type blue/green Chl-a (clear to moderate waters only)
  - Nechad SPM (665)
  - Dogliotti turbidity (665/865)
  - red/NIR SPM index
  - KD2 Kd490
  - Roy & Das Kd_PAR
  - the Bi et al. NDI
  - a QAA variant with 485 nm in place of 443 nm, still to develop
- **Lost:**
  - everything using the red edge: NDCI, Gons, Gilerson, M09, Binding 740, Brezonik CDOM
  - the 443 nm algorithms: OC3, QAAv6, Jiang SPM switching
  - the OWT-blended Chl/SPM of GRSl2bgen

The main consequence is that **Chl-a in productive/eutrophic inland waters (red-edge
algorithms) is not possible**, and CDOM loses its standard algorithms.

**2. OWT (Bi 2024, 10 classes).** Values are the percentage of synthetic members given the
same class as their hyperspectral spectrum:

| band set | noise-free | with spec noise |
|---|---|---|
| S2 (8 water bands) | 95 % | 89 % |
| S2 restricted to TRISHNA bands | 85 % | 78 % |
| TRISHNA SRF (485/555/670/860) | 89 % | 78 % |
| TRISHNA without 860 nm | 78 % | – |

- On TRISHNA, the sub-types 3a/3b (2.0°), 4a/4b (3.2°) and 1/2 (5.4°) are nearly collinear.
  **About 7 OWT groups are realistically separable, not 10.**
- The NIR band matters (+11 points). GRSl2bgen's default `wl_range=slice(350, 800)` drops it.
- On a real GRS L2A image (T33PVQ, turbid lake, 25,700 pixels at 60 m), S2-as-TRISHNA agrees
  with S2-full on **82 %** of pixels. Nearly all the disagreement is OWT 6 classified as 5b.

**3. S2 → TRISHNA band adjustment.**
- On average Green, Red and NIR need a correction of 2–6 %.
- **Blue needs ×1.28 (clear) to ×0.84 (turbid)**, and Red up to ×1.17 in very turbid waters
  (OWT 5b), so the factors depend on the OWT. Emulating TRISHNA by simply dropping S2 bands
  is biased; `emulate.owt_band_adjustment` applies OWT-dependent factors.
- On Lake Chad (20 dates), with S2 algorithm coefficients, the adjustment raises SPM/turbidity
  by 7–16 % and Kd490 by 20–40 %. The algorithms have to be re-tuned for TRISHNA.

**4. Water vapour (ratio method, air mass 3).**
- The 910 nm band absorbs less than 945 nm. Its ratio stays usable at high TCWV, with no
  saturation (ratio 0.37 against 0.03 at W = 60 kg m⁻²), which supports the mission choice.
- However, per unit of W it is 3–4× less sensitive. With the requirement SNRs (WV: 50 at
  Lmin = 30 against S2 B9: 114 at 9 W m⁻² sr⁻¹ µm⁻¹), the noise-only TCWV error is 2–4× larger
  for TRISHNA (land, ρ = 0.3: 0.5–3.5 against 0.1–1.6 kg m⁻²).
- Over water (ρ < 0.01) neither sensor retrieves W from the surface signal. Use land pixels
  plus interpolation, or the coupled path-length model of hgrs.
- Results depend on the actual WV-band SNR: request it from CNES/ISRO.

**5. Semi-analytical inversion (invRrs `MultispectralSimple`: chl, a_cdom, spm, offset).**
- With 4 bands and 4 unknowns TRISHNA is only just determined, with corr(chl, a_cdom) of
  −0.8 to −0.99.
- Chl-a relative σ is 2–16 (not retrievable), against 0.7–9 for S2.
- SPM stays retrievable (σ ≈ 20 %) in low/moderate turbidity.
- What remains on TRISHNA: SPM/turbidity, Kd/transparency, coarse OWT groups, and blue/green
  Chl-a in clear waters. Chl-a beyond that needs constraints such as OWT priors, the
  temperature from the TIR bands, or time series.

**6. HYPERNETS validation (Berre lagoon, 64 S2 matchups, 2021–2023; `5_validation` notebook).**
- In situ TRISHNA/S2 band ratios are 0.92 (Blue), 0.94 (Green), 1.00 (Red) and 1.01 (NIR). The
  OWT-based emulation model matches Blue and Red within ~1 % but is off by 3–5 % for Green/NIR.
- OWT on real spectra: TRISHNA bands recover the hyperspectral class in 93 % of cases (S2: 94 %).
  The synthetic benchmark (78 %) is pessimistic.
- GRS with the TRISHNA band set (aerosol/glint from 865/1610 nm) over-corrects: bias about
  −0.002 sr⁻¹ in the visible, MdAPD 35–60 % against 23–40 % for the S2 band set with the same
  settings. The bias grows where 865/1610 retrieves a higher AOT, because the NIR water signal is
  taken as aerosol. TRISHNA atmospheric correction needs a NIR water-signal treatment.
- The archived standard GRS (20 m) does about twice better than the 60 m re-processing for both
  band sets, so the absolute errors here are pessimistic.

Caveats:
- The members are synthetic (mean ± sd shape perturbations). They should be replaced with
  HYPERNETS/AERONET-OC or the Bi et al. training set.
- The SRFs are modelled.
- The noise is white, at spec level, and S2 and TRISHNA SNRs are specified at different radiances.
- The algorithm coefficients have not been re-tuned for the TRISHNA bands.

## Quick start

```bash
pip install -e .            # + sibling packages: GRSdriver grs GRSl2bgen radcalnet_oc invRrs hgrs
python scripts/prelaunch_assessment.py --inversion
python scripts/l2a_s2_vs_trishna.py <S2_L2Agrs.nc> --crop 0:3600,0:960 --out cmp.nc
pytest -q tests
```

```python
from trishna_tbx import TRISHNA, make_s2, algorithms
from trishna_tbx.emulate import open_s2_l1c_as_trishna      # L1C restricted + 60 m via GRSdriver
from trishna_tbx.water_vapor import RatioRetrieval
prod = open_s2_l1c_as_trishna('S2B_MSIL1C_....SAFE')       # feed to grs.Process().execute(prod, ...)
print(algorithms.coverage([make_s2(), TRISHNA]))
```

VRE gallery notebooks (10 tutorials, executed): `notebooks/` (see `notebooks/README.md`).
See `ROADMAP.md` for the next steps and `CHANGELOG.md` for the history.

## Documentation

https://trishna-tbx.readthedocs.io (sources in `docs/source`). Local build:

```bash
pip install -e ".[docs]"
sphinx-build -b html docs/source docs/build/html
```
