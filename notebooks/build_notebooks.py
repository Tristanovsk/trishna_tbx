"""Generate the TRISHNA notebooks of the Water Quality VRE gallery.

    cd notebooks && python build_notebooks.py           # write the .ipynb files
    python build_notebooks.py --execute                 # write and execute them (outputs kept)

Gallery folders follow the VRE categories (0_basics, 1_importdata, 3_atmosphericcorrection,
4_visualisation) plus a new category 6_trishna_simulations.
"""
import subprocess
import sys
from pathlib import Path

from _nbtools import code, data_cells, header, imports, md, write

HERE = Path(__file__).parent
NOTEBOOKS = {}


def register(path):
    def deco(f):
        NOTEBOOKS[path] = f
        return f
    return deco


# =============================================================================================
@register('0_basics/trishna_vs_sentinel2.ipynb')
def nb00():
    return header(
        'TRISHNA vs Sentinel-2: spectral bands, responses and noise',
        ['No image needed (spectral response functions and reference spectra are packaged).'],
        ['Know the TRISHNA VSWIR bands and how they compare with Sentinel-2/MSI.',
         'See which water-colour bands are missing (443 nm, red edge, 2190 nm) and why it matters.',
         'Understand that "the closest S2 band" is not the TRISHNA band (width, centre).',
         'Compare the radiometric noise of both sensors at 60 m.']) + imports() + [
        md('''
## 1. The TRISHNA VSWIR bands
TRISHNA carries an ISRO VSWIR instrument with 7 bands at 60 m (57-91 m across track)
and the CNES thermal instrument (4 TIR bands). Values from Table 6 of the TRISHNA mission
paper (EGUsphere-2026-4492). Signal-to-noise ratios are specified at the minimum radiance
`Lmin` (W m-2 sr-1 um-1).'''),
        code('''
from trishna_tbx.sensor import TRISHNA_BANDS, S2_BANDS
bands = pd.DataFrame([vars(b) for b in TRISHNA_BANDS]).rename(
    columns={'center': 'centre (nm)', 'fwhm': 'FWHM (nm)', 'lref': 'Lmin', 'snr': 'SNR @ Lmin'})
bands
'''),
        md('''
## 2. Spectral responses
The measured TRISHNA responses are not public yet: `trishna_tbx` models them with
super-Gaussian (near rectangular) shapes. Sentinel-2A responses are the ESA ones (from
GRSdriver). The mean Rrs of three optical water types (Bi et al. 2024) are shown for
context: clear (1), moderately turbid (4a) and very turbid (6).'''),
        code('''
from trishna_tbx import owt
s2 = make_s2()
ref, _ = owt.load_reference('Bi2024')

fig, axs = plt.subplots(2, 1, figsize=(12, 7), sharex=False)
for ax, rng in zip(axs, [(400, 1700), (400, 950)]):
    for b in s2.band_names:
        ax.fill_between(s2.srf.wl, s2.srf.sel(band=b), color='0.6', alpha=0.35, lw=0)
    for b in TRISHNA.band_names:
        ax.plot(TRISHNA.srf.wl, TRISHNA.srf.sel(band=b), lw=2, label=b if rng[1] > 1000 else None)
    ax.set_xlim(rng); ax.set_ylabel('relative response')
axr = axs[1].twinx()
for t, c in zip(['1', '4a', '6'], ['tab:blue', 'tab:olive', 'tab:brown']):
    axr.plot(ref.wl, ref.sel(owt=t), c=c, ls='--', label=f'OWT {t}')
axr.set_ylabel('Rrs (sr$^{-1}$)'); axr.legend(loc='upper right')
axs[0].legend(ncol=7, loc='upper center', bbox_to_anchor=(0.5, 1.2))
axs[0].set_title('TRISHNA (colour) vs Sentinel-2 (grey)', y=1.2)
axs[1].set_xlabel('wavelength (nm)')
fig.tight_layout()
fig.savefig(THUMB_DIR / 'trishna_vs_s2_srf.png', dpi=80)
'''),
        md('''
## 3. Closest Sentinel-2 band is not the TRISHNA band
The response-weighted centre of S2 B2 is 496.5 nm, about 11 nm redder than TRISHNA Blue.
S2 B3/B4 are half as wide as TRISHNA Green/Red. Most importantly, the TRISHNA water-vapour
band (910 nm) has **no S2 equivalent**: S2 B9 is centred at 945 nm, in a much deeper
absorption.'''),
        code('''
def centroid(srf):
    return float((srf * srf.wl).sum() / srf.sum())

def width(srf):
    hi = srf.wl.where(srf >= 0.5 * srf.max(), drop=True)
    return float(hi.max() - hi.min())

rows = []
for b in TRISHNA_BANDS:
    t, s = TRISHNA.srf.sel(band=b.name), s2.srf.sel(band=b.s2_proxy)
    rows.append({'TRISHNA': b.name, 'centre': round(centroid(t), 1), 'FWHM': round(width(t)),
                 'S2 proxy': b.s2_proxy, 'S2 centre': round(centroid(s), 1), 'S2 FWHM': round(width(s)),
                 'shift (nm)': round(centroid(t) - centroid(s), 1)})
pd.DataFrame(rows).set_index('TRISHNA')
'''),
        md('''
## 4. What the missing bands mean for water colour
| missing | typical use over water |
|---|---|
| 443 nm (S2 B1) | chlorophyll/CDOM absorption: OC3, QAA, CDOM algorithms |
| 705, 740, 783 nm (S2 B5-B7) | red-edge Chl-a in productive waters (NDCI, Gons, Gilerson), SPM in very turbid waters (Binding) |
| 2190 nm (S2 B12) | aerosol/glint over turbid waters, masking (GRS uses 1610/2190) |

The available bands keep the blue/green ratio (clear-water Chl-a), red and NIR (SPM,
turbidity), and add a 910 nm water-vapour band and the 1380 nm cirrus band.'''),
        md('''
## 5. Radiometric noise
Noise-equivalent TOA reflectance `NEdrho = pi * L/SNR / (E0 cos(sza))` at sza = 30 deg.
Sentinel-2 bands are degraded to 60 m by averaging 36 (10 m bands) or 9 (20 m bands) pixels,
which reduces white noise by a factor 6 or 3. **Caveat:** the S2 SNR is specified at a
reference radiance (bright) and the TRISHNA one at Lmin (dark): this comparison is indicative.
Rrs noise is roughly `NEdrho / (pi * T)` with T ~ 0.75 the two-way diffuse transmittance.'''),
        code('''
native = {'B1': 60, 'B2': 10, 'B3': 10, 'B4': 10, 'B5': 20, 'B6': 20, 'B7': 20, 'B8': 10, 'B8A': 20,
          'B9': 60, 'B10': 60, 'B11': 20, 'B12': 20}
n_agg = [(60 // native[b]) ** 2 for b in s2.band_names]
noise = pd.concat({
    'TRISHNA (60 m)': TRISHNA.noise_reflectance(30.).assign_coords(band=TRISHNA.band_names).to_pandas(),
    'S2 native': s2.noise_reflectance(30.).assign_coords(band=s2.band_names).to_pandas(),
    'S2 at 60 m': s2.noise_reflectance(30., n_aggregated=n_agg).assign_coords(band=s2.band_names).to_pandas(),
}, axis=1)
noise['Rrs noise (sr-1)'] = noise.max(axis=1) / (np.pi * 0.75)
(noise * 1e4).round(1).rename(columns=lambda c: c + ' x1e4')
'''),
        code('''
fig, ax = plt.subplots(figsize=(10, 4))
x = TRISHNA.centers.values
ax.semilogy(x, noise.loc[TRISHNA.band_names, 'TRISHNA (60 m)'], 'o-', label='TRISHNA (60 m)')
xs = s2.centers.values
ax.semilogy(xs, noise.loc[s2.band_names, 'S2 native'], 's--', c='0.5', label='S2 native resolution')
ax.semilogy(xs, noise.loc[s2.band_names, 'S2 at 60 m'], 's-', c='k', label='S2 averaged to 60 m')
ax.set_xlabel('wavelength (nm)'); ax.set_ylabel('NE$\\\\Delta\\\\rho_{TOA}$'); ax.legend(); ax.grid(alpha=.3)
'''),
        md('''
## Take-home messages
- TRISHNA keeps 4 water-colour bands (485/555/670/860 nm), wider than S2, plus WV (910),
  cirrus (1380) and SWIR (1610) for atmospheric correction and masking.
- Algorithms relying on 443 nm or on the red edge cannot be applied as such.
- Emulating TRISHNA by keeping S2 bands is an approximation (band shift/width): see
  *Emulate a TRISHNA L1C from Sentinel-2* and *Top-of-atmosphere simulations with RadCalNet-OC*.
- At requirement level, TRISHNA noise is comparable to or higher than S2 averaged to 60 m,
  in particular in the blue and WV bands.

## References
- TRISHNA mission paper, EGUsphere-2026-4492 (2026), Table 6.
- Bi, S. et al. (2024), optical water types (10 classes), mean spectra from GRSl2bgen.
- ESA, Sentinel-2 MSI spectral responses (COPE-GSEG-EOPG-TN-15-0007).'''),
    ]


# =============================================================================================
@register('1_importdata/emulate_trishna_l1c.ipynb')
def nb01():
    return header(
        'Emulate a TRISHNA L1C from Sentinel-2',
        ['An L1C Sentinel-2 image (SAFE).'],
        ['Load a Sentinel-2 L1C restricted to the TRISHNA-like bands, at 60 m, with GRSdriver.',
         'Compare with the full S2 image (RGB, TOA spectra).',
         'Add TRISHNA radiometric noise and look at its impact over water.',
         'Prepare the product for the GRS atmospheric correction (next tutorial).'],
        extra_env=['rioxarray'], resources='1 core, 4 GB', runtime='~2 min') + imports('''
from trishna_tbx.emulate import open_s2_l1c_as_trishna, grs_band_set, add_noise, S2_BAND_ORDER
from GRSdriver import Sentinel2Driver
''') + data_cells() + [
        md('''
> This tutorial complements *Simulate TRISHNA data* (Import Data section), which shows the
> band selection on a full Sentinel-2 image: here we add the 60 m resampling, the GRS-ready
> band configuration and the radiometric noise.

## 1. Load the image
`grs_band_set()` gives the S2 bands to load (TRISHNA proxies B2, B3, B4, B8A, B9, B10, B11
plus B12, kept for the GRS masking only) and the bands processed by GRS (B2, B3, B4, B8A, B11).'''),
        code('''
band_idx, tbp_idx = grs_band_set(strict=False)
print('loaded   :', [S2_BAND_ORDER[i] for i in band_idx])
print('processed:', [S2_BAND_ORDER[i] for i in tbp_idx])

tri = open_s2_l1c_as_trishna(S2_L1C, resolution=60, subset=SUBSET)
full = Sentinel2Driver(S2_L1C, resolution=60, resolution_angle=60, subset=SUBSET)
full.load_product()
full = full.prod
tri
'''),
        md('## 2. Plot the product'),
        code('''
def rgb(prod, wls=(665, 560, 490), gain=3.5):
    img = np.stack([prod.bands.sel(wl=w).values for w in wls], -1) * gain
    return np.clip(img ** (1 / 1.8), 0, 1)

fig, axs = plt.subplots(1, 2, figsize=(13, 6.5), sharex=True, sharey=True)
ext = [float(tri.x.min()), float(tri.x.max()), float(tri.y.min()), float(tri.y.max())]
for ax, p, t in zip(axs, [full, tri], ['Sentinel-2 L1C, 60 m', 'TRISHNA-like L1C, 60 m (B4/B3/B2)']):
    ax.imshow(rgb(p), extent=ext); ax.set_title(t)
fig.savefig(THUMB_DIR / 'emulate_trishna_l1c.png', dpi=70)
'''),
        md('''
## 3. Top-of-atmosphere spectra
Spectra of a water pixel (Berre lagoon), a marine pixel (Gulf of Fos) and a land pixel.
The TRISHNA-like product keeps 8 of the 13 S2 bands. Note that the S2 "WV" band is at
945 nm, not at the TRISHNA 910 nm. On this date the water is strongly affected by sunglint
(TOA reflectance ~0.06 at 1610 nm over water): this is what GRS corrects in the next tutorial.'''),
        code('''
points = {'Berre lagoon': (673000, 4814000), 'Gulf of Fos (sea)': (656000, 4802000), 'land': (680000, 4825000)}
fig, ax = plt.subplots(figsize=(10, 4.5))
for (name, (x, y)), c in zip(points.items(), ['tab:blue', 'tab:cyan', 'tab:green']):
    sf = full.bands.sel(x=x, y=y, method='nearest')
    st = tri.bands.sel(x=x, y=y, method='nearest')
    ax.plot(sf.wl.astype(float), sf, 'o--', c=c, alpha=.5, mfc='none', label=f'{name} - S2')
    ax.plot(st.wl.astype(float), st, 's', c=c, ms=8, label=f'{name} - TRISHNA-like')
ax.set_xscale('log'); ax.set_xticks([443, 560, 665, 865, 945, 1375, 1610, 2190])
ax.get_xaxis().set_major_formatter(plt.ScalarFormatter())
ax.set_xlabel('wavelength (nm)'); ax.set_ylabel('TOA reflectance'); ax.legend(fontsize=8, ncol=2)
'''),
        md('''
## 4. Add TRISHNA radiometric noise
White Gaussian noise from the SNR requirement (`trishna_tbx.sensor`). The S2 image at 60 m is
already an average of several pixels and has very little noise; this step mimics the TRISHNA
noise level. Over water, the spatial variability includes sunglint and real water structures:
the noise is estimated from the pixel-to-pixel differences, `std(diff) / sqrt(2)`, which removes
most of the smooth spatial variations.'''),
        code('''
rho = tri.bands.sel(wl=[490, 560, 665, 865]).astype(float)
names = ['Blue', 'Green', 'Red', 'NIR']
rho_t = rho.assign_coords(band=('wl', names)).swap_dims(wl='band')
noisy = add_noise(rho_t, tri.sza, sensor=TRISHNA.subset(names), band_dim='band', seed=1)

box = dict(x=slice(670000, 676000), y=slice(4815000, 4810000))   # Berre lagoon
lagoon = rho_t.sel(band='NIR') < 0.12                                # drop the few land pixels
diff_noise = lambda r: (r.where(lagoon).sel(**box).diff('x') / np.sqrt(2)).std(['x', 'y']).values
stats = pd.DataFrame({
    'mean TOA': rho_t.sel(**box).mean(['x', 'y']).values,
    'noise S2 at 60 m': diff_noise(rho_t),
    'noise TRISHNA-like + noise': diff_noise(noisy),
    'TRISHNA requirement': TRISHNA.subset(names).noise_reflectance(float(tri.sza.mean())).values}, index=names)
stats.round(5)
'''),
        code('''
fig, axs = plt.subplots(1, 2, figsize=(12, 5), sharex=True, sharey=True)
for ax, r, t in zip(axs, [rho_t, noisy], ['S2 B2 at 60 m (Berre lagoon)', 'with TRISHNA Blue noise']):
    im = ax.imshow(r.sel(band='Blue').sel(**box), vmin=0.10, vmax=0.16, cmap='viridis')
    ax.set_title(t)
fig.colorbar(im, ax=axs, shrink=.7, label='TOA reflectance')
'''),
        md('''
## 5. Save for the next steps
The product (`bands`, angles, `wl_to_process`) can be passed directly to
`grs.Process().execute(prod, ...)`: see *Atmospheric correction with TRISHNA bands (GRS)*.'''),
        code('''
out = OUTPUT_DIR / 'trishna_like_l1c_subset.nc'
tri_save = tri.drop_vars(['SRF', 'wl_hr'], errors='ignore')
tri_save = tri_save.assign_coords(wl=tri_save.wl.astype(int))
tri_save.attrs = {k: (str(v) if not isinstance(v, (int, float, str)) else v) for k, v in tri_save.attrs.items()}
tri_save.to_netcdf(out)
print(out)
'''),
        md('''
## Take-home messages
- Emulating a TRISHNA L1C from S2 = band selection + 60 m averaging (+ noise). The shift of the
  blue band (11 nm) and the 945 vs 910 nm water-vapour band are **not** corrected by this step.
- At 60 m the noise of the emulated image is lower than the TRISHNA requirement: add it for
  realistic performance studies.'''),
    ]


# =============================================================================================
GRS_HELPER = '''
import logging
import grs
logging.getLogger().setLevel(logging.WARNING)   # GRS logs at INFO level

def clean_attrs(ds):
    """NetCDF-safe attributes (GRS attributes may hold booleans/lists)."""
    ds = ds.copy()
    for v in ds.variables.values():
        v.attrs = {k: (str(a) if not isinstance(a, (int, float, str)) or isinstance(a, bool) else a)
                   for k, a in v.attrs.items()}
    ds.attrs = {k: (str(a) if not isinstance(a, (int, float, str)) or isinstance(a, bool) else a)
                for k, a in ds.attrs.items()}
    return ds

def run_grs(prod, name, **kw):
    """Run GRS at 60 m (or reload the cached result)."""
    f = OUTPUT_DIR / f'grs_{name}.nc'
    if f.exists():
        return xr.open_dataset(f)
    p = grs.Process()
    p.execute(prod, odir=str(OUTPUT_DIR / name), cams_file=CAMS_FILE, resolution=60, **kw)
    l2 = p.l2a.l2_prod[['Rrs', 'aot550', 'mask', 'flags']].load()
    l2 = l2.assign_coords(wl=l2.wl.astype(int))
    clean_attrs(l2).to_netcdf(f)
    return xr.open_dataset(f)
'''


@register('3_atmosphericcorrection/grs_trishna_bands.ipynb')
def nb02():
    return header(
        'Atmospheric correction with TRISHNA bands (GRS)',
        ['An L1C Sentinel-2 image (SAFE) and the matching CAMS file (GRS input).',
         'The GRS look-up tables (downloaded by `grs` at first use).'],
        ['Run GRS on a Sentinel-2 image restricted to the TRISHNA-like bands.',
         'Compare with GRS run on all the S2 bands: Rrs and aerosol optical thickness.',
         'Quantify the sensitivity of the TRISHNA-band Rrs to an aerosol error.'],
        resources='4 cores, 8 GB', runtime='~25 min for the 3 GRS runs (cached in OUTPUT_DIR)') + imports('''
from GRSdriver import Sentinel2Driver
from trishna_tbx.emulate import open_s2_l1c_as_trishna
''') + data_cells() + [
        md('''
## 1. What changes for GRS with TRISHNA bands
GRS (Harmel et al. 2018) retrieves the aerosol load and the sunglint from the **two
longest-wavelength processed bands**:
- full Sentinel-2: 1610 and 2190 nm (SWIR),
- TRISHNA: 865 and 1610 nm (NIR + SWIR), since there is no 2190 nm band.

The 2190 nm band is still loaded (not processed) because the GRS "high SWIR" mask needs it;
this dependency has to be removed in GRS for real TRISHNA data. Water vapour comes from CAMS
in both cases (see *Water vapour with TRISHNA*).'''),
        code(GRS_HELPER),
        md('## 2. Load the L1C products and run GRS'),
        code('''
prod_tri = open_s2_l1c_as_trishna(S2_L1C, resolution=60, subset=SUBSET)
drv = Sentinel2Driver(S2_L1C, resolution=60, resolution_angle=60, subset=SUBSET)
drv.load_product()
prod_s2 = drv.prod
print('TRISHNA-like processed bands:', prod_tri.attrs['wl_to_process'])
print('S2 processed bands          :', prod_s2.attrs['wl_to_process'])
'''),
        code('''
%%time
l2_tri = run_grs(prod_tri, 'trishna')
l2_s2 = run_grs(prod_s2, 's2')
'''),
        md('''
## 3. Rrs maps
Water pixels are those flagged as water by GRS (`mask == 0`) in both runs.'''),
        code('''
water = (l2_tri.mask == 0) & (l2_s2.mask == 0)
common = [490, 560, 665, 865]
fig, axs = plt.subplots(3, 4, figsize=(16, 11), sharex=True, sharey=True)
for j, w in enumerate(common):
    a = l2_s2.Rrs.sel(wl=w).where(water)
    b = l2_tri.Rrs.sel(wl=w).where(water)
    vmax = float(a.quantile(0.98))
    a.plot(ax=axs[0, j], vmin=0, vmax=vmax, cmap='viridis', add_labels=False)
    b.plot(ax=axs[1, j], vmin=0, vmax=vmax, cmap='viridis', add_labels=False)
    (100 * (b - a) / a).plot(ax=axs[2, j], vmin=-30, vmax=30, cmap='RdBu_r', add_labels=False,
                              cbar_kwargs={'label': '%'})
    axs[0, j].set_title(f'{w} nm')
for i, t in enumerate(['GRS, all S2 bands', 'GRS, TRISHNA bands', 'difference (%)']):
    axs[i, 0].set_ylabel(t)
for ax in axs.ravel():
    ax.set_xticks([]); ax.set_yticks([])
fig.tight_layout()
fig.savefig(THUMB_DIR / 'grs_trishna_bands.png', dpi=60)
'''),
        md('## 4. Band-by-band comparison and aerosol optical thickness'),
        code('''
rows = {}
fig, axs = plt.subplots(1, 5, figsize=(19, 4))
for ax, w in zip(axs, common):
    a = l2_s2.Rrs.sel(wl=w).where(water).values.ravel()
    b = l2_tri.Rrs.sel(wl=w).where(water).values.ravel()
    ok = np.isfinite(a) & np.isfinite(b)
    ax.hexbin(a[ok], b[ok], gridsize=60, bins='log', cmap='magma_r')
    lim = np.nanpercentile(a[ok], [0.5, 99.5]); ax.plot(lim, lim, 'k--', lw=1)
    ax.set_xlim(lim); ax.set_ylim(lim); ax.set_title(f'Rrs {w} nm'); ax.set_xlabel('all S2 bands')
    rows[w] = {'median bias (sr-1)': np.median(b[ok] - a[ok]), 'MdAPD (%)': 100 * np.median(np.abs(b[ok] - a[ok]) / np.abs(a[ok])),
               'r': np.corrcoef(a[ok], b[ok])[0, 1], 'N': ok.sum()}
axs[0].set_ylabel('TRISHNA bands')
a, b = l2_s2.aot550.where(water).values.ravel(), l2_tri.aot550.where(water).values.ravel()
ok = np.isfinite(a) & np.isfinite(b)
axs[4].hexbin(a[ok], b[ok], gridsize=40, bins='log', cmap='magma_r'); axs[4].plot([0, 1], [0, 1], 'k--')
axs[4].set_xlim(0, np.nanpercentile(a, 99)); axs[4].set_ylim(0, np.nanpercentile(a, 99)); axs[4].set_title('AOT(550)')
fig.tight_layout()
pd.DataFrame(rows).T.round(5)
'''),
        code('''
from pyproj import Transformer
from trishna_tbx.io import cams_value
lon, lat = Transformer.from_crs(32631, 4326, always_xy=True).transform(float(l2_tri.x.mean()), float(l2_tri.y.mean()))
aod_cams = float(cams_value(CAMS_FILE, 'aod550', prod_tri.time.values, lon, lat)[0])
pd.Series({'CAMS aod550 (first guess)': aod_cams,
           'GRS AOT, TRISHNA bands (865/1610)': float(l2_tri.aot550.where(water).median()),
           'GRS AOT, all S2 bands (1610/2190)': float(l2_s2.aot550.where(water).median())}).round(3)
'''),
        md('''
On the default scene (26 May 2025, strong sunglint over the Berre lagoon), the two band sets
give very different AOT: the TRISHNA pair (865/1610 nm) stays close to CAMS while the S2 SWIR
pair (1610/2190 nm) gives a much higher value, and the red/NIR Rrs differ accordingly. The aerosol
and the sunglint are retrieved jointly from these two bands: under strong glint the partition
between both is sensitive to the band pair. This single case does not tell which is right:
matchups with AERONET-OC / HYPERNETS are needed (see the *Validation* section).'''),
        md('''
## 5. Sensitivity to the aerosol model and to the sunglint handling
Two GRS options are tested with both band sets:
- `opac_model='MACL_rh70'`: force a maritime aerosol model (by default it is chosen from CAMS);
- `glint_aot=True`: hybrid AOT retrieval for blocks under strong sunglint (recent grs versions only).

For a broader sensitivity study (AOT, aerosol type, H2O, O3) at lower cost, use the
simulations of *Top-of-atmosphere simulations with RadCalNet-OC*.'''),
        code('''
%%time
import inspect
runs = {('TRISHNA', 'default'): l2_tri, ('S2', 'default'): l2_s2,
        ('TRISHNA', 'maritime aerosol'): run_grs(prod_tri, 'trishna_macl', opac_model='MACL_rh70')}
if 'glint_aot' in inspect.signature(grs.Process.execute).parameters:
    runs[('TRISHNA', 'glint_aot')] = run_grs(prod_tri, 'trishna_glintaot', glint_aot=True)
    runs[('S2', 'glint_aot')] = run_grs(prod_s2, 's2_glintaot', glint_aot=True)
else:
    print(f'grs {grs.__version__}: no glint_aot option (needs a recent grs), skipped')
rows = {}
for (bands_, case), l2 in runs.items():
    w = water & (l2.mask == 0)
    r = {f'Rrs {wl}': float(l2.Rrs.sel(wl=wl).where(w).median()) for wl in common}
    r['AOT(550)'] = float(l2.aot550.where(w).median())
    rows[(bands_, case)] = r
pd.DataFrame(rows).T.round(5)
'''),
        md('''
## Take-home messages
- GRS runs on the TRISHNA band set without code change, using 865/1610 nm for the aerosol
  and glint retrieval.
- The aerosol/glint band pair matters: under strong glint the AOT and the red/NIR Rrs can differ
  a lot between the TRISHNA and S2 band sets; in turbid waters the 865 nm water signal is not
  zero and can bias the TRISHNA aerosol retrieval. Validation on in situ matchups is needed.
- For real TRISHNA products, GRS needs: no dependency on 2190 nm for masking, the TRISHNA
  SRF (`wl_true`), and an image-based water vapour (910 nm band).

## References
- Harmel, T. et al. (2018), Sunglint correction of the Multi-Spectral Instrument (MSI)-SENTINEL-2 imagery over inland and sea waters from SWIR bands, RSE 204.'''),
    ]


# =============================================================================================
@register('3_atmosphericcorrection/water_vapor_trishna.ipynb')
def nb04():
    return header(
        'Water vapour with TRISHNA: 860/910 nm vs Sentinel-2 865/945 nm',
        ['An L1C Sentinel-2 image (SAFE) and a CAMS file (for comparison).'],
        ['Understand the two-band ratio retrieval of the total column water vapour (TCWV).',
         'Compare the sensitivity of the TRISHNA (910 nm) and S2 (945 nm) absorption bands.',
         'Estimate the TCWV error due to radiometric noise.',
         'Apply the ratio method to a Sentinel-2 image and build a pseudo-TRISHNA test.'],
        extra_env=['scipy', 'pyproj'], runtime='~3 min') + imports('''
from GRSdriver import Sentinel2Driver
from trishna_tbx.water_vapor import (RatioRetrieval, h2o_optical_thickness, compare_trishna_s2,
                                     noise_error, W_GRID)
from trishna_tbx.io import cams_value
''') + data_cells() + [
        md('''
## 1. Principle
The TOA reflectance in an absorption band (910 or 945 nm) divided by a window band
(860/865 nm) is, for a flat surface reflectance and negligible path radiance,

$$\\frac{R_{abs}}{R_{win}} = \\frac{T_w(\\lambda_{abs}, M W)}{T_w(\\lambda_{win}, M W)}, \\qquad M = \\frac{1}{\\cos\\theta_s} + \\frac{1}{\\cos\\theta_v}$$

The transmittances are computed with the REPTRAN H2O absorption (radcalnet_oc) and integrated
over the band responses.'''),
        code('''
tau = h2o_optical_thickness().sel(wl=slice(820, 1000))
s2 = make_s2()
fig, ax = plt.subplots(figsize=(11, 4.5))
for W, c in zip([10, 30, 50], ['0.7', '0.45', '0.1']):
    T = np.exp(-3 * W * tau).rolling(wl=15, center=True).mean()   # ~3 nm smoothing for display
    ax.plot(tau.wl, T, c=c, lw=1, label=f'T$_w$, W={W} kg m$^{{-2}}$, M=3')
for b, c in [('NIR', 'tab:red'), ('WV', 'tab:red')]:
    ax.fill_between(TRISHNA.srf.wl, TRISHNA.srf.sel(band=b), color=c, alpha=.3, label=f'TRISHNA {b}')
for b in ['B8A', 'B9']:
    ax.fill_between(s2.srf.wl, s2.srf.sel(band=b), color='tab:blue', alpha=.3, label=f'S2 {b}')
ax.set_xlim(820, 1000); ax.set_xlabel('wavelength (nm)'); ax.legend(fontsize=8, loc='center left', bbox_to_anchor=(1.01, 0.5))
fig.tight_layout()
fig.savefig(THUMB_DIR / 'water_vapor_trishna.png', dpi=70)
'''),
        md('## 2. Ratio vs water vapour'),
        code('''
fig, axs = plt.subplots(1, 2, figsize=(13, 4.5))
for s, ls in [('TRISHNA', '-'), ('S2', '--')]:
    r = RatioRetrieval(s)
    for M, c in zip([2.2, 3, 4], ['tab:green', 'tab:orange', 'tab:purple']):
        r.ratio_lut.sel(M=M, method='nearest').plot(ax=axs[0], ls=ls, c=c, label=f'{s}, M={M}')
c = compare_trishna_s2(3.)
for s, ls in [('TRISHNA', '-'), ('S2', '--')]:
    c.dW_for_1pct_ratio_error.sel(sensor=s).plot(ax=axs[1], ls=ls, c='k', label=s)
axs[0].set_ylabel('R$_{abs}$ / R$_{win}$'); axs[0].legend(fontsize=8, ncol=2)
axs[1].set_ylabel('TCWV error for a 1 % ratio error (kg m$^{-2}$)'); axs[1].legend()
for ax in axs: ax.set_xlabel('TCWV W (kg m$^{-2}$)'); ax.set_title('')
'''),
        md('''
The 910 nm band absorbs less: its ratio does not saturate in humid atmospheres (S2 ratio
falls to ~0.03 for W = 60 kg m-2), which motivated the TRISHNA choice. The price is a lower
sensitivity per unit of W, so errors on the ratio (noise, surface spectral slope) cost more.

## 3. Noise-propagated TCWV error
For a flat surface reflectance `rho`, with the requirement SNRs (TRISHNA WV: 50 at Lmin;
S2 B9: 114 at 9 W m-2 sr-1 um-1, native 60 m). Path radiance is ignored, so the dark
cases (rho = 0.01, water) are optimistic.'''),
        code('''
fig, ax = plt.subplots(figsize=(9, 4.5))
for rho, c in zip([0.3, 0.05, 0.01], ['tab:green', 'tab:olive', 'tab:blue']):
    for s, ls in [('TRISHNA', '-'), ('S2', '--')]:
        noise_error(s, rho, m=3., sza=48.19).plot(ax=ax, c=c, ls=ls, label=f'{s}, rho={rho}')
ax.set_yscale('log'); ax.set_ylabel('1-sigma TCWV error (kg m$^{-2}$)'); ax.set_xlabel('TCWV (kg m$^{-2}$)')
ax.legend(fontsize=8, ncol=2); ax.grid(alpha=.3); ax.set_title('')
'''),
        md('''
## 4. Application to a Sentinel-2 image (865/945 nm)
TOA reflectances at 60 m, per-pixel air mass from the GRSdriver angles. The CAMS TCWV at the
image centre is given for comparison.'''),
        code('''
drv = Sentinel2Driver(S2_L1C, band_idx=[2, 8, 9], resolution=60, resolution_angle=60, subset=SUBSET)
drv.load_product()
p = drv.prod
r865, r945 = p.bands.sel(wl=865).astype(float), p.bands.sel(wl=945).astype(float)
vza = p.vza.sel(wl=945)
water = ((p.bands.sel(wl=560) - r865) / (p.bands.sel(wl=560) + r865)) > 0.05

W_s2 = xr.DataArray(RatioRetrieval('S2').retrieve(r865.values, r945.values, p.sza.values, vza.values),
                    dims=('y', 'x'), coords={'y': p.y, 'x': p.x})

from pyproj import Transformer
lon, lat = Transformer.from_crs(32631, 4326, always_xy=True).transform(float(p.x.mean()), float(p.y.mean()))
W_cams = float(cams_value(CAMS_FILE, 'tcwv', p.time.values, lon, lat)[0])
print(f'CAMS TCWV: {W_cams:.1f} kg m-2;  S2 ratio retrieval, land median: {float(W_s2.where(~water).median()):.1f};'
      f' water median: {float(W_s2.where(water).median()):.1f}')
'''),
        code('''
fig, axs = plt.subplots(1, 2, figsize=(13, 5))
W_s2.plot(ax=axs[0], vmin=0, vmax=2 * W_cams, cmap='Spectral_r', cbar_kwargs={'label': 'TCWV (kg m$^{-2}$)'})
axs[0].set_title('S2 865/945 ratio retrieval'); axs[0].set_aspect('equal')
for m, lab, c in [(~water, 'land', 'tab:green'), (water, 'water', 'tab:blue')]:
    axs[1].hist(W_s2.where(m).values.ravel(), bins=np.linspace(0, 2.5 * W_cams, 80), histtype='step', color=c, label=lab)
axs[1].axvline(W_cams, c='k', ls='--', label='CAMS'); axs[1].set_xlabel('TCWV (kg m$^{-2}$)'); axs[1].legend()
'''),
        md('''
Over water, the surface is dark and the signal at 865/945 nm is dominated by the path
radiance, which crosses less water vapour than the surface-reflected light: the decoupled
ratio model underestimates W. **Recommended approach**: retrieve over land/bright pixels and
interpolate over water (W varies slowly), or use the coupled path-length model of hgrs.'''),
        code('''
from scipy import ndimage
land = (~water & np.isfinite(W_s2)).values
idx = ndimage.distance_transform_edt(~land, return_distances=False, return_indices=True)
W_filled = W_s2.copy(data=ndimage.median_filter(W_s2.values[tuple(idx)], size=9))
W_filled.plot(vmin=0, vmax=2 * W_cams, cmap='Spectral_r', figsize=(6.5, 5), cbar_kwargs={'label': 'TCWV (kg m$^{-2}$)'})
plt.gca().set_aspect('equal'); plt.title('land retrieval interpolated over water')
'''),
        md('''
## 5. Pseudo-TRISHNA test
The 910 nm band cannot be emulated from Sentinel-2. We build a pseudo-TRISHNA scene: the
S2-derived TCWV (land) is taken as truth, the 860 nm reflectance is the S2 B8A value, and the
910 nm reflectance is computed with the TRISHNA ratio model plus TRISHNA noise. The same is
done for S2 945 nm with S2 noise, so that both retrievals are compared on the same scene.'''),
        code('''
rng = np.random.default_rng(0)
M = 1 / np.cos(np.radians(p.sza)) + 1 / np.cos(np.radians(vza))
truth = W_filled.where(~water)
res = {}
for s in ['TRISHNA', 'S2']:
    rr = RatioRetrieval(s)
    sens = TRISHNA if s == 'TRISHNA' else s2
    nz = sens.noise_reflectance(float(p.sza.mean())).assign_coords(band=sens.band_names)
    ratio = rr.ratio_lut.interp(W=truth.clip(0.5, 79), M=M.clip(2, 8))
    r_abs = r865 * ratio + float(nz.sel(band=rr.abs)) * rng.standard_normal(r865.shape)
    r_win = r865 + float(nz.sel(band=rr.win)) * rng.standard_normal(r865.shape)
    W = rr.retrieve(r_win.values, r_abs.values, p.sza.values, vza.values)
    res[s] = (W - truth.values)
pd.DataFrame({s: {'bias': np.nanmean(e), 'std': np.nanstd(e), 'p95 |error|': np.nanpercentile(np.abs(e), 95)}
              for s, e in res.items()}).round(2)
'''),
        md('''
## Take-home messages
- 910 nm is less absorbed than 945 nm: no saturation in humid air, but 3-4x less sensitive.
- With the requirement SNR of the TRISHNA WV band (50 at Lmin), the noise-induced TCWV error
  is 2-4x larger than with S2 B9 over bright surfaces. The actual in-flight SNR is key.
- Over water, neither ratio retrieval works directly: use land retrievals + interpolation or
  a coupled radiative transfer model (hgrs).

## References
- Gao, B.-C. & Kaufman, Y. (2003), MODIS near-IR water vapor algorithm, JGR 108.
- TRISHNA mission paper, EGUsphere-2026-4492, section 3.3.3 (need for a water vapor band).'''),
    ]


# =============================================================================================
@register('6_trishna_simulations/toa_simulation_radcalnet_oc.ipynb')
def nb03():
    return header(
        'Top-of-atmosphere simulations with RadCalNet-OC for TRISHNA',
        ['The AERONET-OC database (`aeronet_oc_db`); optionally the hyperspectral invRrs '
         'reconstruction of the AERONET-OC spectra.',
         'The radcalnet_oc look-up tables (`radcalnet_oc_config --download`).'],
        ['Propagate in situ water-leaving reflectances (AERONET-OC, Venice) to the top of the atmosphere.',
         'Convolve the hyperspectral TOA signal with the TRISHNA and Sentinel-2 responses.',
         'Derive Sentinel-2 -> TRISHNA band adjustment factors at TOA.',
         'Analyse the sensitivity of the TRISHNA bands to aerosols, water vapour and aerosol type.'],
        resources='4 cores, 8 GB', runtime='~6 min') + imports('''
from trishna_tbx.simulate import aeronet_oc_input, simulate_toa, to_sensors, sensitivity_grid, build_input
from trishna_tbx.emulate import SBAF
from radcalnet_oc.process import AEROSOL_MODELS
''') + data_cells() + [
        md('''
## 1. In situ water reflectance: AERONET-OC Venice
RadCalNet-OC (radcalnet_oc) takes the water reflectance and the atmospheric state measured
at the site (aerosol optical thickness, water vapour, ozone, NO2, pressure, wind) for each
record. The SeaPRISM multispectral Rrs is replaced by its hyperspectral reconstruction
(invRrs) when available.'''),
        code('''
SITE, PERIOD = 'Venise', slice('2021-05', '2021-09')
hyp = AERONET_HYP_DIR if Path(AERONET_HYP_DIR).exists() else None
inp = aeronet_oc_input(SITE, PERIOD, hyp_dir=hyp, max_records=40)
print(f'{inp.sizes["time"]} records; hyperspectral: {hyp is not None}')
inp[['sza', 'aot550', 'ang_exp', 'tcwv', 'tco3', 'pressure', 'wind', 'chl_aeronet']].to_dataframe().describe().round(3)
'''),
        code('''
fig, ax = plt.subplots(figsize=(10, 4))
cmap = plt.cm.viridis
for i in range(inp.sizes['time']):
    c = cmap(np.log10(inp.chl_aeronet[i]) / 1.5)
    inp.Rrs.isel(time=i).plot(ax=ax, c=c, lw=.7, alpha=.7)
ax.set_xlim(380, 1000); ax.set_title(f'AERONET-OC {SITE}: Rrs (colour: Chl-a)'); ax.set_ylabel('Rrs (sr$^{-1}$)')
'''),
        md('## 2. Simulation at the top of the atmosphere'),
        code('''
%%time
toa = simulate_toa(inp, vza=(5., 20., 35.), azi=(90.,))
bands = to_sensors(toa)
toa
'''),
        code('''
i, v = 0, 1  # record, viewing angle index
fig, ax = plt.subplots(figsize=(12, 5))
toa.Rtoa.isel(time=i, vza=v).squeeze().plot(ax=ax, c='0.6', lw=.5, label='hyperspectral')
for name, mk, c in [('TRISHNA', 'o', 'tab:red'), ('Sentinel-2 A', 's', 'tab:blue')]:
    r = bands[name].Rtoa.isel(time=i, vza=v).squeeze()
    ax.plot(r.wl, r, mk, c=c, ms=8, mfc='none' if mk == 's' else c, label=name)
ax.set_xlim(400, 1700); ax.set_ylim(0, None); ax.legend()
ax.set_title(f'Simulated TOA reflectance, {str(toa.time.values[i])[:16]}, vza={float(toa.vza[v])}')
fig.savefig(THUMB_DIR / 'toa_simulation_radcalnet_oc.png', dpi=70)
'''),
        md('''
## 3. Sentinel-2 -> TRISHNA band adjustment at TOA
Factors fitted on the simulated TOA spectra (all records and angles). Contrary to the water
reflectance (see *Water quality products with TRISHNA*), the TOA factors vary little from one
case to another (narrow 95 % range) because the atmospheric signal dominates; the Blue factor
(~1.08) comes from the bluer TRISHNA band seeing more Rayleigh scattering. The cirrus band has
no factor here: its TOA signal is ~0 for clear-sky simulations. The WV factor (B9 945 nm -> 910 nm) is large and
depends on the water vapour: **a S2 B9 value is not a TRISHNA WV value**.'''),
        code('''
spec = toa.Rtoa.squeeze().stack(sample=[d for d in toa.Rtoa.squeeze().dims if d != 'wl']).transpose('sample', 'wl')
spec = spec.drop_vars([c for c in spec.coords if c not in ('wl',)])
SBAF(spec).table()
'''),
        md('''
## 4. Sensitivity to the aerosol optical thickness and water vapour
One Rrs spectrum is simulated for a grid of AOT(550) and TCWV in a single run. The water
contribution at TOA is `Rtoa - Ratm`.'''),
        code('''
%%time
grid = sensitivity_grid(inp.Rrs.isel(time=0), sza=float(inp.sza[0]),
                        aot=(0.02, 0.05, 0.1, 0.2, 0.4), tcwv=(5., 15., 30., 50.), vza=(20.,), azi=(90.,))
gb = to_sensors(grid, variables=('Rtoa', 'Ratm'))['TRISHNA']
frac = ((gb.Rtoa - gb.Ratm) / gb.Rtoa).squeeze()
'''),
        code('''
fig, axs = plt.subplots(1, 2, figsize=(13, 4.5))
for b in ['Blue', 'Green', 'Red', 'NIR']:
    (100 * frac.sel(band=b, tcwv_sim=15.)).plot(ax=axs[0], marker='o', label=b)
axs[0].set_ylabel('water contribution to TOA (%)'); axs[0].set_xlabel('AOT(550)'); axs[0].legend(); axs[0].set_title('')
for b in ['NIR', 'WV', 'Cirrus', 'SWIR']:
    (gb.Rtoa.sel(band=b, aot550_sim=0.1).squeeze() / gb.Rtoa.sel(band=b, aot550_sim=0.1, tcwv_sim=5.).squeeze()
     ).plot(ax=axs[1], marker='o', label=b)
axs[1].set_ylabel('Rtoa / Rtoa(W=5)'); axs[1].set_xlabel('TCWV (kg m$^{-2}$)'); axs[1].legend(); axs[1].set_title('')
'''),
        md('''
## 5. Sensitivity to the aerosol type
Same record with a single OPAC aerosol model (maritime clean vs urban vs desert).'''),
        code('''
%%time
one = inp.isel(time=[0])
res = {}
for m in ['MACL_rh70', 'URBA_rh70', 'DESE_rh70']:
    comb = [1. if a == m else 0. for a in AEROSOL_MODELS]
    t = simulate_toa(one, vza=(20.,), azi=(90.,), aerosol_combination=comb)
    res[m] = to_sensors(t, sensors=[TRISHNA])['TRISHNA'].Rtoa.squeeze().to_pandas()
df = pd.DataFrame(res)
df['urban - maritime (%)'] = 100 * (df['URBA_rh70'] - df['MACL_rh70']) / df['MACL_rh70']
df.round(4)
'''),
        md('''
## Take-home messages
- RadCalNet-OC propagates any water reflectance (in situ or simulated) to TOA for any band set,
  with the measured atmosphere: the basis for TRISHNA band-adjustment, noise and algorithm tests.
- In the blue, the water signal is only ~10 % of the TOA signal for moderate aerosol loads:
  aerosol type and AOT errors dominate the uncertainty of the TRISHNA Blue Rrs.
- The 910 nm band is sensitive to water vapour without saturating; S2 B9 cannot stand for it.

## References
- RadCalNet-OC / radcalnet_oc documentation (Magellium).
- Zibordi, G. et al. (2009), AERONET-OC: a network for the validation of ocean color primary products, JTECH 26.'''),
    ]


# =============================================================================================
L2A_LOAD = '''
from trishna_tbx.io import list_l2a, open_grs_l2a, water_mask, to_60m
files = list_l2a(L2A_DIR, L2A_TILE)
date, f = min(files, key=lambda tf: abs(tf[0] - np.datetime64(L2A_DATE)))
ds = open_grs_l2a(f, L2A_WINDOW)
rrs = to_60m(ds.Rrs.where(water_mask(ds))).load()
rrs = rrs.where((rrs.sel(wl=[490, 560, 665]) > 0).all('wl'))
print(date, f.split('/')[-1], '| 60 m water pixels:', int(rrs.sel(wl=560).notnull().sum()))
'''


@register('4_visualisation/optical_water_types_trishna.ipynb')
def nb05():
    return header(
        'Optical water types with TRISHNA bands',
        ['A GRS L2A Sentinel-2 product (section 4); sections 1-3 need no image.'],
        ['Classify water reflectances into the 10 optical water types of Bi et al. (2024).',
         'Measure how well the TRISHNA bands separate the water types compared with Sentinel-2.',
         'Define a reduced OWT set adapted to TRISHNA.',
         'Map the OWT on a Sentinel-2 image with all bands and with TRISHNA-like bands.'],
        runtime='~2 min') + imports('''
from trishna_tbx import owt
''') + data_cells() + [
        md('''
## 1. Reference spectra
The Bi et al. (2024) mean spectra (400-900 nm) are **convolved** with each sensor's
responses: interpolating them at the band centres is not adequate for the 60-70 nm wide
TRISHNA bands. Classification uses the spectral angle (SAM), as in GRSl2bgen.'''),
        code('''
s2 = make_s2()
S2W = s2.subset(['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8A'], 'S2 (8 water bands)')
S2T = s2.subset(['B2', 'B3', 'B4', 'B8A'], 'S2 restricted to TRISHNA bands')
TRW = TRISHNA.subset(['Blue', 'Green', 'Red', 'NIR'], 'TRISHNA')
ref, sd = owt.load_reference('Bi2024')

fig, axs = plt.subplots(2, 5, figsize=(17, 6), sharex=True)
for ax, t in zip(axs.ravel(), ref.owt.values):
    ax.plot(ref.wl, ref.sel(owt=t), c='0.4')
    for sens, mk, c in [(S2W, 's', 'tab:blue'), (TRW, 'o', 'tab:red')]:
        cv = owt.convolve_reference(ref.sel(owt=t), sens)
        ax.plot(cv.wl, cv, mk, c=c, mfc='none' if mk == 's' else c)
    ax.set_title(f'OWT {t}')
fig.supxlabel('wavelength (nm)'); fig.supylabel('Rrs (sr$^{-1}$)'); fig.tight_layout()
fig.savefig(THUMB_DIR / 'optical_water_types_trishna.png', dpi=60)
'''),
        md('''
## 2. Separability
Spectral angle (degrees) between the reference types: the smaller the angle, the harder it
is to tell the types apart.'''),
        code('''
fig, axs = plt.subplots(1, 2, figsize=(14, 5.5))
for ax, sens in zip(axs, [S2W, TRW]):
    sep = owt.separability(ref, sens).to_pandas()
    im = ax.imshow(sep, vmin=0, vmax=30, cmap='magma')
    ax.set_xticks(range(10), sep.columns); ax.set_yticks(range(10), sep.index)
    for i in range(10):
        for j in range(10):
            ax.text(j, i, f'{sep.iloc[i, j]:.0f}', ha='center', va='center', fontsize=7,
                    color='w' if sep.iloc[i, j] < 15 else 'k')
    ax.set_title(sens.name)
fig.colorbar(im, ax=axs, shrink=.8, label='spectral angle (deg)')
'''),
        md('''
## 3. Classification benchmark
Synthetic members of each type (mean + standard deviation with random amplitude and
spectral-slope perturbations) are classified with the hyperspectral spectrum (reference)
and with each band set. Noise: requirement SNR converted to Rrs (S2 averaged to 60 m).'''),
        code('''
sensors = {'S2': S2W, 'S2-as-TRISHNA': S2T, 'TRISHNA': TRW,
           'TRISHNA without NIR': TRISHNA.subset(['Blue', 'Green', 'Red'], 'T3')}
k = 1 / (np.pi * 0.75)
noise = {n: s.noise_reflectance(30, n_aggregated=1 if n.startswith('TRISHNA') else 9)
         .assign_coords(band=s.band_names) * k for n, s in sensors.items()}
clean, conf = owt.benchmark(sensors, n=300)
noisy, conf_noisy = owt.benchmark(sensors, n=300, noise=noise)
ax = (100 * pd.concat({'noise-free': clean.loc['all'], 'with noise': noisy.loc['all']}, axis=1)).plot.bar(
    figsize=(9, 4), rot=0)
ax.set_ylabel('agreement with hyperspectral OWT (%)'); ax.set_ylim(50, 100); ax.grid(axis='y', alpha=.3)
(100 * noisy).round(0).T
'''),
        code('''
fig, axs = plt.subplots(1, 2, figsize=(14, 5.5))
for ax, n in zip(axs, ['S2', 'TRISHNA']):
    c = conf_noisy[n]
    ax.imshow(100 * c.values, vmin=0, vmax=100, cmap='Blues')
    ax.set_xticks(range(c.shape[1]), c.columns); ax.set_yticks(range(c.shape[0]), c.index)
    for i in range(c.shape[0]):
        for j in range(c.shape[1]):
            if c.iloc[i, j] > 0.02:
                ax.text(j, i, f'{100 * c.iloc[i, j]:.0f}', ha='center', va='center', fontsize=7)
    ax.set_xlabel(f'{n} class'); ax.set_ylabel('hyperspectral class'); ax.set_title(f'{n} (with noise), %')
'''),
        md('''
## 4. A TRISHNA OWT set
With TRISHNA bands, the sub-types 1/2, 3a/3b and 4a/4b differ mostly in amplitude (angle
< 6 deg). Merging them gives a 7-class set that TRISHNA can resolve much better.'''),
        code('''
print(owt.TRISHNA_MERGE)
members = owt.sample_members(ref, sd, n=300)
hyp, _ = owt.classify(members.sel(wl=slice(400, 900)), ref.sel(wl=slice(400, 900)), dim='wl')
rng = np.random.default_rng(1)
rows = {}
for n in ['S2', 'TRISHNA']:
    sens = sensors[n]
    x = owt.convolve_reference(members.rename(true_owt='owt'), sens).rename(owt='true_owt')
    x = x + noise[n].sel(band=x.band.values).drop_vars('wl', errors='ignore') * xr.DataArray(
        rng.standard_normal(x.shape), dims=x.dims)
    cls, _ = owt.classify(x, owt.convolve_reference(ref, sens).drop_vars('wl'), dim='band')
    rows[n] = {'10 classes': float((cls == hyp).mean()),
               '7 merged classes': float((owt.merge_classes(cls) == owt.merge_classes(hyp)).mean())}
(100 * pd.DataFrame(rows)).round(1)
'''),
        md('''
## 5. OWT map on a Sentinel-2 image
GRS L2A over Lake Chad, averaged to 60 m. The TRISHNA-like classification uses the S2 bands
490/560/665/865 nm with the references convolved with the same S2 bands.'''),
        code(L2A_LOAD),
        code('''
maps = {}
for sens in (S2W, S2T):
    r = owt.convolve_reference(ref, sens)
    x = rrs.sel(wl=r.wl.values, method='nearest').rename(wl='band').assign_coords(band=r.band.values)
    valid = x.notnull().all('band')
    cls, _ = owt.classify(x.where(valid), r.drop_vars('wl'), dim='band')
    maps[sens.name] = cls.where(valid)

labels = list(ref.owt.values)
to_idx = lambda c: xr.apply_ufunc(np.vectorize(lambda v: labels.index(v) if isinstance(v, str) else np.nan), c)
fig, axs = plt.subplots(1, 3, figsize=(17, 7), sharey=True)
cm = plt.get_cmap('tab10', 10)
for ax, (n, c) in zip(axs, maps.items()):
    im = ax.imshow(to_idx(c), cmap=cm, vmin=-.5, vmax=9.5, interpolation='nearest'); ax.set_title(n)
a, b = maps[S2W.name], maps[S2T.name]
agree = (a == b).where(a.notnull())
axs[2].imshow(agree, cmap='RdYlGn', vmin=0, vmax=1, interpolation='nearest')
axs[2].set_title(f'agreement: {100 * float(agree.mean()):.1f} %  (merged classes: '
                 f'{100 * float((owt.merge_classes(a) == owt.merge_classes(b)).where(a.notnull()).mean()):.1f} %)')
cb = fig.colorbar(im, ax=axs[:2], shrink=.6, ticks=range(10)); cb.ax.set_yticklabels(labels)
for ax in axs: ax.set_xticks([]); ax.set_yticks([])
'''),
        code('''
ok = a.notnull().values
pd.crosstab(a.values[ok], b.values[ok], rownames=['S2 (8 bands)'], colnames=['S2-as-TRISHNA'])
'''),
        md('''
## Take-home messages
- With the 4 TRISHNA water bands, about 9 OWT out of 10 are recovered without noise, ~8 with
  requirement noise (S2: 9.5 and 9). The 860 nm band is essential.
- The sub-types 1/2, 3a/3b and 4a/4b are hardly separable: a 7-class TRISHNA OWT set is
  more robust.
- The OWT remain a useful qualitative product (clear / productive / turbid waters) and a
  switch for algorithm selection.

## References
- Bi, S., Hieronymi, M. et al. (2024), optical water types of inland and coastal waters (10 classes), pyOWT.
- Spyrakos, E. et al. (2018), Optical types of inland and coastal waters, L&O 63.'''),
    ]


# =============================================================================================
@register('4_visualisation/wq_products_trishna.ipynb')
def nb06():
    return header(
        'Which water quality products survive with TRISHNA bands?',
        ['A GRS L2A Sentinel-2 product.'],
        ['List the water quality algorithms that TRISHNA bands can feed.',
         'Compute the products from a Sentinel-2 image with all bands and with TRISHNA-like bands.',
         'Quantify the impact of the band differences (adjustment factors, OWT-dependent) and of the TRISHNA noise.'],
        runtime='~2 min') + imports('''
from trishna_tbx import algorithms, owt
''') + data_cells() + [
        md('''
## 1. Algorithm coverage
Each algorithm declares the wavelengths it needs; a sensor covers it if it has a band within
15 nm of each. This avoids the silent "nearest band" substitution (443 -> 485 nm,
705 -> 670 nm) of generic codes.'''),
        code('''
s2 = make_s2()
cov = algorithms.coverage([s2, TRISHNA])
cov.style.map(lambda v: 'background-color: #c8e6c9' if v is True else
              ('background-color: #ffcdd2' if v is False else ''))
'''),
        md('## 2. Products from a Sentinel-2 image'),
        code(L2A_LOAD),
        code('''
S2T = s2.subset(['B2', 'B3', 'B4', 'B8A'], 'S2-as-TRISHNA')
prod_s2 = algorithms.run_all(rrs, s2)
rrs_t = rrs.sel(wl=[490, 560, 665, 865])
prod_t = algorithms.run_all(rrs_t, S2T)
print('S2 only :', sorted(set(prod_s2) - set(prod_t)))
print('both    :', sorted(prod_t))
'''),
        code('''
show = [('SPM_Nechad10', 'SPM (g m$^{-3}$)', 0, 60), ('TUR_Dogliotti15', 'turbidity (FNU)', 0, 40),
        ('Kd490_KD2', 'Kd490 (m$^{-1}$)', 0, 2), ('Chla_OC2', 'Chl-a OC2 (mg m$^{-3}$)', 0, 20),
        ('Chla_index_NDCI', 'NDCI (S2 only)', -0.2, 0.3)]
fig, axs = plt.subplots(1, 5, figsize=(20, 7), sharey=True)
for ax, (v, lab, lo, hi) in zip(axs, show):
    src = prod_t if v in prod_t else prod_s2
    im = ax.imshow(src[v], vmin=lo, vmax=hi, cmap='Spectral_r', interpolation='nearest')
    fig.colorbar(im, ax=ax, orientation='horizontal', pad=.03, label=lab)
    ax.set_xticks([]); ax.set_yticks([])
fig.suptitle(f'{L2A_TILE} {str(date)[:10]}, 60 m')
fig.savefig(THUMB_DIR / 'wq_products_trishna.png', dpi=60)
'''),
        md('''
OC2 is designed for clear waters: in these turbid waters its values are not meaningful (it is
shown because it is the only Chl-a algorithm left with TRISHNA bands). The red-edge NDCI,
sensitive to Chl-a in turbid/productive waters, is lost.

## 3. From "S2 restricted to TRISHNA bands" to TRISHNA
On the same S2 values, covered algorithms give identical products. Two effects make TRISHNA
different: (1) the band responses (adjustment factors, which depend on the water type),
(2) the radiometric noise.'''),
        code('''
from trishna_tbx.emulate import owt_band_adjustment
# OWT of each pixel (TRISHNA-like bands), then OWT-dependent S2 -> TRISHNA factors
rrs_tri, cls, factors = owt_band_adjustment(rrs_t, S2T)
factors.round(3)
'''),
        code('''
prod_tri = algorithms.run_all(rrs_tri, S2T)

# TRISHNA noise (requirement SNR, TOA -> Rrs)
sig = (TRISHNA.subset(['Blue', 'Green', 'Red', 'NIR']).noise_reflectance(30.) / (np.pi * 0.75))
sig = xr.DataArray(sig.values, dims='wl', coords={'wl': rrs_t.wl})
rng = np.random.default_rng(0)
rrs_noisy = rrs_tri + sig * xr.DataArray(rng.standard_normal(rrs_tri.shape), dims=rrs_tri.dims)
prod_noisy = algorithms.run_all(rrs_noisy, S2T)

rows = {}
for v in prod_t:
    a, b, c = prod_t[v], prod_tri[v], prod_noisy[v]
    rel = lambda u: 100 * np.nanmedian(np.abs(u - a) / np.abs(a))
    rows[v] = {'median S2': float(a.median()), 'band adjustment: MdAPD (%)': rel(b),
               '+ noise: MdAPD (%)': rel(c)}
pd.DataFrame(rows).T.round(2)
'''),
        code('''
fig, axs = plt.subplots(1, 3, figsize=(16, 4.8))
for ax, (v, lab, lo, hi) in zip(axs, show[:3]):
    a, c = prod_t[v].values.ravel(), prod_noisy[v].values.ravel()
    ok = np.isfinite(a) & np.isfinite(c)
    ax.hexbin(a[ok], c[ok], gridsize=60, bins='log', extent=(lo, hi, lo, hi), cmap='magma_r')
    ax.plot([lo, hi], [lo, hi], 'k--'); ax.set_xlabel(f'{lab}, S2 bands'); ax.set_ylabel('TRISHNA-like + noise')
fig.tight_layout()
'''),
        md('''
## Take-home messages
- Retrievable with TRISHNA bands: SPM / turbidity (red, NIR), diffuse attenuation (Kd490,
  Kd_PAR), OWT, and blue/green Chl-a in clear waters only.
- Lost: red-edge Chl-a (productive/eutrophic waters), CDOM, QAA-based IOPs and Secchi depth
  as implemented (443 nm); a QAA variant using 485 nm has to be developed.
- The band differences are not negligible: TRISHNA Red (670/60 nm) is wider and redder than
  S2 B4, so in very turbid waters (OWT 5b) it sees ~15 % more signal, and SPM/turbidity
  increase accordingly; blue-band products (OC2, KD2) change by 20-30 %. Algorithm
  coefficients must be re-tuned for the TRISHNA bands (`algorithms.retune_ocx`, simulations).
- Requirement noise adds scatter but, at 60 m, less than the band effect.'''),
    ]


# =============================================================================================
@register('6_trishna_simulations/inversion_information_content.ipynb')
def nb07():
    return header(
        'Semi-analytical inversion with TRISHNA bands: what can be retrieved?',
        ['No image needed (invRrs forward model).'],
        ['Simulate water reflectances from concentrations (Chl-a, CDOM, SPM) with invRrs.',
         'Integrate them over TRISHNA and Sentinel-2 bands.',
         'Map the information content (expected retrieval uncertainty) of both sensors.',
         'Invert noisy synthetic TRISHNA spectra.'],
        runtime='~3 min') + imports('''
import invRrs.inverse as inv
from invRrs import utils as u
from trishna_tbx.inversion import trishna_forward_model, information_content
''') + [
        md('''
## 1. Forward model
`MultispectralSimple` (invRrs): Chl-a (phytoplankton absorption/backscattering), CDOM
absorption at 440 nm, SPM (Mie, normalised) and a spectrally flat offset (residual glint /
atmospheric correction error). Rrs is computed at 1 nm and integrated over each band,
weighted by the downwelling irradiance.'''),
        code('''
aphy, bbphy = u.data().load_iop_phyto('Bricaud1995')
SZA = 30.
models = {'TRISHNA': trishna_forward_model('MultispectralSimple', aphy_star=aphy, bbphy_star=bbphy, sza=SZA),
          'S2': inv.MultispectralSimple(aphy, bbphy, SZA, sat='s2a', band_idx=[0, 1, 2, 3, 4, 5, 6, 8])}
hyper = inv.MultispectralSimple(aphy, bbphy, SZA, wl_=np.arange(400., 901., 2.))
for n, m in models.items():
    print(n, np.round(m.wl), m.param_names)
'''),
        code('''
cases = {'clear': [0.5, 0.05, 0.2, 0.], 'productive': [20., 0.5, 2., 0.], 'turbid': [3., 0.5, 30., 0.],
         'CDOM-rich': [2., 2., 2., 0.]}
fig, axs = plt.subplots(1, 4, figsize=(18, 4))
for ax, (n, x) in zip(axs, cases.items()):
    ax.plot(hyper.wl, hyper.forward_model(np.array(x), output='Rrs'), c='0.5', label='hyperspectral')
    for (s, m), mk in zip(models.items(), ['o', 's']):
        ax.plot(m.wl, m.forward_model(np.array(x), output='Rrs'), mk, ms=7, label=s, mfc='none' if mk == 's' else None)
    ax.set_title(f'{n}: chl={x[0]}, a_cdom={x[1]}, spm={x[2]}'); ax.set_xlabel('wavelength (nm)')
axs[0].set_ylabel('Rrs (sr$^{-1}$)'); axs[0].legend()
fig.tight_layout()
fig.savefig(THUMB_DIR / 'inversion_information_content.png', dpi=60)
'''),
        md('''
## 2. Information content
Linearised (Gauss-Newton) posterior standard deviation of each unknown for an Rrs uncertainty
of 5 % + 2e-4 sr-1 per band. A relative standard deviation > 1 means the parameter is not
retrievable at this state. With 4 bands and 4 unknowns, TRISHNA is exactly determined.'''),
        code('''
%%time
chl = np.logspace(-1, 2, 9)
spm = np.logspace(-1, 2, 9)
A_CDOM = 0.3
out = {n: np.full((2, chl.size, spm.size), np.nan) for n in models}
for n, m in models.items():
    for i, c in enumerate(chl):
        for j, s in enumerate(spm):
            x = [c, A_CDOM, s, 0.]
            r = np.asarray(m.forward_model(np.array(x), output='Rrs'))
            ic = information_content(m, x, 0.05 * r + 2e-4)
            out[n][:, i, j] = ic['std'].values[[0, 2]] / np.array([c, s])
'''),
        code('''
fig, axs = plt.subplots(2, 2, figsize=(12, 10), sharex=True, sharey=True)
for col, n in enumerate(models):
    for row, p in enumerate(['Chl-a', 'SPM']):
        ax = axs[row, col]
        im = ax.pcolormesh(spm, chl, out[n][row], norm=plt.matplotlib.colors.LogNorm(0.05, 10), cmap='RdYlGn_r')
        cs = ax.contour(spm, chl, out[n][row], levels=[0.3, 1.], colors='k', linewidths=1)
        ax.clabel(cs, fmt='%.1f'); ax.set_xscale('log'); ax.set_yscale('log')
        ax.set_title(f'{n}: relative std of {p} (a_cdom={A_CDOM})')
for ax in axs[-1]: ax.set_xlabel('SPM (g m$^{-3}$)')
for ax in axs[:, 0]: ax.set_ylabel('Chl-a (mg m$^{-3}$)')
fig.colorbar(im, ax=axs, shrink=.6, label='relative standard deviation')
'''),
        md('''
## 3. Inversion of noisy synthetic spectra
Monte-Carlo: 100 noise realisations of a productive, moderately turbid case.'''),
        code('''
%%time
truth = np.array([10., 0.5, 5., 0.])
rng = np.random.default_rng(0)
res = {}
for n, m in models.items():
    r0 = np.asarray(m.forward_model(truth, output='Rrs'))
    sigma = 0.05 * r0 + 2e-4
    sols = []
    for k in range(100):
        y = r0 + sigma * rng.standard_normal(r0.size)
        sol = m.call_solver_scipy(y, sigma=sigma)
        sols.append(sol.x)
    res[n] = pd.DataFrame(sols, columns=m.param_names)
fig, axs = plt.subplots(1, 3, figsize=(15, 4))
for ax, (p, t) in zip(axs, zip(['chl', 'a_cdom', 'spm_norm'], truth)):
    ax.boxplot([res[n][p] for n in models], labels=list(models)); ax.axhline(t, c='r', ls='--')
    ax.set_yscale('log'); ax.set_title(p)
pd.concat({n: (100 * (r.iloc[:, :3].quantile([0.16, 0.5, 0.84]) / truth[:3] - 1)).round(0)
           for n, r in res.items()}, axis=1)
'''),
        md('''
## Take-home messages
- With TRISHNA bands, Chl-a and CDOM are strongly correlated (no 443 nm, no red edge): Chl-a
  is only weakly constrained, except at high concentrations.
- SPM remains retrievable from low to high turbidity, slightly less precisely than with S2.
- Priors are needed for Chl-a: OWT-based constraints, regional CDOM-salinity/temperature
  relationships, or the temporal information (TRISHNA revisit and thermal bands).'''),
    ]


# =============================================================================================
@register('4_visualisation/datacube_trishna.ipynb')
def nb08():
    return header(
        'Time series of TRISHNA-like water products (datacube)',
        ['A time series of GRS L2A Sentinel-2 products on one tile.'],
        ['Build a 60 m datacube of TRISHNA-like water products from a Sentinel-2 L2A time series.',
         'Compare the time series of products available with TRISHNA and with all S2 bands.',
         'Follow the optical water types through time.'],
        resources='2 cores, 8 GB', runtime='~5 min') + imports('''
from trishna_tbx import algorithms, owt
from trishna_tbx.io import list_l2a, open_grs_l2a, water_mask, to_60m
from trishna_tbx.emulate import owt_band_adjustment
''') + data_cells() + [
        md('''
## 1. Build the datacube
For each date: water pixels, 60 m averaging, then products from all S2 bands and from
TRISHNA-like Rrs (S2 490/560/665/865 nm with OWT-dependent S2 -> TRISHNA band adjustment and
requirement noise), and OWT with both band sets.'''),
        code('''
files = list_l2a(L2A_DIR, L2A_TILE)
print(len(files), 'dates:', str(files[0][0])[:10], '->', str(files[-1][0])[:10])
s2 = make_s2()
S2W = s2.subset(['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8A'], 'S2')
S2T = s2.subset(['B2', 'B3', 'B4', 'B8A'], 'S2-as-TRISHNA')
ref, _ = owt.load_reference('Bi2024')
sig = TRISHNA.subset(['Blue', 'Green', 'Red', 'NIR']).noise_reflectance(30.) / (np.pi * 0.75)
rng = np.random.default_rng(0)

def classify(rrs_, sens):
    r = owt.convolve_reference(ref, sens)
    x = rrs_.sel(wl=r.wl.values, method='nearest').rename(wl='band').assign_coords(band=r.band.values)
    v = x.notnull().all('band')
    c, _ = owt.classify(x.where(v), r.drop_vars('wl'), dim='band')
    return owt.merge_classes(c.where(v))

cubes = []
for t, f in files:
    ds = open_grs_l2a(f, L2A_WINDOW)
    rrs = to_60m(ds.Rrs.where(water_mask(ds))).load()
    rrs = rrs.where((rrs.sel(wl=[490, 560, 665]) > 0).all('wl'))
    rt, _, _ = owt_band_adjustment(rrs.sel(wl=[490, 560, 665, 865]), S2T)
    rt = rt + xr.DataArray(sig.values, dims='wl', coords={'wl': rt.wl}) * xr.DataArray(
        rng.standard_normal(rt.shape), dims=rt.dims)
    pt = algorithms.run_all(rt, S2T)
    ps = algorithms.run_all(rrs, s2)
    cube = xr.Dataset({
        'SPM_trishna': pt.SPM_Nechad10, 'TUR_trishna': pt.TUR_Dogliotti15, 'Kd490_trishna': pt.Kd490_KD2,
        'SPM_s2': ps.SPM_Nechad10, 'TUR_s2': ps.TUR_Dogliotti15, 'Kd490_s2': ps.Kd490_KD2,
        'NDCI_s2': ps.Chla_index_NDCI,
        'owt_trishna': classify(rrs, S2T).astype(str), 'owt_s2': classify(rrs, S2W).astype(str)})
    cubes.append(cube.expand_dims(time=[t]))
cube = xr.concat(cubes, 'time')
cube
'''),
        code('''
f = OUTPUT_DIR / f'datacube_trishna_{L2A_TILE}.nc'
cube.to_netcdf(f); print(f)
'''),
        md('## 2. Time series over the lake'),
        code('''
valid = cube.SPM_s2.notnull().sum(['x', 'y'])
keep = valid > 0.2 * valid.max()   # dates with enough water pixels (clouds)
c = cube.sel(time=keep)
fig, axs = plt.subplots(4, 1, figsize=(12, 12), sharex=True)
for ax, (v, lab) in zip(axs, [('SPM', 'SPM (g m$^{-3}$)'), ('TUR', 'turbidity (FNU)'), ('Kd490', 'Kd490 (m$^{-1}$)')]):
    for src, col in [('s2', 'tab:blue'), ('trishna', 'tab:red')]:
        q = c[f'{v}_{src}'].quantile([0.1, 0.5, 0.9], dim=['x', 'y'])
        ax.fill_between(c.time, q.sel(quantile=.1), q.sel(quantile=.9), color=col, alpha=.15)
        ax.plot(c.time, q.sel(quantile=.5), 'o-', c=col, label='S2' if src == 's2' else 'TRISHNA-like + noise')
    ax.set_ylabel(lab); ax.legend(fontsize=8)
q = c.NDCI_s2.quantile([0.1, 0.5, 0.9], dim=['x', 'y'])
axs[3].fill_between(c.time, q.sel(quantile=.1), q.sel(quantile=.9), color='tab:green', alpha=.15)
axs[3].plot(c.time, q.sel(quantile=.5), 'o-', c='tab:green'); axs[3].set_ylabel('NDCI (S2 only)')
fig.tight_layout()
fig.savefig(THUMB_DIR / 'datacube_trishna.png', dpi=60)
'''),
        md('## 3. Optical water types through time'),
        code('''
def occurrence(v):
    df = c[v].to_dataframe(name='owt').reset_index()
    df = df[df.owt.isin(sorted(set(owt.TRISHNA_MERGE.values())))]
    return pd.crosstab(df.time, df.owt, normalize='index')

fig, axs = plt.subplots(1, 2, figsize=(16, 4.5), sharey=True)
for ax, (v, t) in zip(axs, [('owt_s2', 'S2 (8 bands)'), ('owt_trishna', 'TRISHNA-like')]):
    o = occurrence(v); o.index = o.index.strftime('%m-%d')
    o.plot.bar(stacked=True, ax=ax, width=.9, colormap='tab10', legend=ax is axs[1]); ax.set_title(t)
axs[1].legend(bbox_to_anchor=(1, 1), title='OWT (merged)'); axs[0].set_ylabel('fraction of water pixels')
'''),
        md('''
## 4. Maps
SPM from TRISHNA-like bands for a few dates.'''),
        code('''
sel = c.isel(time=np.linspace(0, c.sizes['time'] - 1, 6).astype(int))
g = sel.SPM_trishna.plot(col='time', col_wrap=6, vmin=0, vmax=60, cmap='Spectral_r', figsize=(18, 6),
                         cbar_kwargs={'label': 'SPM (g m$^{-3}$)'})
for ax in g.axs.ravel(): ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel(''); ax.set_ylabel('')
'''),
        md('''
## Take-home messages
- The temporal dynamics of SPM, turbidity and Kd490 are preserved with TRISHNA-like bands,
  but with a systematic offset due to the band differences (here SPM/turbidity +7 to +16 %,
  Kd490 +20 to +40 % with S2 coefficients): coefficients must be re-tuned for TRISHNA before
  merging S2 and TRISHNA time series.
- Phytoplankton dynamics seen by the red-edge index (NDCI) cannot be followed with TRISHNA.
- TRISHNA revisit (3 acquisitions every 8 days, day and night for the TIR) will densify these
  series, but with variable viewing angles (up to 38 deg): angular effects and sunglint have
  to be handled (GRS).'''),
    ]


# =============================================================================================
HYP_DATA = """import os
from pathlib import Path

# ---- edit here (or set the environment variables) ------------------------------------
HYP_FILE = os.environ.get('TRISHNA_TBX_HYPERNETS', '/DATA/hypernets/BEFR/BEFR_hypernets.nc')  # Berre lagoon
L1C_SUBSET_DIR = os.environ.get('TRISHNA_TBX_L1C_SUBSETS',
                                '/data/satellite/Sentinel-2/subset/hypernets/L1C/BEFR')  # S2 L1C around the site
CAMS_DIR = os.environ.get('TRISHNA_TBX_CAMS_DIR', '/data/cams/world')   # cams_forecast_YYYY-MM.nc
MAX_DT = '1h'           # maximum time difference HYPERNETS - Sentinel-2
HALF_WINDOW = 1         # 3x3 pixels at 60 m around the site
MAX_MATCHUPS = None     # None: all (about 1 h of GRS processing the first time, then cached)
# optional reference: archived standard GRS L2A (all bands, 20 m) of the same acquisitions (None to skip)
L2A_ARCHIVE_DIR = os.environ.get('TRISHNA_TBX_L2A_ARCHIVE', '/data/satellite/Sentinel-2/subset/hypernets/L2A_v221/BEFR')
OUTPUT_DIR = Path(os.environ.get('TRISHNA_TBX_OUTPUT', './output')).resolve()
# ---------------------------------------------------------------------------------------
CACHE_DIR = OUTPUT_DIR / 'hypernets_matchups'
CACHE_DIR.mkdir(parents=True, exist_ok=True)
"""


@register('5_validation/validation_hypernets_trishna.ipynb')
def nb09():
    return header(
        'Validation with HYPERNETS: Sentinel-2 vs TRISHNA bands',
        ['HYPERNETS water reflectance (L2A) of a site (here BEFR, Berre lagoon, from the WATERHYPERNET portal).',
         'Sentinel-2 L1C subsets around the site (NetCDF, GRSdriver format) and the CAMS files.',
         'The GRS look-up tables.'],
        ['Convolve the hyperspectral HYPERNETS reflectance with the TRISHNA and Sentinel-2 bands.',
         'Measure the real S2 -> TRISHNA band differences and the OWT classification on in situ spectra.',
         'Validate GRS run with TRISHNA-like bands and with all S2 bands on the same matchups.',
         'Validate the water quality products available with TRISHNA bands.'],
        resources='4 cores, 8 GB', runtime='~1 h the first time (GRS on all matchups), then < 2 min') + imports("""
from trishna_tbx import algorithms, owt
from trishna_tbx.validation import (load_hypernets, convolve, find_matchups, process_matchup, metrics,
                                    find_l2a, extract_l2a_window)
from trishna_tbx.emulate import owt_band_adjustment
from pyproj import Transformer
""") + [md('## Data\n> HYPERNETS data: https://waterhypernet.org (account needed). The L1C subsets can be '
           'produced with GRSdriver (`Sentinel2Driver(..., subset=...)` then `prod.to_netcdf`).'),
            code(HYP_DATA)] + [
        md("""
> This tutorial extends *Validation L2A* (Validation section) to the TRISHNA bands.

## 1. HYPERNETS at the Berre lagoon
HYPERNETS measures the water-leaving reflectance every 30 min from a fixed platform
(HYPSTAR radiometer, 400-900 nm in the distributed L2A). `Rrs` is the reflectance divided
by pi, with the NIR similarity correction (Ruddick et al. 2006). Basic QC: quality flag 0,
solar zenith angle < 70 deg."""),
        code("""
hyp = load_hypernets(HYP_FILE)
site = hyp.attrs['site_id']
lon, lat = hyp.attrs['site_longitude'], hyp.attrs['site_latitude']
site_xy = Transformer.from_crs(4326, 32631, always_xy=True).transform(lon, lat)
print(f'{site}: {hyp.sizes["time"]} records after QC, {str(hyp.time.values[0])[:10]} -> {str(hyp.time.values[-1])[:10]}')
print(f'site: {lat:.4f} N, {lon:.4f} E  ->  UTM 31N {site_xy[0]:.0f}, {site_xy[1]:.0f}')
"""),
        code("""
fig, axs = plt.subplots(1, 2, figsize=(16, 4.5), gridspec_kw={'width_ratios': [2, 1]})
for w, c in [(490, 'tab:blue'), (560, 'tab:green'), (665, 'tab:red'), (860, 'k')]:
    hyp.Rrs.sel(wl=w, method='nearest').plot(ax=axs[0], ls='', marker='.', ms=2, c=c, label=f'{w} nm')
axs[0].legend(markerscale=6); axs[0].set_title(''); axs[0].set_ylabel('Rrs (sr$^{-1}$)')
sub = hyp.Rrs.isel(time=slice(None, None, 40))
for i in range(sub.sizes['time']):
    sub.isel(time=i).plot(ax=axs[1], lw=.6, alpha=.6)
for b in TRISHNA.subset(['Blue', 'Green', 'Red', 'NIR']).bands:
    axs[1].axvspan(b.center - b.fwhm / 2, b.center + b.fwhm / 2, color='0.85', zorder=0)
axs[1].set_title('spectra (grey: TRISHNA bands)'); axs[1].set_ylabel('')
fig.tight_layout()
fig.savefig(THUMB_DIR / 'validation_hypernets_trishna.png', dpi=60)
"""),
        md("""
## 2. HYPERNETS seen by TRISHNA and Sentinel-2
The hyperspectral Rrs is band-averaged (weighted by the solar irradiance) with each sensor
response. The ratio TRISHNA / S2 is the **real** band adjustment factor for this site, to be
compared with the factors derived from the Bi et al. (2024) OWT spectra. The TRISHNA WV band
(910 nm) is outside the 400-900 nm range of the data."""),
        code("""
s2 = make_s2()
S2W = s2.subset(['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8A'], 'S2')
S2T = s2.subset(['B2', 'B3', 'B4', 'B8A'], 'S2-as-TRISHNA')
TRW = TRISHNA.subset(['Blue', 'Green', 'Red', 'NIR'], 'TRISHNA')
h_s2 = convolve(hyp.Rrs, S2W)
h_tri = convolve(hyp.Rrs, TRW)
pairs = [('Blue', 'B2'), ('Green', 'B3'), ('Red', 'B4'), ('NIR', 'B8A')]
ratio = pd.DataFrame({t: (h_tri.sel(band=t) / h_s2.sel(band=b)).to_pandas() for t, b in pairs})
ratio.describe(percentiles=[.05, .5, .95]).T[['5%', '50%', '95%']].round(3)
"""),
        code("""
# OWT of the in situ spectra (hyperspectral reference classification)
ref, _ = owt.load_reference('Bi2024')
hyp_owt, _ = owt.classify(hyp.Rrs.interp(wl=ref.wl.sel(wl=slice(400, 900))),
                          ref.sel(wl=slice(400, 900)), dim='wl')
# factors predicted from the OWT mean spectra (as used to emulate TRISHNA)
_, _, f_table = owt_band_adjustment(h_s2.sel(band=['B2', 'B3', 'B4', 'B8A'])
                                    .drop_vars('wl').assign_coords(band=[490, 560, 665, 865]).rename(band='wl'), S2T)
obs = ratio.assign(owt=hyp_owt.values).groupby('owt').median()
obs.columns = [490, 560, 665, 865]
cmp = pd.concat({'HYPERNETS (median)': obs, 'OWT model': f_table.loc[obs.index]}, axis=1)
cmp.insert(0, 'N', ratio.assign(owt=hyp_owt.values).groupby('owt').size())
cmp.round(3)
"""),
        md("""
## 3. Optical water types from in situ spectra
Real-data version of the benchmark of *Optical water types with TRISHNA bands*: the class
obtained with the full spectrum is the reference."""),
        code("""
rows = {}
for sens, conv in [(S2W, h_s2), (TRW, h_tri)]:
    r = owt.convolve_reference(ref, sens)
    cls, _ = owt.classify(conv, r.drop_vars('wl'), dim='band')
    rows[sens.name] = {'10 classes': float((cls.values == hyp_owt.values).mean()),
                       '7 merged classes': float((owt.merge_classes(cls).values == owt.merge_classes(hyp_owt).values).mean())}
print('OWT occurrence (hyperspectral):', hyp_owt.to_series().value_counts().to_dict())
(100 * pd.DataFrame(rows)).round(1)
"""),
        md("""
## 4. Matchups with Sentinel-2
S2 L1C subsets with a HYPERNETS record within `MAX_DT`. GRS is run at 60 m on each subset
twice: with the TRISHNA-like bands (aerosol/glint from 865/1610 nm) and with all the S2
bands (1610/2190 nm), with the same GRS version and settings. The median of the valid water
pixels in a 3x3 window (180 m) around the platform is kept. Results are cached in
`CACHE_DIR`."""),
        code("""
mu = find_matchups(hyp.time.values, L1C_SUBSET_DIR, MAX_DT)
if MAX_MATCHUPS:
    mu = mu.iloc[:MAX_MATCHUPS]
print(len(mu), 'matchups within', MAX_DT)
"""),
        code("""
%%time
res = []
for i, row in mu.iterrows():
    r = process_matchup(row.l1c, CAMS_DIR, CACHE_DIR, site_xy, half=HALF_WINDOW)
    if r is not None:
        res.append(r.expand_dims(matchup=[i]))
sat = xr.concat(res, 'matchup', join='outer')
print(f'{sat.sizes["matchup"]} matchups processed')
"""),
        md("""
Reference: the archived standard GRS products (GRS v2.2.1, all bands, native 20 m, 180 m window)
of the same acquisitions, to check that the 60 m processing used here does not degrade the
results."""),
        code("""
if L2A_ARCHIVE_DIR and Path(L2A_ARCHIVE_DIR).exists():
    arch = []
    for i in sat.matchup.values:
        f = find_l2a(L2A_ARCHIVE_DIR, mu.loc[i, 's2_time'])
        if f:
            arch.append(extract_l2a_window(f, site_xy, half_m=90, name='S2 archive (20 m)').expand_dims(matchup=[i]))
    print(f'{len(arch)} archived products found')
    sat = xr.concat([sat.reset_coords(drop=True), xr.concat(arch, 'matchup').reset_coords(drop=True)], 'bandset', join='outer')
BANDSETS = list(sat.bandset.values)
print(BANDSETS)
"""),
        code("""
# in situ values at the matchup times, quality screening of the satellite windows
hm = mu.loc[sat.matchup.values]
ins_s2 = h_s2.sel(time=hm.hyp_time.values).assign_coords(time=sat.matchup.values).rename(time='matchup')
ins_tri = h_tri.sel(time=hm.hyp_time.values).assign_coords(time=sat.matchup.values).rename(time='matchup')
cv560 = sat.Rrs_std.sel(wl=560) / sat.Rrs.sel(wl=560)
good = (sat.n_valid >= 5).all('bandset') & (cv560 < 0.2).all('bandset') & sat.Rrs.sel(wl=560).notnull().all('bandset')
good = good.drop_vars('wl', errors='ignore')
print(f'{int(good.sum())} matchups kept (>= 5 valid water pixels, CV(560) < 20 %)')
sat_g, s2_g, tri_g = sat.sel(matchup=good), ins_s2.sel(matchup=good), ins_tri.sel(matchup=good)
"""),
        md("""
## 5. Rrs validation
All GRS runs use the same S2 measurements: they are compared with HYPERNETS convolved with
the **S2** responses. Differences between the `S2` and `TRISHNA` rows come only from the
atmospheric correction band set (aerosol/glint bands); the last row is the archived reference."""),
        code("""
wls = [(490, 'B2'), (560, 'B3'), (665, 'B4'), (865, 'B8A')]
fig, axs = plt.subplots(len(BANDSETS), 4, figsize=(18, 4.5 * len(BANDSETS)), squeeze=False)
stats = {}
for j, (w, b) in enumerate(wls):
    x = s2_g.sel(band=b).values
    lim = [min(0, np.nanmin(x)), 1.3 * np.nanmax(x)]
    for i, bs in enumerate(BANDSETS):
        y = sat_g.Rrs.sel(bandset=bs, wl=w).values
        ax = axs[i, j]
        ax.errorbar(x, y, yerr=sat_g.Rrs_std.sel(bandset=bs, wl=w).values, fmt='o', ms=4, alpha=.7, lw=.6)
        ax.plot(lim, lim, 'k--', lw=1); ax.set_xlim(lim); ax.set_ylim(lim)
        ax.set_title(f'{w} nm, GRS: {bs}'); ax.set_xlabel('HYPERNETS (S2 bands)')
        stats[(bs, w)] = metrics(x, y)
for ax in axs[:, 0]:
    ax.set_ylabel('GRS Rrs (sr$^{-1}$)')
fig.tight_layout()
pd.DataFrame(stats).T.round(5)
"""),
        code("""
fig, axs = plt.subplots(1, 2, figsize=(12, 4.5))
axs[0].scatter(sat_g.aot550.sel(bandset='S2'), sat_g.aot550.sel(bandset='TRISHNA'), s=15)
m = float(sat_g.aot550.max()); axs[0].plot([0, m], [0, m], 'k--')
axs[0].set_xlabel('AOT(550), all S2 bands'); axs[0].set_ylabel('AOT(550), TRISHNA bands')
err = {bs: (sat_g.Rrs.sel(bandset=bs, wl=665) - s2_g.sel(band='B4')).values for bs in ['S2', 'TRISHNA']}
axs[1].scatter(sat_g.aot550.sel(bandset='S2') - sat_g.aot550.sel(bandset='TRISHNA'), err['S2'] - err['TRISHNA'], s=15)
axs[1].axhline(0, c='k', lw=.5); axs[1].axvline(0, c='k', lw=.5)
axs[1].set_xlabel('AOT difference (S2 - TRISHNA bands)'); axs[1].set_ylabel('Rrs(665) error difference (S2 - TRISHNA bands)')
fig.tight_layout()
"""),
        md("""
### What a TRISHNA validation would see
With real TRISHNA data, HYPERNETS would be convolved with the TRISHNA responses. On the
emulated data, the GRS-TRISHNA Rrs (S2 bands) is converted to TRISHNA bands with the OWT-dependent
factors, then compared with HYPERNETS convolved with the TRISHNA responses."""),
        code("""
rt = sat_g.Rrs.sel(bandset='TRISHNA', wl=[490, 560, 665, 865]).transpose('wl', 'matchup')
rt_adj, _, _ = owt_band_adjustment(rt, S2T)
stats_t = {}
for (w, _), t in zip(wls, ['Blue', 'Green', 'Red', 'NIR']):
    stats_t[(t, 'S2 bands (no adjustment)')] = metrics(tri_g.sel(band=t).values, rt.sel(wl=w).values)
    stats_t[(t, 'OWT band adjustment')] = metrics(tri_g.sel(band=t).values, rt_adj.sel(wl=w).values)
pd.DataFrame(stats_t).T[['N', 'bias', 'MdAPD (%)', 'r']].round(4)
"""),
        md("""
## 6. Water quality products
Products computed from HYPERNETS (S2 bands) and from GRS with both band sets, for the
algorithms available with TRISHNA bands."""),
        code("""
def products(rrs_band, names):
    r = rrs_band.sel(band=names).drop_vars('wl').assign_coords(band=[490, 560, 665, 865]).rename(band='wl')
    return algorithms.run_all(r.transpose('wl', ...), S2T)

p_ins = products(s2_g, ['B2', 'B3', 'B4', 'B8A'])
p_sat = {bs: algorithms.run_all(sat_g.Rrs.sel(bandset=bs, wl=[490, 560, 665, 865]).transpose('wl', ...), S2T)
         for bs in BANDSETS}
rows = {}
for v in ['SPM_Nechad10', 'TUR_Dogliotti15', 'Kd490_KD2', 'Kd_par_RD22', 'Chla_OC2']:
    for bs in BANDSETS:
        rows[(v, f'GRS {bs}')] = metrics(p_ins[v].values, p_sat[bs][v].values)
pd.DataFrame(rows).T[['N', 'bias', 'MdAPD (%)', 'slope', 'r']].round(3)
"""),
        md("""
## Take-home messages
With the default data (HYPERNETS BEFR, Berre lagoon, 2021-2023: 80 S2 matchups within 1 h,
64 kept after screening):
- **In situ band differences**: TRISHNA / S2 ratios of the in situ Rrs are ~0.92 (Blue), ~0.94
  (Green), ~1.00 (Red), ~1.01 (NIR). The OWT-based model used to emulate TRISHNA reproduces Blue and
  Red within ~1 %, but is off by 3-5 % for Green and NIR: use site data such as HYPERNETS when available.
- **OWT on real spectra**: TRISHNA bands give the hyperspectral class in ~93 % of the cases (S2: ~94 %),
  much better than the synthetic benchmark: real spectra at one site are less variable than the
  synthetic members.
- **Atmospheric correction**: with the same GRS version and settings, the TRISHNA band set
  (aerosol/glint from 865/1610 nm) has a larger negative bias than the S2 band set
  (1610/2190 nm): about -0.002 sr-1 in the visible, MdAPD 35-60 % vs 23-40 %, with similar
  correlations. The AOT figure shows why: the matchups where the 865/1610 pair retrieves a
  higher AOT than the SWIR pair are those with the most negative Rrs error. In these moderately
  turbid waters the 865 nm water signal is not zero and is partly taken as aerosol
  (over-correction): a TRISHNA processing needs a NIR water-signal correction (similarity
  spectrum, iterative scheme) or a stronger weight on 1610 nm.
- **Processing setup**: the archived standard GRS products (20 m, GRS v2.2.1) perform about
  twice better (MdAPD 15-22 %) than the 60 m re-processing used here for both band sets: the
  TRISHNA vs S2 comparison is fair, but absolute errors are pessimistic.
- **Band adjustment** changes the TRISHNA validation statistics by only a few %: the atmospheric
  correction error dominates.
- **Products**: SPM / turbidity inherit the visible/red bias (MdAPD ~63 % with TRISHNA bands,
  ~41 % with S2 bands, ~23 % archive); Kd_PAR (Roy & Das) is robust (~7 %); Kd490 and OC2 are poor
  at this site with all band sets.
- The 865 nm Rrs is not validated (no skill for any band set) and the TRISHNA WV band (910 nm) is
  outside the distributed HYPERNETS range.

## References
- HYPERNETS / WATERHYPERNET: https://waterhypernet.org, processor https://hypernets-processor.readthedocs.io
- Ruddick, K. et al. (2006), Seaborne measurements of near infrared water-leaving reflectance: the similarity spectrum, L&O 51.
- Harmel, T. et al. (2018), Sunglint correction of the MSI-SENTINEL-2 imagery over inland and sea waters from SWIR bands, RSE 204."""),
    ]



def main():
    execute = '--execute' in sys.argv
    only = [a for a in sys.argv[1:] if not a.startswith('--')]
    for rel, build in NOTEBOOKS.items():
        if only and not any(o in rel for o in only):
            continue
        path = HERE / rel
        path.parent.mkdir(parents=True, exist_ok=True)
        write(path, build())
        print('written', rel)
        if execute:
            extra = ['--allow-errors'] if '--allow-errors' in sys.argv else []
            r = subprocess.run(['jupyter', 'nbconvert', '--to', 'notebook', '--execute', '--inplace',
                                '--ExecutePreprocessor.timeout=-1', *extra, path.name], cwd=path.parent,
                               capture_output=True, text=True)
            print('  executed' if r.returncode == 0 else '  FAILED\n' + r.stderr[-3000:])


if __name__ == '__main__':
    main()
