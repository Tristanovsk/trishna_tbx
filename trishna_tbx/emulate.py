"""Emulate TRISHNA-like L1C / L2A products from Sentinel-2.

Two levels of emulation:

1. **band restriction** (``restrict_bands``): keep only the S2 bands that have a TRISHNA
   counterpart (B2, B3, B4, B8A, B9, B10, B11), degrade to 60 m. This is what the
   CNES proposal asks for ("prendre une image S2 et retirer les bandes non couvertes").
2. **spectral band adjustment** (``SBAF``): correct the S2 band values for the different
   centre/width of the TRISHNA bands (S2 B3 is 36 nm wide vs 70 nm for TRISHNA Green,
   B9 is centred at 945 nm vs 910 nm for the TRISHNA WV band...). Factors are fitted on
   hyperspectral libraries (OWT mean spectra, HYPERNETS, radcalnet_oc TOA simulations).

Optionally add TRISHNA radiometric noise (``add_noise``) from the SNR spec.

GRS integration notes (grs/product.py, grs/mask.py):
- aerosol/glint use the *last two processed bands* -> with the TRISHNA set this is
  B8A (865) + B11 (1610): consistent with TRISHNA (NIR + SWIR 1610).
- the ``high_swir`` flag requires an exact 2190 nm band: TRISHNA has none.
  ``grs_band_set(strict=False)`` keeps B12 loaded for masking only (not processed);
  ``strict=True`` drops it and requires the grs patch described in ROADMAP.md.
- water vapour comes from CAMS ``tcwv``: use ``trishna_tbx.water_vapor`` to replace it by
  an image-based 860/910 (TRISHNA) or 865/945 (S2) retrieval.
"""
import numpy as np
import xarray as xr

from .sensor import S2_TO_TRISHNA, TRISHNA, make_s2, solar_irradiance

#: Order of the Sentinel-2 bands in GRSdriver (``BAND_ID``).
S2_BAND_ORDER = ['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B8A', 'B9', 'B10', 'B11', 'B12']
TRISHNA_S2_IDX = [S2_BAND_ORDER.index(b) for b in S2_TO_TRISHNA]  # [1, 2, 3, 8, 9, 10, 11]
#: Nominal wavelength (nm) of each Sentinel-2 band in the GRS products.
S2_NOMINAL_WL = dict(B1=443, B2=490, B3=560, B4=665, B5=705, B6=740, B7=783, B8=842, B8A=865,
                     B9=945, B10=1375, B11=1610, B12=2190)


def grs_band_set(strict=False):
    """Sentinel-2 band indices of a TRISHNA-like processing with GRSdriver and grs.

    The processed bands are the S2 proxies of the TRISHNA water bands (B2, B3, B4, B8A, B11):
    grs retrieves the aerosol and the sunglint from the last two, 865 and 1610 nm, as TRISHNA
    would. B9 and B10 are loaded for the water vapour and the cirrus only.

    :param strict: if False, B12 (2190 nm) is also loaded, for the ``high_swir`` flag of grs
        only (not processed); if True, B12 is dropped (needs a grs version without this flag)
    :return: (``band_idx``, ``band_tbp_idx``), indices in the GRSdriver band order
        :data:`S2_BAND_ORDER`, for ``Sentinel2Driver(band_idx=..., band_tbp_idx=...)``
    """
    tbp = [S2_BAND_ORDER.index(b) for b in ('B2', 'B3', 'B4', 'B8A', 'B11')]
    load = list(TRISHNA_S2_IDX)
    if not strict:
        load = load + [S2_BAND_ORDER.index('B12')]  # masking only (grs high_swir flag)
    return sorted(load), tbp


def open_s2_l1c_as_trishna(safe, resolution=60, strict=False, **kw):
    """Open a Sentinel-2 L1C SAFE restricted to the TRISHNA-like bands (GRSdriver).

    :param safe: path of the ``.SAFE`` folder
    :param resolution: pixel size of the output (m)
    :param strict: see :func:`grs_band_set`
    :param kw: other arguments of ``GRSdriver.Sentinel2Driver`` (e.g. ``subset``)
    :return: L1C :class:`xarray.Dataset` ready for ``grs.Process().execute(prod, ...)``, with the
        attribute ``emulated_sensor``
    """
    from GRSdriver import Sentinel2Driver
    band_idx, tbp = grs_band_set(strict)
    drv = Sentinel2Driver(safe, band_idx=band_idx, band_tbp_idx=tbp, resolution=resolution,
                          resolution_angle=60, **kw)
    drv.load_product()
    prod = drv.prod
    prod.attrs['emulated_sensor'] = 'TRISHNA (S2 band restriction)'
    return prod


def restrict_bands(ds, var='Rrs', wl_dim='wl', keep=None, resolution=None, coarsen_from=None):
    """Restrict a Sentinel-2 L2A dataset to the TRISHNA-like bands.

    :param ds: L2A :class:`xarray.Dataset` (GRS-like: ``var`` on ``wl``, ``y``, ``x``)
    :param var: name of the reflectance variable
    :param wl_dim: name of the wavelength dimension
    :param keep: nominal S2 wavelengths to keep (nm); default: the S2 proxies of the TRISHNA water
        bands, 490, 560, 665, 865 and 1610 nm
    :param resolution: output pixel size (m), with ``coarsen_from``
    :param coarsen_from: native pixel size (m): the image is block-averaged from ``coarsen_from``
        to ``resolution`` (e.g. 20 -> 60 m)
    :return: :class:`xarray.Dataset` with the TRISHNA band names as coordinate ``band``
    """
    keep = keep or [S2_NOMINAL_WL[b] for b in ('B2', 'B3', 'B4', 'B8A', 'B11')]
    out = ds.sel({wl_dim: keep}, method='nearest')
    if resolution and coarsen_from and resolution > coarsen_from:
        f = int(round(resolution / coarsen_from))
        out = out.coarsen(x=f, y=f, boundary='trim').mean()
    names = {S2_NOMINAL_WL[b]: S2_TO_TRISHNA[b] for b in S2_TO_TRISHNA}
    out = out.assign_coords(band=(wl_dim, [names.get(int(w), str(int(w))) for w in out[wl_dim].values]))
    out.attrs['emulated_sensor'] = 'TRISHNA (S2 band restriction)'
    return out


class SBAF:
    r"""Spectral band adjustment factors (SBAF) from Sentinel-2 to TRISHNA bands.

    For each pair of bands, :math:`X_{TRISHNA} = a\, X_{S2}` with :math:`a` the median of the
    ratios of the band-averaged spectra of a hyperspectral library (multiplicative factor
    fitted on relative differences, so that the bright spectra do not dominate). A band is
    fitted only if the library covers it, and if at least 5 spectra are above ``floor``;
    otherwise its factor is NaN.

    Use Rrs libraries for L2A products and radcalnet_oc ``Rtoa`` simulations for L1C.

    :param spectra: hyperspectral library, :class:`xarray.DataArray` (sample, ``wl_dim``)
    :param wl_dim: name of the wavelength dimension (nm)
    :param s2: Sentinel-2 :class:`~trishna_tbx.sensor.Sensor`; default :func:`~trishna_tbx.sensor.make_s2`
    :param trishna: TRISHNA :class:`~trishna_tbx.sensor.Sensor`
    :param floor: minimum band value used in the fit

    :ivar coef: dict, per TRISHNA band: S2 band, factor ``a``, offset ``b`` (0), median absolute
        relative difference before (``mdape_raw``) and after (``mdape_adjusted``) the adjustment,
        2.5-97.5 % interval of the ratios and number of spectra
    """

    def __init__(self, spectra, wl_dim='wl', s2=None, trishna=TRISHNA, floor=1e-4):
        s2 = s2 or make_s2()
        e0 = solar_irradiance()
        wl_max = float(spectra[wl_dim].max())
        pairs = {b.s2_proxy: b.name for b in trishna.bands
                 if b.s2_proxy and b.center + b.fwhm / 2 <= wl_max}  # library must cover the band
        x = s2.subset(list(pairs)).convolve(spectra, wl_dim, weight=e0)
        y = trishna.subset(list(pairs.values())).convolve(spectra, wl_dim, weight=e0)
        self.coef = {}
        for bs2, btri in pairs.items():
            xs, ys = x.sel(band=bs2).values.ravel(), y.sel(band=btri).values.ravel()
            ok = np.isfinite(xs) & np.isfinite(ys) & (ys > floor) & (xs > floor)
            if ok.sum() < 5:  # too dark (e.g. NIR of clear waters): no factor
                self.coef[btri] = dict(s2=bs2, a=np.nan, b=0., mdape_raw=np.nan, mdape_adjusted=np.nan,
                                       p95_ratio=None, n=int(ok.sum()))
                continue
            # multiplicative factor, fitted on relative differences (bright spectra do not dominate)
            a = float(np.median(ys[ok] / xs[ok]))
            b = 0.
            rel = np.median(np.abs(a * xs[ok] - ys[ok]) / ys[ok])
            raw = np.median(np.abs(xs[ok] - ys[ok]) / ys[ok])
            self.coef[btri] = dict(s2=bs2, a=a, b=b, mdape_raw=raw, mdape_adjusted=rel,
                                   p95_ratio=np.percentile(ys[ok] / xs[ok], [2.5, 97.5]).round(3), n=int(ok.sum()))

    def table(self):
        """Factors and statistics as a :class:`pandas.DataFrame` (one row per TRISHNA band)."""
        import pandas as pd
        return pd.DataFrame(self.coef).T

    def apply(self, ds_restricted, var='Rrs', band_dim='wl'):
        """Apply the factors to a band-restricted dataset (output of :func:`restrict_bands`).

        :param ds_restricted: :class:`xarray.Dataset` with the coordinate ``band`` (TRISHNA names)
        :param var: variable to adjust
        :param band_dim: dimension of the bands
        :return: adjusted copy of ``ds_restricted``
        """
        out = ds_restricted.copy()
        for i, name in enumerate(out['band'].values):
            if name in self.coef:
                c = self.coef[name]
                out[var][{band_dim: i}] = c['a'] * out[var][{band_dim: i}] + c['b']
        return out


def add_noise(rho_toa, sza, sensor=TRISHNA, band_dim='band', seed=None):
    """Add Gaussian TOA reflectance noise from the SNR specification (white, per pixel).

    :param rho_toa: TOA reflectance, :class:`xarray.DataArray` with the dimension ``band_dim``
    :param sza: solar zenith angle (deg); its mean is used
    :param sensor: :class:`~trishna_tbx.sensor.Sensor` giving the noise
        (:meth:`~trishna_tbx.sensor.Sensor.noise_reflectance`), same bands as ``rho_toa``
    :param band_dim: dimension of the bands
    :param seed: seed of the random generator
    :return: noisy reflectance
    """
    rng = np.random.default_rng(seed)
    sigma = sensor.noise_reflectance(sza=float(np.nanmean(sza)))
    sigma = sigma.rename(band=band_dim).assign_coords({band_dim: rho_toa[band_dim]})
    eps = xr.DataArray(rng.standard_normal(rho_toa.shape), dims=rho_toa.dims, coords=rho_toa.coords)
    return rho_toa + sigma * eps


def owt_band_adjustment(rrs, s2_sensor, n=200, seed=0):
    """Sentinel-2 to TRISHNA Rrs adjustment with factors depending on the optical water type.

    The blue factor ranges from about 1.28 (clear waters) to 0.84 (very turbid): a single global
    factor is biased. The factors (:class:`SBAF`) are fitted for each Bi et al. (2024) type on
    synthetic members (:func:`trishna_tbx.owt.sample_members`); each pixel gets the factors of its
    type (spectral angle on the ``s2_sensor`` bands). Where a band is too dark for a factor (NIR
    of clear waters), the factor fitted on all types is used.

    :param rrs: :class:`xarray.DataArray` (``wl``, ...) on the S2 nominal wavelengths of the
        ``s2_sensor`` bands
    :param s2_sensor: Sentinel-2 :class:`~trishna_tbx.sensor.Sensor` restricted to the bands of
        ``rrs`` (e.g. B2, B3, B4, B8A)
    :param n: number of synthetic members per type
    :param seed: seed of the random generator
    :return: (adjusted ``rrs``, optical water type of each pixel, :class:`pandas.DataFrame` of the
        factors per type)
    """
    import pandas as pd
    from . import owt
    ref, sd = owt.load_reference('Bi2024')
    members = owt.sample_members(ref, sd, n=n, seed=seed)
    tri_of = {b.s2_proxy: b.name for b in TRISHNA.bands}

    def factors(m):
        spec = m.stack(s=[d for d in m.dims if d != 'wl']).drop_vars(['s', 'true_owt', 'member'], errors='ignore')
        t = SBAF(spec.transpose('s', 'wl')).table()
        return {S2_NOMINAL_WL[b]: t.loc[tri_of[b], 'a'] for b in s2_sensor.band_names}

    f_all = factors(members)
    table = {'all types': f_all}
    for t in ref.owt.values:
        f = factors(members.sel(true_owt=[t]))
        table[t] = {w: (f_all[w] if np.isnan(a) else a) for w, a in f.items()}
    r = owt.convolve_reference(ref, s2_sensor)
    x = rrs.sel(wl=r.wl.values, method='nearest').rename(wl='band').assign_coords(band=r.band.values)
    valid = x.notnull().all('band')
    cls, _ = owt.classify(x.where(valid), r.drop_vars('wl'), dim='band')
    cls = cls.where(valid)
    fac = xr.concat([xr.apply_ufunc(np.vectorize(lambda c, w=w: table[c][w] if isinstance(c, str) else np.nan,
                                                 otypes=[float]), cls) for w in rrs.wl.values], dim=rrs.wl)
    return rrs * fac, cls, pd.DataFrame(table).T


def l1c_dataset_as_trishna(ds, resolution=60, native=None, strict=False, trishna=True):
    """L1C dataset of GRSdriver/grs to a TRISHNA-like (or full Sentinel-2) dataset at 60 m.

    Used for the L1C subsets saved as NetCDF (e.g. around HYPERNETS or AERONET-OC sites), which
    grs processes on their native grid. The float variables (bands, angles) are block-averaged,
    the masks (``surfwater``, ``cloud_*``) take the majority of the block.

    :param ds: L1C :class:`xarray.Dataset` (``bands`` on ``wl``, ``y``, ``x``)
    :param resolution: output pixel size (m)
    :param native: pixel size of ``ds`` (m); default: from the ``x`` coordinate
    :param strict: see :func:`grs_band_set`
    :param trishna: if True, keep the TRISHNA-like bands and set ``wl_to_process`` to 490, 560,
        665, 865 and 1610 nm; if False, keep all the bands (same 60 m grid, for a fair S2 vs
        TRISHNA comparison)
    :return: :class:`xarray.Dataset`
    """
    if native is None:
        native = float(abs(ds.x[1] - ds.x[0]))
    f = int(round(resolution / native))
    out = ds
    if trishna:
        band_idx, tbp_idx = grs_band_set(strict)
        keep = [S2_NOMINAL_WL[S2_BAND_ORDER[i]] for i in band_idx]
        out = out.sel(wl=[w for w in out.wl.values if int(w) in keep])
        out.attrs = dict(ds.attrs, wl_to_process=np.array([S2_NOMINAL_WL[S2_BAND_ORDER[i]] for i in tbp_idx]),
                         emulated_sensor='TRISHNA (S2 band restriction)')
    if f > 1:
        cont = [v for v in out.data_vars if {'y', 'x'} <= set(out[v].dims) and out[v].dtype.kind == 'f']
        cat = [v for v in out.data_vars if {'y', 'x'} <= set(out[v].dims) and v not in cont]
        coarse = out[cont].coarsen(x=f, y=f, boundary='trim').mean()
        for v in cat:
            coarse[v] = (out[v].astype(float).coarsen(x=f, y=f, boundary='trim').mean() >= 0.5).astype(out[v].dtype)
        for v in out.data_vars:
            if v not in coarse:
                coarse[v] = out[v]
        coarse.attrs = out.attrs
        out = coarse
    return out
