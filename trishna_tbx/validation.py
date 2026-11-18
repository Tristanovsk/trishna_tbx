"""Matchup validation with HYPERNETS (WATERHYPERNET) for Sentinel-2 and TRISHNA-like bands.

HYPERNETS gives hyperspectral Rrs (400-900 nm) every ~30 min: it can be convolved with any
band set, so the *same* in situ measurement validates S2 and TRISHNA-like products. The
TRISHNA WV band (910 nm) is outside the 400-900 nm range of the distributed L2A data.

Workflow:
1. ``load_hypernets`` (QC) -> ``convolve`` to S2 / TRISHNA bands;
2. ``find_matchups`` between HYPERNETS records and S2 L1C subsets (NetCDF, GRSdriver format);
3. ``process_matchup`` runs GRS at 60 m on the L1C subset with TRISHNA-like and with all S2
   bands (cached), and extracts the window statistics around the site;
4. ``metrics`` per band.
"""
import glob
import re
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

from .sensor import solar_irradiance


def load_hypernets(path, quality_flag=0, max_sza=70., wl_max=900.):
    """HYPERNETS water L2A time series (one NetCDF file per site) with a basic quality control.

    :param path: NetCDF file of the stacked HYPERNETS L2A data (``Rrs`` on ``time``, ``wl``)
    :param quality_flag: value of ``quality_flag`` kept (None: all records)
    :param max_sza: maximum solar zenith angle (deg)
    :param wl_max: maximum wavelength kept (nm)
    :return: :class:`xarray.Dataset` of the records with a valid Rrs at 560 nm, with the site
        attributes
    """
    h = xr.open_dataset(path)
    ok = h.Rrs.sel(wl=560, method='nearest').notnull() & (h.solar_zenith_angle < max_sza)
    if quality_flag is not None and 'quality_flag' in h:
        ok = ok & (h.quality_flag == quality_flag)
    h = h.sel(time=ok.values).sel(wl=slice(None, wl_max))
    return h


def convolve(rrs, sensor, bands=None):
    """Band-averaged Rrs (weighted by the solar irradiance) of hyperspectral data.

    :param rrs: :class:`xarray.DataArray` (``time``, ``wl``)
    :param sensor: :class:`~trishna_tbx.sensor.Sensor`
    :param bands: names of the bands; default all
    :return: :class:`xarray.DataArray` (``time``, ``band``)
    """
    s = sensor.subset(bands) if bands else sensor
    return s.convolve(rrs, 'wl', weight=solar_irradiance())


def find_matchups(hyp_time, l1c_dir, max_dt='1h'):
    """Sentinel-2 L1C subsets with a HYPERNETS record within ``max_dt``.

    :param hyp_time: times of the HYPERNETS records
    :param l1c_dir: folder of the L1C subsets (NetCDF, searched recursively, names
        ``*_MSIL1C_<YYYYMMDDTHHMMSS>_N*.nc``)
    :param max_dt: maximum time difference (pandas time string)
    :return: :class:`pandas.DataFrame`: ``s2_time``, ``hyp_time``, ``dt_min``, ``l1c``
    """
    files = sorted(glob.glob(f'{l1c_dir}/**/*MSIL1C*.nc', recursive=True))
    ht = pd.DatetimeIndex(hyp_time)
    rows = []
    for f in files:
        t = pd.Timestamp(re.search(r'_(\d{8}T\d{6})_N', Path(f).name).group(1))
        i = np.argmin(np.abs(ht - t))
        dt = abs(ht[i] - t)
        if dt <= pd.Timedelta(max_dt):
            rows.append({'s2_time': t, 'hyp_time': ht[i], 'dt_min': dt.total_seconds() / 60, 'l1c': f})
    return pd.DataFrame(rows)


def _site_window(l2, x, y, half=1):
    ix = int(np.argmin(np.abs(l2.x.values - x)))
    iy = int(np.argmin(np.abs(l2.y.values - y)))
    return l2.isel(x=slice(ix - half, ix + half + 1), y=slice(iy - half, iy + half + 1))


def process_matchup(l1c_file, cams_dir, cache_dir, site_xy, half=1, resolution=60, crop_m=15000, grs_kw=None):
    """Run GRS on one L1C subset with the TRISHNA-like and with all S2 bands; site statistics.

    The subset is cropped to ``crop_m`` around the site and averaged to ``resolution`` once;
    the TRISHNA-like band set is taken from this averaged product
    (:func:`~trishna_tbx.emulate.l1c_dataset_as_trishna`), so that both band sets have the same
    grid. A failure of GRS (e.g. no water, fully cloudy subset) is logged and the band set is
    skipped. The result is cached in ``cache_dir`` (``<l1c name>_win<n>_crop<km>km.nc``).

    :param l1c_file: L1C subset (NetCDF of GRSdriver)
    :param cams_dir: folder of the monthly CAMS files ``cams_forecast_YYYY-MM.nc``
    :param cache_dir: folder of the results
    :param site_xy: (x, y) of the site in the projection of the subset (m)
    :param half: half size of the window (pixels): ``2 * half + 1`` pixels wide
    :param resolution: pixel size of the processing (m)
    :param crop_m: size of the crop around the site (m)
    :param grs_kw: other arguments of ``grs.Process().execute``
    :return: :class:`xarray.Dataset` (``bandset``, ``wl``): median and standard deviation of Rrs
        on the water pixels of the window, number of water pixels ``n_valid``, median ``aot550``;
        None if GRS failed for both band sets
    """
    import logging
    import grs
    from .emulate import l1c_dataset_as_trishna
    logging.getLogger().setLevel(logging.WARNING)
    cache = Path(cache_dir) / (Path(l1c_file).stem + f'_win{2 * half + 1}_crop{crop_m / 1000:g}km.nc')
    if cache.exists():
        return xr.open_dataset(cache).load()
    ds = xr.open_dataset(l1c_file)
    x0, y0 = site_xy
    sl_y = slice(y0 + crop_m / 2, y0 - crop_m / 2) if ds.y[0] > ds.y[-1] else slice(y0 - crop_m / 2, y0 + crop_m / 2)
    ds = ds.sel(x=slice(x0 - crop_m / 2, x0 + crop_m / 2), y=sl_y)
    t = pd.Timestamp(ds.time.values)
    cams = f'{cams_dir}/cams_forecast_{t:%Y-%m}.nc'
    full = l1c_dataset_as_trishna(ds, resolution=resolution, trishna=False).load()
    out = []
    for name, tri in (('TRISHNA', True), ('S2', False)):
        prod = l1c_dataset_as_trishna(full, resolution=resolution, native=resolution, trishna=True) if tri else full
        p = grs.Process()
        try:
            res = p.execute(prod, odir=str(Path(cache_dir) / 'tmp'), cams_file=cams, **(grs_kw or {}))
        except Exception as exc:  # e.g. fully cloudy subset
            logging.warning(f'{Path(l1c_file).name} ({name}): GRS failed: {exc}')
            continue
        if res is None or getattr(p, 'no_water', False):
            continue
        l2 = p.l2a.l2_prod[['Rrs', 'aot550', 'mask']]
        win = _site_window(l2, *site_xy, half=half).load()
        water = win.mask == 0
        r = win.Rrs.where(water)
        out.append(xr.Dataset({
            'Rrs': r.median(['x', 'y']), 'Rrs_std': r.std(['x', 'y']),
            'n_valid': water.sum(), 'aot550': win.aot550.where(water).median()},
        ).assign_coords(wl=r.wl.astype(int).values).expand_dims(bandset=[name]))
    if not out:
        return None
    res = xr.concat(out, 'bandset', join='outer').assign_coords(s2_time=t)
    for v in res.variables.values():  # GRS attributes (booleans...) are not NetCDF-safe
        v.attrs = {}
    res.attrs = {'l1c': str(l1c_file), 'window': f'{2 * half + 1}x{2 * half + 1} at {resolution} m'}
    cache.parent.mkdir(parents=True, exist_ok=True)
    res.to_netcdf(cache)
    return res


def metrics(x, y):
    """Validation statistics of satellite values against in situ values.

    :param x: in situ values (1-D)
    :param y: satellite values (1-D); the pairs with a NaN are dropped
    :return: :class:`pandas.Series`: ``N``, median bias, median absolute percentage difference
        ``MdAPD (%)``, ``RMSD``, slope and intercept of the linear fit, correlation ``r`` (only
        ``N`` if fewer than 3 pairs)
    """
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    if ok.sum() < 3:
        return pd.Series({'N': ok.sum()})
    slope, inter = np.polyfit(x, y, 1)
    return pd.Series({'N': ok.sum(), 'bias': np.median(y - x), 'MdAPD (%)': 100 * np.median(np.abs(y - x) / np.abs(x)),
                      'RMSD': np.sqrt(np.mean((y - x) ** 2)), 'slope': slope, 'intercept': inter,
                      'r': np.corrcoef(x, y)[0, 1]})


def find_l2a(l2a_dir, s2_time, layout='%Y/%m/%d'):
    """Archived GRS L2A product of the same Sentinel-2 acquisition.

    Only the date folder is searched (fast on large archives).

    :param l2a_dir: root of the archive
    :param s2_time: acquisition time
    :param layout: date folders (strftime format); default ``<l2a_dir>/YYYY/MM/DD/<name>/<name>.nc``
        with ``<name>`` containing ``MSIL2A_<YYYYMMDDTHHMMSS>``
    :return: path, or None if not found
    """
    t = pd.Timestamp(s2_time)
    stamp = t.strftime('%Y%m%dT%H%M%S')
    files = [f for f in glob.glob(f'{l2a_dir}/{t.strftime(layout)}/*MSIL2A_{stamp}*/*MSIL2A_{stamp}*.nc')
             if not f.endswith('_anc.nc')]
    return sorted(files)[0] if files else None


def extract_l2a_window(path, site_xy, half_m=90, name='S2 archive'):
    """Site statistics, as in :func:`process_matchup`, from an existing GRS L2A product.

    :param path: GRS L2A NetCDF file
    :param site_xy: (x, y) of the site in the projection of the product (m)
    :param half_m: half size of the window (m), on the native grid of the product
    :param name: name given to the ``bandset`` coordinate
    :return: :class:`xarray.Dataset` (``bandset``, ``wl``)
    """
    ds = xr.open_dataset(path)
    x0, y0 = site_xy
    sl_y = slice(y0 + half_m, y0 - half_m) if ds.y[0] > ds.y[-1] else slice(y0 - half_m, y0 + half_m)
    win = ds.sel(x=slice(x0 - half_m, x0 + half_m), y=sl_y).load()
    water = (win['mask'] == 0) if 'mask' in win else win.Rrs.sel(wl=560).notnull()
    r = win.Rrs.where(water)
    out = xr.Dataset({'Rrs': r.median(['x', 'y']), 'Rrs_std': r.std(['x', 'y']),
                      'n_valid': water.sum(), 'aot550': win.aot550.where(water).median()})
    out = out.drop_vars([c for c in out.coords if c != 'wl']).assign_coords(wl=out.wl.astype(int).values)
    for v in out.variables.values():
        v.attrs = {}
    return out.expand_dims(bandset=[name])
