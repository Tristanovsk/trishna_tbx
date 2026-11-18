"""Readers for GRS L2A products (beam/snap-compliant ``Rrs_<wl>`` or stacked ``Rrs``) and CAMS."""
import glob
import os
import re

import numpy as np
import xarray as xr


def parse_window(window):
    """Pixel window given as a string to a dict of slices for :meth:`xarray.Dataset.isel`.

    :param window: ``'y=0:3600,x=0:960'``, or ``'0:3600,0:960'`` (y then x), or None
    :return: dict, e.g. ``{'y': slice(0, 3600), 'x': slice(0, 960)}`` (empty for None)
    """
    if window is None:
        return {}
    out = {}
    for i, part in enumerate(window.split(',')):
        if '=' in part:
            dim, rng = part.split('=')
        else:
            dim, rng = ('y', 'x')[i], part
        a, b = rng.split(':')
        out[dim.strip()] = slice(int(a), int(b))
    return out


def open_grs_l2a(path, window=None, chunks={}, level=0):
    """Open a GRS L2A product (NetCDF, Zarr or Zarr pyramid) with a stacked ``Rrs`` variable.

    The products of the beam/SNAP profile (one variable ``Rrs_<wl>`` per band) are stacked.
    For a Zarr pyramid of GRS v3 (levels ``0``, ``1``... in the store, ``0`` at full resolution),
    the level ``level`` is read.

    :param path: NetCDF file or Zarr store
    :param window: pixel window, see :func:`parse_window`
    :param chunks: dask chunks of :func:`xarray.open_dataset`
    :param level: level of a Zarr pyramid
    :return: :class:`xarray.Dataset` with ``Rrs`` (``wl``, ``y``, ``x``)
    """
    path = str(path).rstrip('/')
    if path.endswith('.zarr'):
        ds = xr.open_zarr(path)
        if not ds.data_vars and os.path.isdir(os.path.join(path, str(level))):  # pyramid
            ds = xr.open_zarr(os.path.join(path, str(level)))
    else:
        ds = xr.open_dataset(path, chunks=chunks)
    ds = ds.isel(parse_window(window))
    if 'Rrs' not in ds:
        wls = sorted(int(v.split('_')[1]) for v in ds.data_vars if re.fullmatch(r'Rrs_\d+', v))
        ds['Rrs'] = xr.concat([ds[f'Rrs_{w}'] for w in wls], dim=xr.DataArray(wls, dims='wl', name='wl'))
        ds = ds.drop_vars([f'Rrs_{w}' for w in wls])
    ds['Rrs'] = ds.Rrs.assign_coords(wl=ds.wl.astype(float)).transpose('wl', 'y', 'x')
    return ds


def water_mask(ds):
    """Water pixels of a GRS L2A product.

    ``mask == 0`` for the recent versions of GRS, otherwise ``ndwi > 0`` (beam profile), otherwise
    Rrs(1610) < 0.001; always with a valid Rrs at 560 nm.

    :param ds: L2A dataset of :func:`open_grs_l2a`
    :return: boolean :class:`xarray.DataArray` (``y``, ``x``)
    """
    valid = ds.Rrs.sel(wl=560, method='nearest').notnull()
    if 'mask' in ds:
        return (ds['mask'] == 0) & valid
    if 'ndwi' in ds:
        return (ds['ndwi'] > 0) & valid
    return (ds.Rrs.sel(wl=1610, method='nearest') < 1e-3) & valid


def to_60m(rrs, native=20.):
    """Block-average an image to 60 m (TRISHNA pixel size at nadir).

    :param rrs: :class:`xarray.DataArray` or :class:`xarray.Dataset` with ``x`` and ``y``
    :param native: native pixel size (m)
    :return: averaged image
    """
    f = int(round(60 / native))
    return rrs.coarsen(x=f, y=f, boundary='trim').mean() if f > 1 else rrs


def list_l2a(l2a_dir, tile):
    """GRS L2A files of a tile, sorted by acquisition time.

    :param l2a_dir: folder of the products (``<l2a_dir>/*<tile>*/*<tile>*.nc``)
    :param tile: Sentinel-2 tile, e.g. ``'T33PVQ'``
    :return: list of (:class:`numpy.datetime64`, path)
    """
    files = sorted(glob.glob(f'{l2a_dir}/*{tile}*/*{tile}*[0-9].nc'))
    stamps = [re.search(r'_(\d{8}T\d{6})_', f.split('/')[-1]).group(1) for f in files]
    times = [np.datetime64(f'{s[:4]}-{s[4:6]}-{s[6:8]}T{s[9:11]}:{s[11:13]}:{s[13:15]}') for s in stamps]
    return sorted(zip(times, files))


def cams_value(cams_file, var, time, lon, lat):
    """CAMS field at the closest valid time, interpolated linearly in longitude and latitude.

    :param cams_file: CAMS forecast file of grs (``valid_time`` on two dimensions)
    :param var: variable, e.g. ``'aod550'`` or ``'tcwv'``
    :param time: date and time
    :param lon: longitude(s) (deg, -180-180 or 0-360)
    :param lat: latitude(s) (deg)
    :return: :class:`numpy.ndarray`
    """
    c = xr.open_dataset(cams_file)[var]
    vt = c.valid_time.values
    i, j = np.unravel_index(np.argmin(np.abs(vt - np.datetime64(time))), vt.shape)
    f = c.isel({c.valid_time.dims[0]: i, c.valid_time.dims[1]: j})
    lon = np.where(np.asarray(lon) > 180, np.asarray(lon) - 360, lon)
    return f.interp(longitude=xr.DataArray(np.atleast_1d(lon)), latitude=xr.DataArray(np.atleast_1d(lat))).values
