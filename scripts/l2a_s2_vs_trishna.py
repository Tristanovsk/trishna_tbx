"""Compare water products from a GRS S2 L2A image: full S2 bands vs TRISHNA-like bands.

usage: python scripts/l2a_s2_vs_trishna.py <GRS_L2A.nc> [--native 20] [--out out.nc] [--crop y0:y1,x0:x1]

Steps: stack Rrs_<wl> -> Rrs(wl), water mask, block-average to 60 m, then
 - S2 full: OWT (Bi2024, S2 SRF) + all S2-covered algorithms
 - S2-as-TRISHNA: keep 490/560/665/865 (+1610), OWT on those bands, TRISHNA-covered algorithms
Agreement statistics are printed; the 60 m product is written to --out.
"""
import argparse

import numpy as np
import xarray as xr

from trishna_tbx import algorithms, owt
from trishna_tbx.sensor import make_s2


def open_grs_l2a(path, crop=None):
    ds = xr.open_dataset(path, chunks={})
    if crop:  # 'y0:y1,x0:x1' in pixels
        (y0, y1), (x0, x1) = [map(int, c.split(':')) for c in crop.split(',')]
        ds = ds.isel(y=slice(y0, y1), x=slice(x0, x1))
    if 'Rrs' not in ds:  # beam / snap-compliant profile
        wls = [int(v.split('_')[1]) for v in ds.data_vars if v.startswith('Rrs_') and v.split('_')[1].isdigit()]
        ds['Rrs'] = xr.concat([ds[f'Rrs_{w}'] for w in wls], dim=xr.DataArray(wls, dims='wl', name='wl'))
    return ds


def water(ds):
    # GRS L2A: 'mask' (0 = water) in recent versions; older 'beam' products are NaN outside water
    if 'mask' in ds:
        return ds['mask'] == 0
    if 'ndwi' in ds:  # GRS ndwi = (R490 - R842)/(R490 + R842) on Rrs: water > 0
        return (ds['ndwi'] > 0) & ds.Rrs.sel(wl=560, method='nearest').notnull()
    return (ds.Rrs.sel(wl=1610, method='nearest') < 1e-3) & (ds.Rrs.sel(wl=560, method='nearest') > 0)


def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument('l2a')
    p.add_argument('--native', type=float, default=20.)
    p.add_argument('--crop', default=None, help="pixel window 'y0:y1,x0:x1'")
    p.add_argument('--out', default=None)
    a = p.parse_args(argv)

    ds = open_grs_l2a(a.l2a, a.crop)
    rrs = ds.Rrs.where(water(ds)).transpose('wl', 'y', 'x')
    f = int(round(60 / a.native))
    rrs = rrs.coarsen(x=f, y=f, boundary='trim').mean().load()  # 60 m
    rrs = rrs.where((rrs.sel(wl=[490, 560, 665]) > 0).all('wl'))
    print(f'60 m water pixels: {int(rrs.sel(wl=560).notnull().sum())}')

    s2 = make_s2()
    full = s2.subset(['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8A'], 'S2-full')
    tri_like = s2.subset(['B2', 'B3', 'B4', 'B8A'], 'S2-as-TRISHNA')
    mean, _ = owt.load_reference('Bi2024')
    out = {}
    for sens in (full, tri_like):
        ref = owt.convolve_reference(mean, sens)
        x = rrs.sel(wl=ref.wl.values, method='nearest').rename(wl='band').assign_coords(band=ref.band.values)
        cls, _ = owt.classify(x.where(x.notnull().all('band')), ref.drop_vars('wl'), dim='band')
        out[f'owt_{sens.name}'] = cls.where(x.notnull().all('band'))
    a_ = out['owt_S2-full'].values.ravel()
    b_ = out['owt_S2-as-TRISHNA'].values.ravel()
    ok = np.array([isinstance(v, str) for v in a_]) & np.array([isinstance(v, str) for v in b_])
    print(f'OWT agreement S2-full vs S2-as-TRISHNA: {100 * np.mean(a_[ok] == b_[ok]):.1f}% ({ok.sum()} px)')
    import pandas as pd
    print(pd.crosstab(a_[ok], b_[ok], rownames=['S2-full'], colnames=['S2-as-TRISHNA']))

    prod_full = algorithms.run_all(rrs, s2)
    prod_tri = algorithms.run_all(rrs.sel(wl=[490, 560, 665, 865]), tri_like)
    print('\nproducts S2 only   :', sorted(set(prod_full) - set(prod_tri)))
    print('products both      :', sorted(set(prod_tri)))
    q = lambda v: np.nanpercentile(v, [10, 50, 90]).round(3)
    for v in prod_tri:
        print(f'  {v:22s} p10/50/90 = {q(prod_tri[v].values)}')
    if 'Chla_index_NDCI' in prod_full:
        print(f'  {"NDCI (S2 only)":22s} p10/50/90 = {q(prod_full["Chla_index_NDCI"].values)}')
    if a.out:
        xr.merge([prod_full, xr.Dataset(out).astype(str)]).to_netcdf(a.out)


if __name__ == '__main__':
    main()
