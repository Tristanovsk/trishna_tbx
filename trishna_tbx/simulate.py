"""TOA simulations for TRISHNA and S2 with radcalnet_oc.

Water signal sources:
- in situ: AERONET-OC (``radcalnet_oc.aeronet_oc.Aeronet.open_site``), HYPERNETS (reader TBD);
- synthetic: invRrs forward model on a grid of (chl, spm, a_cdom), or OWT mean spectra.

The hyperspectral TOA output (350-2500 nm) is convolved with both sensors' SRF, giving
paired TRISHNA / S2 band values for: SBAF fitting (``emulate.SBAF``), atmospheric
correction sensitivity (AOT, aerosol model, H2O, O3), water vapour algorithm testing and
re-tuning of band-ratio algorithms (``algorithms.retune_ocx``).
"""
import numpy as np
import pandas as pd
import xarray as xr

from .sensor import TRISHNA, make_s2


def build_input(rrs, sza, aot550=0.1, tcwv=15., tco3=7e-3, tcno2=1e-5, pressure=1013.25, wind=5.,
                doy=180, ang_exp=1.0):
    """Input dataset of radcalnet_oc from hyperspectral Rrs and atmospheric parameters.

    The scalars are broadcast on ``time``; the units are those of radcalnet_oc. ``aot865`` is
    computed from ``aot550`` and ``ang_exp`` (radcalnet_oc copies it to the output; the spectral
    dependence of the aerosol comes from its aerosol models).

    :param rrs: remote-sensing reflectance (sr-1), :class:`xarray.DataArray` (``time`` or any
        sample dimension, ``wl`` in nm)
    :param sza: solar zenith angle (deg)
    :param aot550: aerosol optical thickness at 550 nm
    :param tcwv: total column water vapour (kg m-2)
    :param tco3: total column ozone (kg m-2)
    :param tcno2: total column nitrogen dioxide (kg m-2)
    :param pressure: surface pressure (hPa)
    :param wind: wind speed at 10 m (m s-1)
    :param doy: day of year
    :param ang_exp: Angström exponent
    :return: :class:`xarray.Dataset` (``time``, ``wl``)
    """
    if 'time' not in rrs.dims:
        rrs = rrs.rename({rrs.dims[0]: 'time'})
    n = rrs.sizes['time']
    f = lambda v: ('time', np.broadcast_to(np.asarray(v, float), (n,)).copy())
    return xr.Dataset({'Rrs': rrs, 'sza': f(sza), 'aot550': f(aot550), 'tcwv': f(tcwv),
                       'tco3': f(tco3), 'tcno2': f(tcno2), 'pressure': f(pressure), 'wind': f(wind),
                       'day_of_year': f(doy), 'ang_exp': f(ang_exp),
                       # aot865/ang_exp are only copied to the output by radcalnet_oc (aerosol
                       # spectral dependence comes from the aerosol_combination)
                       'aot865': f(np.asarray(aot550) * (865. / 550.) ** -np.asarray(ang_exp))})


DU_O3, DU_NO2 = 2.14485e-5, 2.055482e-5  # kg m-2 per Dobson unit


def aeronet_oc_input(site, period, level='lev20', hyp_dir=None, rrs_var='Rrs_Lwn_fQ', max_records=None):
    """Input dataset of radcalnet_oc from the measurements of an AERONET-OC site.

    Water signal: the hyperspectral invRrs reconstruction of the AERONET-OC spectra when
    ``hyp_dir`` is given, otherwise the multispectral ``rrs_var`` interpolated between the
    SeaPRISM bands (crude around 700 nm). Atmosphere: AOT at 550 nm (log-linear interpolation of
    the sun-photometer AOT), Angström exponent (440-870 nm), water vapour, ozone, NO2, pressure and
    wind of each record.

    :param site: AERONET-OC site, e.g. ``'Venise'``
    :param period: str or slice for ``.sel(time=...)``, e.g. ``'2021-06'``
    :param level: AERONET-OC level, ``'lev20'`` or ``'lev15'``
    :param hyp_dir: folder of the invRrs reconstructions
        (``<hyp_dir>/<site>/<level>/<site>_<level>_YYYY_MM_v2.nc``)
    :param rrs_var: multispectral Rrs variable used without ``hyp_dir``
    :param max_records: maximum number of records
    :return: :class:`xarray.Dataset` of :func:`build_input`, plus ``chl_aeronet``
    """
    import glob
    from radcalnet_oc.aeronet_oc import Aeronet
    ds = Aeronet.open_site(site, level).sel(time=period).load()
    ok = ds[rrs_var].notnull().sum('wl') > 5
    ds = ds.sel(time=ok)
    if max_records:
        ds = ds.isel(time=slice(0, max_records))
    if hyp_dir:
        files = sorted(glob.glob(f'{hyp_dir}/{site}/{level}/{site}_{level}_*_v2.nc'))
        hyp = xr.open_mfdataset(files).Rrs.reindex(time=ds.time, method='nearest',
                                                   tolerance=np.timedelta64(1, 'm')).load()
        found = hyp.notnull().any('wl')
        rrs = hyp.sel(time=found).fillna(0.)  # NaN only beyond ~1050 nm (black water)
        ds = ds.sel(time=rrs.time)
    else:
        rrs = ds[rrs_var].dropna('wl', how='all').interpolate_na('wl')
    aot = ds['aot'].dropna('wl', how='all')
    la = np.log(aot.where(aot > 0))
    aot550 = np.exp(la.interp(wl=550.))  # log(aot) linear in wl: adequate for 500-675 nm
    aot550 = aot550.fillna(np.exp(la.interpolate_na('wl').interp(wl=550.)))
    alpha = -(np.log(aot.sel(wl=870, method='nearest')) - np.log(aot.sel(wl=440, method='nearest'))) / np.log(870 / 440)
    doy = ds.time.dt.dayofyear
    inp = build_input(rrs.drop_vars([c for c in rrs.coords if c not in ('time', 'wl')]),
                      sza=ds.sza.values, aot550=aot550.values, tcwv=10 * ds.tcwv.values,
                      tco3=DU_O3 * ds.tco3.values, tcno2=DU_NO2 * ds.tno2.values,
                      pressure=ds.pressure.values, wind=ds.wind.fillna(5.).values, doy=doy.values,
                      ang_exp=alpha.values)
    return inp.assign_coords(time=ds.time.values).assign(chl_aeronet=('time', ds.chl.values))


def simulate_toa(input_db, vza=(5., 20., 35.), azi=(90.,), aerosol_combination=None, **kw):
    """Run radcalnet_oc.

    :param input_db: input dataset, see :func:`build_input`
    :param vza: viewing zenith angles (deg)
    :param azi: relative azimuth angles (deg)
    :param aerosol_combination: proportions of the aerosol models of radcalnet_oc; default: its
        configuration
    :param kw: other arguments of :class:`radcalnet_oc.process.Process`
    :return: hyperspectral TOA :class:`xarray.Dataset` (``Rtoa``, ``Ratm``, ``Ed``, ``E0``,
        ``Rrs``...), 350-2500 nm
    """
    from radcalnet_oc.process import Process
    vza = xr.DataArray(list(vza), dims='vza', coords={'vza': list(vza)})
    azi = xr.DataArray(list(azi), dims='azi', coords={'azi': list(azi)})
    p = Process(input_db=input_db, vza=vza, azi=azi, **kw)
    p.execute() if aerosol_combination is None else p.execute(aerosol_combination=aerosol_combination)
    return p.radcalnet_db


def to_sensors(toa, sensors=None, variables=('Rtoa', 'Ratm', 'Rrs')):
    """Band-average the hyperspectral simulation for each sensor (radcalnet_oc ``Spectral``).

    :param toa: output of :func:`simulate_toa`
    :param sensors: list of :class:`~trishna_tbx.sensor.Sensor`; default TRISHNA and Sentinel-2A
    :param variables: variables to convolve
    :return: dict sensor name -> :class:`xarray.Dataset` on ``band``
    """
    from radcalnet_oc.lut import Spectral
    sensors = sensors or [TRISHNA, make_s2()]
    out = {}
    for s in sensors:
        spec = Spectral.from_srf(s.srf, wl_dim='wl', band_dim='band', central_wl=s.centers.values)
        conv = {}
        for v in variables:
            if v in toa:
                c = spec.convolve(toa[v], name=v)
                c = c.squeeze('variable', drop=True) if 'variable' in c.dims else c
                if 'band' in c.coords and 'band' not in c.dims:  # band names as dimension
                    c = c.swap_dims(wl='band')
                conv[v] = c
        out[s.name] = xr.Dataset(conv)
    return out


def sensitivity_grid(rrs, sza=30., aot=(0.05, 0.1, 0.2, 0.4), tcwv=(5., 15., 30., 50.), **kw):
    """Simulations of one Rrs spectrum for all the combinations of AOT and water vapour.

    The combinations are run at once (one radcalnet_oc run).

    :param rrs: Rrs (sr-1), :class:`xarray.DataArray` on ``wl``
    :param sza: solar zenith angle (deg)
    :param aot: aerosol optical thicknesses at 550 nm
    :param tcwv: total columns of water vapour (kg m-2)
    :param kw: arguments of :func:`simulate_toa`
    :return: output of radcalnet_oc with ``time`` replaced by (``aot550_sim``, ``tcwv_sim``)
    """
    import itertools
    combos = list(itertools.product(aot, tcwv))
    rr = xr.concat([rrs] * len(combos), dim='time').assign_coords(time=np.arange(len(combos)))
    inp = build_input(rr, sza, aot550=[c[0] for c in combos], tcwv=[c[1] for c in combos])
    toa = simulate_toa(inp, **kw)
    idx = pd.MultiIndex.from_tuples(combos, names=['aot550_sim', 'tcwv_sim'])
    return toa.assign_coords(xr.Coordinates.from_pandas_multiindex(idx, 'time')).unstack('time')
