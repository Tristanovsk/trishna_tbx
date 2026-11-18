"""Optical Water Types with a band-limited sensor.

Differences with GRSl2bgen.owt (on purpose, for the TRISHNA study):
- reference spectra are *convolved* with each sensor SRF (GRSl2bgen interpolates them at
  band centres, which is wrong for 60-70 nm wide TRISHNA bands);
- the band set is explicit (GRSl2bgen default ``wl_range=slice(350, 800)`` silently drops
  the 860 nm band, leaving TRISHNA with only 3 bands);
- a Monte-Carlo benchmark measures how often a sensor recovers the class given by the
  hyperspectral spectrum (the "truth"), with optional radiometric noise.

Reference sets come from the GRSl2bgen data folder:
  Bi2024 (10 types, 400-900 nm), Spyrakos2018 (13 inland types, 400-800 nm).
"""
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

from .sensor import solar_irradiance

#: Labels of the 10 optical water types of Bi et al. (2024).
BI2024_LABELS = ['1', '2', '3a', '3b', '4a', '4b', '5a', '5b', '6', '7']


def _data_dir():
    import GRSl2bgen
    return Path(GRSl2bgen.__file__).parent / 'data'


def load_reference(database='Bi2024', param='m_Rrs'):
    """Reference spectra of the optical water types (data of GRSl2bgen).

    :param database: ``'Bi2024'`` (Bi et al. 2024, 10 types, 400-900 nm) or ``'Spyrakos2018'``
        (Spyrakos et al. 2018, 13 inland types, 400-800 nm, standardized spectra)
    :param param: variable of the Bi et al. (2024) table, ``'m_Rrs'`` (mean Rrs) by default
    :return: (mean spectra, standard deviation), :class:`xarray.DataArray` (``owt``, ``wl``); the
        standard deviation is None for ``'Spyrakos2018'``
    """
    if database == 'Bi2024':
        df = pd.read_csv(_data_dir() / 'Bi_etal_2024_OWT_mean_spec_v01.csv')
        ds = df.set_index(['type', 'wavelen']).to_xarray().rename(type='owt', wavelen='wl')
        ds = ds.assign_coords(owt=BI2024_LABELS[:ds.sizes['owt']])
        sd = 'sd_' + param.split('_', 1)[1]
        return ds[param].rename('Rrs'), ds[sd].rename('sd')
    if database == 'Spyrakos2018':
        df = pd.read_csv(_data_dir() / 'Spyrakos_et_al_2018_OWT_inland_mean_standardised.csv', index_col=0)
        da = xr.DataArray(df.values, dims=('owt', 'wl'),
                          coords={'owt': df.index.astype(str), 'wl': df.columns.astype(float)})
        return da.rename('Rrs'), None
    raise ValueError(database)


#: Merge of the Bi et al. (2024) types for TRISHNA (7 classes). The sub-types that are nearly
#: collinear on the TRISHNA bands, with a spectral angle below 6 deg (1/2, 3a/3b, 4a/4b), are merged.
TRISHNA_MERGE = {'1': '1-2', '2': '1-2', '3a': '3', '3b': '3', '4a': '4', '4b': '4',
                 '5a': '5a', '5b': '5b', '6': '6', '7': '7'}


def merge_classes(cls, mapping=TRISHNA_MERGE):
    """Map optical water type labels to a merged class set.

    :param cls: :class:`xarray.DataArray` of labels (str, NaN allowed)
    :param mapping: dict label -> merged label, default :data:`TRISHNA_MERGE` (7 classes)
    :return: :class:`xarray.DataArray` of merged labels
    """
    return xr.apply_ufunc(np.vectorize(lambda v: mapping.get(v, v) if isinstance(v, str) else v,
                                       otypes=[object]), cls)


def sam(spectra, ref, dim='band'):
    r"""Spectral angle between spectra and reference spectra.

    :math:`\alpha = \arccos\left(\sum_k x_k r_k \,/\, \sqrt{\sum_k x_k^2 \sum_k r_k^2}\right)`

    :param spectra: :class:`xarray.DataArray` (..., ``dim``)
    :param ref: reference spectra, :class:`xarray.DataArray` (``owt``, ``dim``)
    :param dim: spectral dimension
    :return: angle (rad), :class:`xarray.DataArray` (..., ``owt``)
    """
    num = (spectra * ref).sum(dim)
    den = np.sqrt((spectra ** 2).sum(dim) * (ref ** 2).sum(dim))
    return np.arccos(np.clip(num / den, -1, 1))


def classify(spectra, ref, dim='band', distance='sam'):
    """Optical water type of each spectrum: the closest reference spectrum.

    :param spectra: :class:`xarray.DataArray` (..., ``dim``)
    :param ref: reference spectra (``owt``, ``dim``) on the same bands, e.g. from
        :func:`convolve_reference`
    :param dim: spectral dimension
    :param distance: ``'sam'`` (spectral angle, :func:`sam`, insensitive to the amplitude) or
        any other value for the Euclidean distance of the log spectra (amplitude-aware)
    :return: (type of each spectrum, distance to every type)
    """
    if distance == 'sam':
        d = sam(spectra, ref, dim)
    else:  # euclidean on log spectra (amplitude-aware)
        d = np.sqrt(((np.log(spectra.clip(1e-6)) - np.log(ref.clip(1e-6))) ** 2).sum(dim))
    return d.idxmin('owt'), d


def convolve_reference(ref, sensor, wl_max=900.):
    """Reference spectra averaged over the bands of a sensor (weighted by the solar irradiance).

    The bands that are not fully inside the spectral range of the references (and below
    ``wl_max``) are dropped.

    :param ref: reference spectra (``owt``, ``wl``)
    :param sensor: :class:`~trishna_tbx.sensor.Sensor`
    :param wl_max: maximum wavelength of the bands kept (nm)
    :return: :class:`xarray.DataArray` (``owt``, ``band``)
    """
    names = [b.name for b in sensor.bands if b.center + b.fwhm / 2 <= min(wl_max, float(ref.wl.max()))]
    return sensor.subset(names).convolve(ref, 'wl', weight=solar_irradiance())


def separability(ref, sensor):
    """Spectral angle between the reference types once averaged over the sensor bands.

    A small angle (a few degrees) means that the two types are hardly separable with the
    sensor.

    :param ref: reference spectra (``owt``, ``wl``)
    :param sensor: :class:`~trishna_tbx.sensor.Sensor`
    :return: angle (deg), :class:`xarray.DataArray` (``owt_i``, ``owt_j``)
    """
    conv = convolve_reference(ref, sensor)
    a = conv.rename(owt='owt_i')
    b = conv.rename(owt='owt_j')
    return np.degrees(sam(a, b, 'band'))


def sample_members(mean, sd, n=200, seed=0, slope=True):
    r"""Synthetic members of each optical water type.

    :math:`R_{rs} = \overline{R_{rs}} + \sigma\, \left(z_1 + z_2\, (\lambda - 600)/200\right)`,
    :math:`z_i \sim \mathcal{N}(0, 0.5)`, clipped at :math:`10^{-5}` sr\ :sup:`-1`.

    The spectrally varying term changes the spectral *shape* (not only the amplitude), otherwise
    the spectral angle would be trivially invariant. This is a crude proxy, to be replaced by the
    Bi et al. (2024) training spectra or by in situ HYPERNETS/AERONET-OC spectra.

    :param mean: mean spectra (``owt``, ``wl``)
    :param sd: standard deviation (``owt``, ``wl``)
    :param n: number of members per type
    :param seed: seed of the random generator
    :param slope: if False, :math:`z_2 = 0` (amplitude perturbations only)
    :return: :class:`xarray.DataArray` (``true_owt``, ``member``, ``wl``)
    """
    rng = np.random.default_rng(seed)
    z1 = xr.DataArray(rng.normal(0, .5, (mean.sizes['owt'], n)), dims=('owt', 'member'))
    z2 = xr.DataArray(rng.normal(0, .5 if slope else 0, (mean.sizes['owt'], n)), dims=('owt', 'member'))
    x = (mean.wl - 600.) / 200.
    return (mean + sd * (z1 + z2 * x)).clip(min=1e-5).rename(owt='true_owt')


def benchmark(sensors, database='Bi2024', n=200, noise=None, seed=0):
    """Agreement of the optical water type of each sensor with the hyperspectral type.

    The synthetic members (:func:`sample_members`) are classified on their hyperspectral
    spectrum (400-900 nm, the "truth") and on the bands of each sensor.

    :param sensors: dict name -> :class:`~trishna_tbx.sensor.Sensor`
    :param database: reference set, see :func:`load_reference`
    :param n: number of members per type
    :param noise: optional dict name -> Rrs noise (1 sigma, sr-1), :class:`xarray.DataArray` on
        ``band``, added to the band values
    :param seed: seed of the random generator
    :return: (:class:`pandas.DataFrame` of the fraction of agreement per true type and sensor,
        with a row ``'all'``; dict name -> confusion matrix normalized by row)
    """
    mean, sd = load_reference(database)
    members = sample_members(mean, sd, n=n, seed=seed)
    ref_hyp = mean.sel(wl=slice(400, 900))
    hyp_cls, _ = classify(members.sel(wl=slice(400, 900)), ref_hyp, dim='wl')
    rng = np.random.default_rng(seed + 1)
    rows, conf = {}, {}
    for name, sens in sensors.items():
        ref = convolve_reference(mean, sens)
        x = convolve_reference(members.rename(true_owt='owt'), sens).rename(owt='true_owt')
        if noise and name in noise:
            eps = xr.DataArray(rng.standard_normal(x.shape), dims=x.dims, coords=x.coords)
            sigma = noise[name].sel(band=x.band.values).drop_vars('wl', errors='ignore')
            x = x + sigma * eps
        cls, _ = classify(x, ref, dim='band')
        ok = (cls == hyp_cls)
        rows[name] = ok.mean('member').to_pandas()
        conf[name] = pd.crosstab(hyp_cls.values.ravel(), cls.values.ravel(),
                                 rownames=['hyperspectral'], colnames=[name], normalize='index')
    df = pd.DataFrame(rows)
    df.loc['all'] = df.mean()
    return df, conf
