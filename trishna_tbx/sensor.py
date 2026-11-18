"""Sensor definitions: TRISHNA VSWIR, Sentinel-2 MSI, and band-restricted S2 subsets.

TRISHNA VSWIR specification from Table 6 of
"The TRISHNA Mission: A New Era of High-Resolution Thermal..." (EGUsphere-2026-4492):

    band    centre  FWHM  Lmin  Lmax   SNR (specified at Lmin)
    Blue      485    70    30   663.5   50
    Green     555    70    30   600.5  100
    Red       670    60    20   486.1  100
    NIR       860    40    30   325.3  100
    WV        910    30    30   277     50
    Cirrus   1380    30     6    55.0   25
    SWIR     1610   100     4    77.9  100     (radiances in W m-2 sr-1 um-1)

The real TRISHNA spectral responses are not public yet: we use a super-Gaussian
(near rectangular, the shape used in the mission study) whose exponent can be changed.
Replace ``TRISHNA.srf`` with the measured SRF as soon as ISRO delivers it.
"""
from dataclasses import dataclass, field
from functools import lru_cache
from importlib import resources
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

#: Common spectral grid of the spectral responses (nm), 350-2500 nm by 1 nm.
WL = np.arange(350., 2500.5, 1.)


@dataclass
class Band:
    """Spectral band of a sensor.

    :param name: band name (``'Blue'``, ``'B8A'``...)
    :param center: central wavelength (nm)
    :param fwhm: full width at half maximum (nm)
    :param snr: signal-to-noise ratio specified at ``lref``
    :param lref: radiance at which ``snr`` is specified (W m-2 sr-1 um-1): Lmin for TRISHNA,
        Lref for Sentinel-2
    :param s2_proxy: closest Sentinel-2 band, used to emulate the band from S2 images
    :param role: main use of the band over water
    """
    name: str
    center: float  # nm
    fwhm: float  # nm
    snr: float = np.nan
    lref: float = np.nan  # radiance at which SNR is specified (W m-2 sr-1 um-1)
    s2_proxy: str = None  # closest S2 band used for emulation
    role: str = ''


#: TRISHNA VSWIR bands (Table 6 of the mission paper) and their Sentinel-2 proxies.
TRISHNA_BANDS = [
    Band('Blue', 485, 70, 50, 30, 'B2', 'water colour, AOT'),
    Band('Green', 555, 70, 100, 30, 'B3', 'water colour'),
    Band('Red', 670, 60, 100, 20, 'B4', 'SPM/turbidity, vegetation'),
    Band('NIR', 860, 40, 100, 30, 'B8A', 'turbid waters, glint/aerosol, WV reference'),
    Band('WV', 910, 30, 50, 30, 'B9', 'water vapour absorption (S2 B9 is at 945 nm!)'),
    Band('Cirrus', 1380, 30, 25, 6, 'B10', 'cirrus detection'),
    Band('SWIR', 1610, 100, 100, 4, 'B11', 'aerosol/glint over water, masking'),
]

#: Sentinel-2 MSI bands (S2A central wavelengths, L1C SNR at Lref of the MSI requirement document).
S2_BANDS = [
    Band('B1', 443, 21, 129, 129), Band('B2', 490, 66, 154, 128), Band('B3', 560, 36, 168, 128),
    Band('B4', 665, 31, 142, 108), Band('B5', 705, 15, 117, 74.5), Band('B6', 740, 15, 89, 68),
    Band('B7', 783, 20, 105, 67), Band('B8', 842, 106, 174, 103), Band('B8A', 865, 21, 72, 52.5),
    Band('B9', 945, 20, 114, 9), Band('B10', 1375, 31, 50, 6), Band('B11', 1610, 91, 100, 4),
    Band('B12', 2190, 175, 100, 1.5),
]

#: Sentinel-2 band -> TRISHNA band, the bands kept to emulate TRISHNA from Sentinel-2.
S2_TO_TRISHNA = {b.s2_proxy: b.name for b in TRISHNA_BANDS}


def super_gaussian(wl, center, fwhm, expon=2.):
    """Super-Gaussian spectral response, normalized to 1 at its maximum.

    :param wl: wavelengths (nm)
    :param center: central wavelength (nm)
    :param fwhm: full width at half maximum (nm)
    :param expon: exponent; 1 gives a Gaussian, large values a boxcar
    :return: response on ``wl``
    """
    sigma = fwhm / (2. * (2. * np.log(2.)) ** (1. / (2. * expon)))
    return np.exp(-0.5 * np.abs((wl - center) / sigma) ** (2 * expon))


@dataclass
class Sensor:
    """Multispectral sensor: list of bands and their spectral response functions.

    :param name: sensor name
    :param bands: list of :class:`Band`
    :param srf: spectral response functions, :class:`xarray.DataArray` of dimensions
        (``band``, ``wl``), wavelengths in nm, maximum normalized to 1
    """
    name: str
    bands: list
    srf: xr.DataArray = field(default=None, repr=False)  # dims (band, wl)

    @property
    def band_names(self):
        """Names of the bands."""
        return [b.name for b in self.bands]

    @property
    def centers(self):
        """Central wavelengths (nm), :class:`xarray.DataArray` on ``band``."""
        return xr.DataArray([b.center for b in self.bands], dims='band', coords={'band': self.band_names})

    def subset(self, names, name=None):
        """Sensor restricted to some bands.

        :param names: names of the bands to keep (the order of the sensor is kept)
        :param name: name of the new sensor; default ``<name>_subset``
        :return: :class:`Sensor`
        """
        bands = [b for b in self.bands if b.name in names]
        return Sensor(name or f'{self.name}_subset', bands, self.srf.sel(band=[b.name for b in bands]))

    def convolve(self, spectrum, wl_dim='wl', weight=None):
        r"""Band-average a spectrum with the spectral response functions.

        :math:`X_k = \int X\, w\, S_k\, d\lambda \,/\, \int w\, S_k\, d\lambda`

        :param spectrum: :class:`xarray.DataArray` with a wavelength dimension (nm), any other
            dimensions
        :param wl_dim: name of the wavelength dimension
        :param weight: spectral weight :math:`w`, :class:`xarray.DataArray` on wavelength, e.g. the
            solar irradiance :func:`solar_irradiance` for reflectances; ``None`` (flat) for
            radiances or transmittances
        :return: :class:`xarray.DataArray` with ``band`` in place of ``wl_dim`` and the central
            wavelengths as coordinate ``wl``
        """
        srf = self.srf.interp(wl=spectrum[wl_dim].values).fillna(0.).rename(wl=wl_dim)
        if weight is not None:
            srf = srf * weight.interp({weight.dims[0]: spectrum[wl_dim].values}).values
        out = (spectrum * srf).sum(wl_dim) / srf.sum(wl_dim)
        return out.assign_coords(wl=('band', self.centers.values))

    def noise_reflectance(self, sza=30., d2=1., n_aggregated=1):
        r"""1-sigma noise of the TOA reflectance from the SNR specification.

        :math:`\sigma_\rho = \pi\, (L_{ref}/SNR) \,/\, (E_0 \cos\theta_s / d^2) \,/\, \sqrt{n}`

        S2 SNR is specified at Lref (bright) whereas TRISHNA SNR is specified at Lmin (dark): use
        the radiance-dependent model :meth:`noise_at_radiance` for a fair dark-water comparison.

        :param sza: solar zenith angle (deg)
        :param d2: squared Earth-Sun distance (AU²)
        :param n_aggregated: number of native pixels averaged per output pixel, assuming
            uncorrelated noise (S2 10 m -> 60 m: 36, 20 m -> 60 m: 9)
        :return: :class:`xarray.DataArray` on ``band``
        """
        e0 = self.convolve(solar_irradiance())  # W m-2 nm-1 -> x1000 per um
        lref = xr.DataArray([b.lref for b in self.bands], dims='band')
        snr = xr.DataArray([b.snr for b in self.bands], dims='band')
        n = xr.DataArray(np.broadcast_to(n_aggregated, len(self.bands)), dims='band')
        return np.pi * (lref / snr) / (e0 * 1e3 * np.cos(np.radians(sza)) / d2) / np.sqrt(n)

    def noise_at_radiance(self, L, read_fraction=0.5):
        r"""Radiance-dependent noise (read + shot noise).

        :math:`\sigma(L) = (L_{ref}/SNR)\, \sqrt{f + (1-f)\, L/L_{ref}}`

        :param L: radiance (W m-2 sr-1 um-1), :class:`xarray.DataArray` on ``band``
        :param read_fraction: fraction :math:`f` of the noise variance at :math:`L_{ref}` that does
            not depend on the signal (unknown before launch: sensitivity parameter)
        :return: radiance noise (W m-2 sr-1 um-1)
        """
        lref = xr.DataArray([b.lref for b in self.bands], dims='band', coords={'band': self.band_names})
        snr = xr.DataArray([b.snr for b in self.bands], dims='band', coords={'band': self.band_names})
        return (lref / snr) * np.sqrt(read_fraction + (1 - read_fraction) * L / lref)


@lru_cache
def solar_irradiance():
    """Extraterrestrial solar irradiance of Thuillier et al. (2003), from radcalnet_oc.

    Without radcalnet_oc, a 5778 K black body is used (approximate).

    :return: :class:`xarray.DataArray` (W m-2 nm-1) on the 1 nm grid :data:`WL` (350-2500 nm)
    """
    try:
        import radcalnet_oc
        f = Path(radcalnet_oc.__file__).parent / 'data/auxdata/ref_atlas_thuillier3.nc'
        ds = xr.open_dataset(f)
        e0 = ds['data'].squeeze().rename(wavelength='wl')
        e0 = e0.where(e0 > 0).dropna('wl').drop_vars('time', errors='ignore')
        return e0.interp(wl=WL)
    except Exception:  # pragma: no cover - approximate fallback (Planck 5778 K)
        h, c, k = 6.626e-34, 2.998e8, 1.381e-23
        lam = WL * 1e-9
        b = 2 * h * c ** 2 / lam ** 5 / (np.exp(h * c / (lam * k * 5778)) - 1)
        return xr.DataArray(b * 6.8e-5 * 1e-9, dims='wl', coords={'wl': WL})


def make_trishna(expon=4., bands=TRISHNA_BANDS):
    """TRISHNA VSWIR with super-Gaussian spectral responses (measured ones not public yet).

    :param expon: exponent of the super-Gaussian (:func:`super_gaussian`)
    :param bands: list of :class:`Band`, default :data:`TRISHNA_BANDS`
    :return: :class:`Sensor`
    """
    srf = np.stack([super_gaussian(WL, b.center, b.fwhm, expon) for b in bands])
    srf = xr.DataArray(srf, dims=('band', 'wl'), coords={'band': [b.name for b in bands], 'wl': WL})
    return Sensor('TRISHNA', list(bands), srf)


@lru_cache
def _s2_srf_table(platform='S2A'):
    import GRSdriver
    f = Path(GRSdriver.__file__).parent / 'rsr/data/Sentinel-2_MSI_Spectral_Responses.xlsx'
    df = pd.read_excel(f, sheet_name=f'Spectral Responses ({platform})').set_index('SR_WL')
    df.columns = [c.split('_')[-1] for c in df.columns]
    return df


def make_s2(platform='S2A', bands=S2_BANDS):
    """Sentinel-2 MSI with the ESA spectral responses (file of GRSdriver).

    Without GRSdriver, the responses are approximated by Gaussians.

    :param platform: ``'S2A'``, ``'S2B'`` (or ``'S2C'`` if in the ESA file)
    :param bands: list of :class:`Band`, default :data:`S2_BANDS`
    :return: :class:`Sensor`
    """
    try:
        df = _s2_srf_table(platform)
        srf = xr.DataArray(df[[b.name for b in bands]].values.T, dims=('band', 'wl'),
                           coords={'band': [b.name for b in bands], 'wl': df.index.values.astype(float)})
        srf = srf.interp(wl=WL).fillna(0.)
    except Exception:  # GRSdriver not installed: Gaussian approximation
        srf = xr.DataArray(np.stack([super_gaussian(WL, b.center, b.fwhm, 1.) for b in bands]),
                           dims=('band', 'wl'), coords={'band': [b.name for b in bands], 'wl': WL})
    return Sensor(f'Sentinel-2 {platform[-1]}', list(bands), srf)


def make_s2_trishna_like(platform='S2A'):
    """Sentinel-2 restricted to the bands that have a TRISHNA counterpart.

    These are the bands kept in an emulated TRISHNA L1C or L2A (B2, B3, B4, B8A, B9, B10, B11).

    :param platform: Sentinel-2 platform
    :return: :class:`Sensor` named ``'S2-as-TRISHNA'``
    """
    s2 = make_s2(platform)
    return s2.subset(list(S2_TO_TRISHNA), name='S2-as-TRISHNA')


#: TRISHNA VSWIR sensor (super-Gaussian responses, exponent 4), see :func:`make_trishna`.
TRISHNA = make_trishna()
