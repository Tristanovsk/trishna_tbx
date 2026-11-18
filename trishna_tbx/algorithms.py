"""Water-quality algorithm registry with explicit band requirements.

Each entry declares the wavelengths it needs. ``coverage(sensor)`` tells whether the sensor
has a band within ``tol`` nm of each, which avoids the silent ``sel(method='nearest')``
substitution of GRSl2bgen (on TRISHNA 443 -> 485 and 705 -> 670 would be used without
warning).

Input convention: ``Rrs`` DataArray (sr-1) with a ``wl`` dimension in nm (GRS L2A).
Coefficients are the published ones for S2/MERIS/OCTS bands: for TRISHNA they must be
re-tuned (see ``retune`` and ``trishna_tbx.simulate``), because the TRISHNA bands are
2x wider than S2 B3/B4 and shifted (485 vs 490, 670 vs 665, 860 vs 865).
"""
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
import xarray as xr

NEAREST = dict(method='nearest')


def _b(Rrs, wl):
    return Rrs.sel(wl=wl, **NEAREST, drop=True)


@dataclass
class Algo:
    """Water-quality algorithm and the wavelengths it needs.

    :param name: algorithm name (``'OC2'``, ``'Nechad10'``...)
    :param product: retrieved product (``'Chla'``, ``'SPM'``, ``'Kd490'``...)
    :param wls: wavelengths needed (nm)
    :param func: function of ``Rrs`` (None when not implemented: listed for the coverage table)
    :param reference: publication
    :param water: water types for which the algorithm is designed
    """
    name: str
    product: str
    wls: list
    func: callable = field(repr=False)
    reference: str = ''
    water: str = ''

    def run(self, Rrs, sensor=None, strict=True, tol=15.):
        """Run the algorithm.

        :param Rrs: :class:`xarray.DataArray` (sr-1) with a ``wl`` dimension (nm)
        :param sensor: :class:`~trishna_tbx.sensor.Sensor` of ``Rrs``
        :param strict: raise :class:`ValueError` if the sensor has no band within ``tol`` of a
            needed wavelength (instead of using the nearest band without warning)
        :param tol: tolerance (nm)
        :return: :class:`xarray.DataArray` named ``<product>_<name>``
        """
        if strict and sensor is not None and not self.covered(sensor, tol):
            raise ValueError(f'{self.name} needs {self.wls} nm, not available on {sensor.name}')
        return self.func(Rrs).rename(f'{self.product}_{self.name}')

    def covered(self, sensor, tol=15.):
        """True if the sensor has a band within ``tol`` nm of each needed wavelength."""
        c = np.array([b.center for b in sensor.bands])
        return all(np.min(np.abs(c - w)) <= tol for w in self.wls)


# ---------------------------------------------------------------- implementations
def ocx(ratio, coef):
    r"""OCx polynomial, :math:`10^{\sum_i c_i\, x^i}` with :math:`x` the log10 of a band ratio."""
    return 10 ** sum(c * ratio ** i for i, c in enumerate(coef))


def oc2(Rrs, coef=(0.2236, -1.8296, 1.9094, -2.9481, -0.1718)):  # NASA OC2 (OCTS), GRSl2bgen
    """Chlorophyll-a (mg m-3) with the blue/green ratio Rrs(490)/Rrs(560) (O'Reilly et al., 2000).

    Coefficients of OCTS (GRSl2bgen); NaN outside the valid range 0.001-100 mg m-3.
    """
    chl = ocx(np.log10(_b(Rrs, 490) / _b(Rrs, 560)), coef)
    return chl.where((chl > 0.001) & (chl < 100.))  # NASA valid range


def kd490_kd2(Rrs, coef=(-0.8641, -1.6549, 2.0112, -2.5174, -1.1035), kw=0.0166):  # KD2E
    r"""Diffuse attenuation coefficient at 490 nm (m-1), KD2 with the MERIS KD2E coefficients.

    :math:`K_d(490) = K_w + \mathrm{OCx}(R_{rs}(490)/R_{rs}(560))`; NaN beyond 6.4 m-1 (NASA
    valid range, the polynomial diverges beyond).
    """
    kd = kw + ocx(np.log10(_b(Rrs, 490) / _b(Rrs, 560)), coef)
    return kd.where(kd <= 6.4)  # NASA KD2 valid range (polynomial diverges beyond)


def kdpar_rd22(Rrs, a=(3.09, -90.17)):  # Roy & Das 2022
    r"""Diffuse attenuation coefficient of PAR (m-1), :math:`a_0\, e^{a_1 (R_{rs}(490) - R_{rs}(665))}`
    (Roy & Das, 2022, inland waters)."""
    return a[0] * np.exp(a[1] * (_b(Rrs, 490) - _b(Rrs, 665)))


def spm_nechad10(Rrs, wl=665, A=355.85, C=0.1728):  # Nechad et al. 2010, rho_w = pi Rrs
    r"""Suspended particulate matter (g m-3), Nechad et al. (2010) at 665 nm.

    :math:`SPM = A \rho_w / (1 - \rho_w / C)`, :math:`\rho_w = \pi R_{rs}`; NaN where
    :math:`\rho_w / C \geq 0.9` (saturation).
    """
    r = np.pi * _b(Rrs, wl)
    return (A * r / (1 - r / C)).where(r / C < 0.9)  # diverges when rho_w -> C (saturation)


def turb_dogliotti15(Rrs, low=0.05, high=0.07):  # Dogliotti et al. 2015 (645 & 859 nm)
    r"""Turbidity (FNU), Dogliotti et al. (2015): red algorithm, NIR algorithm in turbid waters.

    The red (665 nm) and NIR (865 nm) estimates are blended linearly for :math:`\rho_w(665)`
    between ``low`` and ``high``.
    """
    r_red, r_nir = np.pi * _b(Rrs, 665), np.pi * _b(Rrs, 865)
    t_red = (228.1 * r_red / (1 - r_red / 0.1641)).where(r_red / 0.1641 < 0.9)
    t_nir = (3078.9 * r_nir / (1 - r_nir / 0.2112)).where(r_nir / 0.2112 < 0.9)
    w = ((r_red - low) / (high - low)).clip(0, 1)
    blend = (1 - w) * t_red + w * t_nir
    return xr.where(w <= 0, t_red, xr.where(w >= 1, t_nir, blend))  # a NaN branch only where used


def ndci(Rrs):  # Mishra & Mishra 2012 (index only)
    """Normalized difference chlorophyll index (705, 665 nm), Mishra & Mishra (2012); not on TRISHNA."""
    return (_b(Rrs, 705) - _b(Rrs, 665)) / (_b(Rrs, 705) + _b(Rrs, 665))


def green_red_index(Rrs):
    """(G-R)/(G+R): Bi et al. (2024) NDI, TRISHNA-compatible proxy for turbidity/OWT."""
    return (_b(Rrs, 560) - _b(Rrs, 665)) / (_b(Rrs, 560) + _b(Rrs, 665))


def red_nir_ratio(Rrs):
    """R(860)/R(665): high-SPM indicator, saturation-free branch of SPM in very turbid waters."""
    return _b(Rrs, 865) / _b(Rrs, 665)


#: Algorithms and the wavelengths they need (those with ``func`` None are listed for the coverage only).
REGISTRY = [
    Algo('OC2', 'Chla', [490, 560], oc2, 'O\'Reilly et al. 2000 (OCTS coefs)', 'clear/moderate'),
    Algo('OC3', 'Chla', [443, 490, 560], None, 'O\'Reilly et al. 2000', 'clear'),
    Algo('NDCI', 'Chla_index', [665, 705], ndci, 'Mishra & Mishra 2012', 'productive'),
    Algo('Gons', 'Chla', [665, 705, 783], None, 'Gons et al. 2005', 'productive'),
    Algo('Gilerson2', 'Chla', [665, 705], None, 'Gilerson et al. 2010', 'productive'),
    Algo('Nechad10', 'SPM', [665], spm_nechad10, 'Nechad et al. 2010', 'moderate turbid'),
    Algo('Dogliotti15', 'TUR', [665, 865], turb_dogliotti15, 'Dogliotti et al. 2015', 'all turbid'),
    Algo('RedNIR', 'SPM_index', [665, 865], red_nir_ratio, '', 'very turbid'),
    Algo('Binding10', 'SPM', [740], None, 'Binding et al. 2010', 'very turbid'),
    Algo('KD2', 'Kd490', [490, 560], kd490_kd2, 'NASA KD2 / MERIS KD2E coefs', 'clear'),
    Algo('RD22', 'Kd_par', [490, 665], kdpar_rd22, 'Roy & Das 2022', 'inland'),
    Algo('QAAv6', 'IOPs/Zsd', [443, 490, 560, 665], None, 'Lee et al. 2002/2015', 'all'),
    Algo('QAA-485', 'IOPs/Zsd', [490, 560, 665], None, 'QAA with 485 replacing 443 (to develop)', 'all'),
    Algo('Brezonik15', 'aCDOM', [490, 740], None, 'Brezonik et al. 2015', 'inland'),
    Algo('NDI', 'OWT_index', [560, 665], green_red_index, 'Bi et al. 2024', 'all'),
]


def coverage(sensors, tol=15.):
    """Algorithms that each sensor can run.

    :param sensors: list of :class:`~trishna_tbx.sensor.Sensor`
    :param tol: tolerance on the wavelengths (nm)
    :return: :class:`pandas.DataFrame`, index (product, algorithm), the wavelengths needed and one
        boolean column per sensor
    """
    rows = {}
    for a in REGISTRY:
        rows[(a.product, a.name)] = {'needs (nm)': ', '.join(map(str, a.wls)),
                                    **{s.name: a.covered(s, tol) for s in sensors}}
    return pd.DataFrame(rows).T


def run_all(Rrs, sensor, tol=15.):
    """Run every implemented algorithm of :data:`REGISTRY` that the sensor covers.

    :param Rrs: :class:`xarray.DataArray` (sr-1) with a ``wl`` dimension (nm)
    :param sensor: :class:`~trishna_tbx.sensor.Sensor` of ``Rrs``
    :param tol: tolerance on the wavelengths (nm)
    :return: :class:`xarray.Dataset`, one variable per algorithm
    """
    out = {}
    for a in REGISTRY:
        if a.func is not None and a.covered(sensor, tol):
            da = a.run(Rrs, sensor, tol=tol)
            out[da.name] = da
    return xr.Dataset(out)


def retune_ocx(rrs_sensor, chl, blue, green, order=4):
    """Fit the coefficients of an OCx polynomial for a new pair of bands.

    :param rrs_sensor: band-averaged Rrs (simulated or in situ), :class:`xarray.DataArray` with a
        ``band`` dimension
    :param chl: chlorophyll-a concentration of each spectrum (mg m-3)
    :param blue: name of the blue band
    :param green: name of the green band
    :param order: order of the polynomial
    :return: coefficients, for :func:`ocx` (increasing powers)
    """
    x = np.log10(rrs_sensor.sel(band=blue) / rrs_sensor.sel(band=green)).values.ravel()
    y = np.log10(np.asarray(chl)).ravel()
    ok = np.isfinite(x) & np.isfinite(y)
    return np.polynomial.polynomial.polyfit(x[ok], y[ok], order)
