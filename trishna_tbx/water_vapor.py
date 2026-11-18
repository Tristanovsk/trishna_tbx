"""Atmospheric water vapour (TCWV) from a window / absorption band pair.

TRISHNA: NIR 860 nm (window) / WV 910 nm (absorption)
S2:      B8A 865 nm (window) / B9 945 nm (absorption)

Ratio model (decoupled, as in Sen2Cor/MAJA-type "APDA" retrievals):

    R_abs / R_win = c * Tw_abs(M * W) / Tw_win(M * W)

with M = 1/cos(sza) + 1/cos(vza) the geometric air mass, W the TCWV (kg m-2) and c the
surface continuum ratio (c = 1 for a flat surface spectrum, the S2 assumption).

Over water the surface is dark and the path radiance dominates, which makes the
decoupled model underestimate W by a factor 3-5 (see hgrs.hgrs_kernel.WaterVapor docs).
For water pixels use ``coupled=True``, which uses the hgrs/radcalnet_oc path-length
gas coupling (needs AOT and a path reflectance first guess), or retrieve W on land/bright
pixels and interpolate over water (recommended operational approach).

Gas absorption: REPTRAN-fine H2O optical thickness per kg m-2 from radcalnet_oc
(``data/lut/gases/lut_abs_opt_thickness_normalized.nc``).
"""
from functools import lru_cache
from pathlib import Path

import numpy as np
import xarray as xr

from .sensor import TRISHNA, make_s2, solar_irradiance

#: Water vapour nodes of the ratio look-up table (kg m-2, i.e. 0.1 g cm-2 per unit).
W_GRID = np.linspace(0., 80., 161)
#: Air mass nodes of the ratio look-up table.
M_GRID = np.linspace(2., 8., 25)
#: Window and absorption bands of each sensor.
PAIRS = {'TRISHNA': ('NIR', 'WV'), 'S2': ('B8A', 'B9')}


@lru_cache
def h2o_optical_thickness():
    """Water vapour optical thickness per kg m-2 (REPTRAN fine) of radcalnet_oc (or grs).

    :return: :class:`xarray.DataArray` on ``wl`` (nm)
    """
    try:
        import radcalnet_oc
        f = Path(radcalnet_oc.__file__).parent / 'data/lut/gases/lut_abs_opt_thickness_normalized.nc'
    except ImportError:  # same table is shipped with grs
        import grs
        f = Path(grs.__file__).parent / 'data/lut/gases/lut_abs_opt_thickness_normalized.nc'
    return xr.open_dataset(f)['h2o'].load()


def band_transmittance(sensor, bands, w=W_GRID, m=M_GRID):
    r"""Band-averaged water vapour transmittance, weighted by the SRF and the solar irradiance.

    :math:`T_k(W, M) = \int e^{-\tau(\lambda) W M} E_0 S_k\, d\lambda \,/\, \int E_0 S_k\, d\lambda`,
    computed on the native REPTRAN grid between 800 and 1000 nm.

    :param sensor: :class:`~trishna_tbx.sensor.Sensor`
    :param bands: names of the bands
    :param w: total column water vapour (kg m-2)
    :param m: geometric air mass
    :return: :class:`xarray.DataArray` (``band``, ``W``, ``M``)
    """
    tau = h2o_optical_thickness()
    wl = tau.wl.values
    sel = (wl > 800) & (wl < 1000)
    tau = tau[sel]
    srf = sensor.srf.sel(band=list(bands)).interp(wl=tau.wl.values).fillna(0.)
    srf = srf * solar_irradiance().interp(wl=tau.wl.values).values
    W = xr.DataArray(w, dims='W', coords={'W': w})
    M = xr.DataArray(m, dims='M', coords={'M': m})
    T = np.exp(-tau * W * M)
    return (T * srf).sum('wl') / srf.sum('wl')


class RatioRetrieval:
    """Two-band ratio retrieval of the total column water vapour.

    The look-up table of the ratio :math:`T_{abs}/T_{win}` is computed on :data:`W_GRID` and
    :data:`M_GRID` with :func:`band_transmittance`.

    :param sensor: ``'TRISHNA'`` (860/910 nm) or ``'S2'`` (865/945 nm)

    :ivar ratio_lut: ratio look-up table (``W``, ``M``)
    :ivar T: band transmittances (``band``, ``W``, ``M``)
    """

    def __init__(self, sensor='TRISHNA'):
        self.name = sensor
        sens = TRISHNA if sensor == 'TRISHNA' else make_s2()
        self.win, self.abs = PAIRS[sensor]
        T = band_transmittance(sens, [self.win, self.abs])
        self.ratio_lut = (T.sel(band=self.abs) / T.sel(band=self.win)).drop_vars(['band', 'wl'], errors='ignore')
        self.T = T

    def forward(self, w, m, c=1.):
        """Ratio of the absorption band to the window band.

        :param w: total column water vapour (kg m-2)
        :param m: geometric air mass
        :param c: surface continuum ratio (1 for a flat surface spectrum)
        :return: ratio
        """
        return c * self.ratio_lut.interp(W=w, M=m)

    def retrieve(self, r_win, r_abs, sza, vza, c=1.):
        """Water vapour of each pixel from the band ratio.

        The ratio curve is interpolated at the air mass of the pixel, then inverted (it decreases
        monotonically with W). Ratios out of the table give NaN.

        :param r_win: reflectance of the window band (860 or 865 nm)
        :param r_abs: reflectance of the absorption band (910 or 945 nm)
        :param sza: solar zenith angle (deg)
        :param vza: viewing zenith angle (deg)
        :param c: surface continuum ratio
        :return: total column water vapour (kg m-2), :class:`numpy.ndarray`
        """
        m = 1. / np.cos(np.radians(sza)) + 1. / np.cos(np.radians(vza))
        ratio = np.asarray(r_abs / r_win / c)
        m = np.broadcast_to(np.asarray(m), ratio.shape)
        lut = self.ratio_lut.transpose('M', 'W').values
        out = np.full(ratio.shape, np.nan)
        for i in np.ndindex(ratio.shape):
            if not np.isfinite(ratio[i]):
                continue
            j = np.clip(np.searchsorted(M_GRID, m[i]), 1, len(M_GRID) - 1)
            f = (m[i] - M_GRID[j - 1]) / (M_GRID[j] - M_GRID[j - 1])
            curve = (1 - f) * lut[j - 1] + f * lut[j]
            out[i] = np.interp(ratio[i], curve[::-1], W_GRID[::-1], left=np.nan, right=np.nan)
        return out

    def sensitivity(self, m=3.):
        """Sensitivity of the ratio to the water vapour.

        :param m: geometric air mass
        :return: :class:`xarray.Dataset` on ``W``: ratio, its derivative ``dratio_dW`` and the
            error on W (kg m-2) for a 1 % error on the ratio
        """
        r = self.ratio_lut.sel(M=m, method='nearest')
        dr = r.differentiate('W')
        return xr.Dataset({'ratio': r, 'dratio_dW': dr, 'dW_for_1pct_ratio_error': np.abs(0.01 * r / dr)})


def compare_trishna_s2(m=3.):
    """Sensitivity to the water vapour of the TRISHNA 860/910 and S2 865/945 nm ratios.

    :param m: geometric air mass
    :return: :class:`xarray.Dataset` of :meth:`RatioRetrieval.sensitivity`, dimension ``sensor``
    """
    return xr.concat([RatioRetrieval(s).sensitivity(m) for s in ('TRISHNA', 'S2')],
                     dim=xr.DataArray(['TRISHNA', 'S2'], dims='sensor'))


def noise_error(sensor='TRISHNA', rho_surface=0.3, w=W_GRID, m=3., sza=30., n_aggregated=1):
    r"""Error (1 sigma) on the water vapour due to the radiometric noise only.

    :math:`\sigma_r / r = \sqrt{(\sigma_{abs}/R_{abs})^2 + (\sigma_{win}/R_{win})^2}`,
    :math:`\sigma_W = \sigma_r \,/\, |\partial r / \partial W|`, for a flat surface reflectance
    and without the path radiance (valid over land and bright targets, not over dark water).

    :param sensor: ``'TRISHNA'`` or ``'S2'``
    :param rho_surface: surface reflectance
    :param w: total column water vapour (kg m-2)
    :param m: geometric air mass
    :param sza: solar zenith angle (deg)
    :param n_aggregated: number of native pixels averaged per output pixel
    :return: error on W (kg m-2), :class:`xarray.DataArray` on ``W``
    """
    from .sensor import TRISHNA, make_s2
    rr = RatioRetrieval(sensor)
    sens = TRISHNA if sensor == 'TRISHNA' else make_s2()
    nz = sens.noise_reflectance(sza=sza, n_aggregated=n_aggregated).assign_coords(band=sens.band_names)
    T = rr.T.sel(M=m, method='nearest').interp(W=w)
    r_abs, r_win = rho_surface * T.sel(band=rr.abs), rho_surface * T.sel(band=rr.win)
    rel = np.sqrt((nz.sel(band=rr.abs) / r_abs) ** 2 + (nz.sel(band=rr.win) / r_win) ** 2)
    ratio = rr.ratio_lut.sel(M=m, method='nearest').interp(W=w)
    return (rel * ratio / np.abs(ratio.differentiate('W'))).drop_vars('band', errors='ignore')
