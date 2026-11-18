"""Semi-analytical Rrs inversion (invRrs) for TRISHNA bands.

invRrs integrates the forward Rrs over the sensor RSR (E_d weighted) but only knows the
Py6S sensors (oli, s2a, s2b, s3a/b, meris, modis). ``trishna_forward_model`` injects the
TRISHNA (or any ``trishna_tbx.sensor.Sensor``) responses without modifying invRrs.

With 4 water bands (485/555/670/860) only ``MultispectralSimple`` (chl, a_cdom, spm_norm,
offset) is well posed; even then chl and a_cdom are strongly correlated without 443 nm.
Use ``information_content`` to quantify that from the Jacobian.
"""
import numpy as np
import xarray as xr

from .sensor import TRISHNA

#: TRISHNA bands used for the inversion of the water reflectance.
WATER_BANDS = ['Blue', 'Green', 'Red', 'NIR']


class SensorRSR:
    """Spectral responses of a :class:`~trishna_tbx.sensor.Sensor` in the format of invRrs.

    Same interface as ``invRrs.spectral_integration.spectral_responses``, without Py6S.

    :param sensor: :class:`~trishna_tbx.sensor.Sensor`
    :param bands: names of the bands; default all
    """

    def __init__(self, sensor, bands=None):
        from invRrs.spectral_integration import spectral_responses
        self._impl = spectral_responses.__new__(spectral_responses)  # skip Py6S loading
        sens = sensor.subset(bands) if bands else sensor
        self._impl.rsrs = np.array([sens.srf.sel(band=b).where(sens.srf.sel(band=b) > 1e-4, drop=True)
                                    .drop_vars('band') for b in sens.band_names], dtype=object)
        self.rsrs = self._impl.rsrs

    def convolution(self, wl, param):
        """Band-average ``param`` given on ``wl`` (method of invRrs)."""
        return self._impl.convolution(wl, param)

    def xr_convolution(self, xarr, wl_sat):
        """Band-average an :class:`xarray.DataArray` (method of invRrs)."""
        return self._impl.xr_convolution(xarr, wl_sat)


def trishna_forward_model(model='MultispectralSimple', sensor=TRISHNA, bands=WATER_BANDS, sza=30., **kw):
    """invRrs inversion class on the TRISHNA bands (or the bands of any sensor).

    The spectral responses of ``sensor`` are injected in ``invRrs.forward`` while the class is
    built (invRrs only knows the Py6S sensors); invRrs is not modified.

    :param model: inversion class of ``invRrs.inverse``, e.g. ``'MultispectralSimple'``
    :param sensor: :class:`~trishna_tbx.sensor.Sensor`
    :param bands: names of the bands used
    :param sza: solar zenith angle (deg)
    :param kw: other arguments of the class (``aphy_star``, ``bbphy_star``...)
    :return: instance of the invRrs class
    """
    import invRrs.forward as fwd
    import invRrs.inverse as inv
    cls = getattr(inv, model)
    orig = fwd.spectral
    fwd.spectral = lambda sat, band_idx=None: SensorRSR(sensor, bands)
    try:
        aphy_star, bbphy_star = kw.pop('aphy_star', None), kw.pop('bbphy_star', None)
        return cls(aphy_star, bbphy_star, sza, sat='trishna', **kw)
    finally:
        fwd.spectral = orig


def information_content(model, x, sigma):
    r"""Information content of the inversion: posterior uncertainty of the unknowns.

    Linearized, with a flat prior: :math:`S = (J^T S_\epsilon^{-1} J)^{-1}`, with :math:`J` the
    Jacobian of the forward model at ``x`` (finite differences) and :math:`S_\epsilon` the
    diagonal noise covariance. A large standard deviation, or a correlation close to ±1, means
    that the unknown is not retrievable.

    :param model: invRrs instance, e.g. from :func:`trishna_forward_model`
    :param x: state of the unknowns (``model.param_names``)
    :param sigma: Rrs noise of each band (sr-1)
    :return: :class:`xarray.Dataset`: standard deviation ``std`` (``param``) and correlation
        matrix ``corr`` (``param``, ``param2``)
    """
    from scipy.optimize import approx_fprime
    f = lambda p: np.asarray(model.forward_model(np.asarray(p, float), output='Rrs'))
    J = np.array([approx_fprime(np.asarray(x, float), lambda p, i=i: f(p)[i], 1e-6 * np.abs(x) + 1e-9)
                  for i in range(len(model.wl))])
    S = np.linalg.pinv(J.T @ np.diag(1 / np.asarray(sigma) ** 2) @ J)
    std = np.sqrt(np.diag(S))
    return xr.Dataset({'std': ('param', std), 'corr': (('param', 'param2'), S / np.outer(std, std))},
                      coords={'param': list(model.param_names), 'param2': list(model.param_names)})
