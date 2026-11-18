.. _methods:

Methods
=======

This page gives the equations of each step of the assessment, as implemented.

Notation
--------

.. list-table::
   :header-rows: 1
   :widths: 25 55 20

   * - Symbol
     - Definition
     - Unit
   * - :math:`\lambda`
     - wavelength
     - nm
   * - :math:`S_k(\lambda)`
     - spectral response of band :math:`k`, maximum 1
     - --
   * - :math:`E_0(\lambda)`
     - extraterrestrial solar irradiance (Thuillier et al., 2003)
     - W m\ :sup:`-2` nm\ :sup:`-1`
   * - :math:`R_{rs}`
     - remote-sensing reflectance
     - sr\ :sup:`-1`
   * - :math:`\rho_w = \pi R_{rs}`
     - water reflectance
     - --
   * - :math:`\theta_s,\ \theta_v`
     - solar and viewing zenith angles
     - deg
   * - :math:`M = 1/\cos\theta_s + 1/\cos\theta_v`
     - geometric air mass
     - --
   * - :math:`W`
     - total column water vapour (TCWV)
     - kg m\ :sup:`-2`

.. _methods-sensors:

Sensor models
-------------

The TRISHNA VSWIR bands (centre, width, SNR) are those of Table 6 of the mission paper. Their spectral
responses are not public yet: they are modelled by super-Gaussians,

.. math::
   :label: supergauss

   S_k(\lambda) = \exp\left(-\frac{1}{2}\left|\frac{\lambda - \lambda_k}{\sigma_k}\right|^{2p}\right),
   \qquad \sigma_k = \frac{\mathrm{FWHM}_k}{2\,(2\ln 2)^{1/(2p)}},

with :math:`p = 4` (nearly rectangular; :math:`p = 1` is a Gaussian). The Sentinel-2 responses are
the ESA ones (file of GRSdriver).

Band values are averaged with the solar irradiance as weight, which is the right definition for
reflectances (:meth:`Sensor.convolve <trishna_tbx.sensor.Sensor.convolve>`):

.. math::
   :label: convolve

   R_k = \frac{\int R(\lambda)\, E_0(\lambda)\, S_k(\lambda)\, d\lambda}{\int E_0(\lambda)\, S_k(\lambda)\, d\lambda}.

The noise of the TOA reflectance comes from the SNR specified at the radiance :math:`L_{ref}`
(Lmin for TRISHNA, Lref for Sentinel-2):

.. math::
   :label: noise

   \sigma_\rho = \frac{\pi\, L_{ref}/SNR}{E_{0,k} \cos\theta_s / d^2}\, \frac{1}{\sqrt{n}},
   \qquad \sigma(L) = \frac{L_{ref}}{SNR}\, \sqrt{f + (1 - f)\, L / L_{ref}},

with :math:`n` the number of native pixels averaged per output pixel (S2 20 m to 60 m: 9) and
:math:`f` the fraction of signal-independent noise of the radiance-dependent model. The two sensors
specify their SNR at different radiances: the dark-water comparison depends on :math:`f`, unknown
before launch.

.. _methods-emulation:

Emulation of TRISHNA from Sentinel-2
------------------------------------

**Band restriction.** The S2 bands that have a TRISHNA counterpart are kept (B2, B3, B4, B8A, B9, B10,
B11) and the image is block-averaged to 60 m (float variables averaged, masks by majority). In grs,
the processed bands are B2, B3, B4, B8A and B11: the aerosol and the sunglint are retrieved from the
last two processed bands, 865 and 1610 nm, as for TRISHNA (S2 uses 1610 and 2190 nm). B12 is loaded
for the ``high_swir`` mask of grs only (:func:`~trishna_tbx.emulate.grs_band_set`).

**Spectral band adjustment.** Band restriction does not correct the differences in centre and width
(TRISHNA Blue is 11 nm bluer than B2, Green and Red are twice wider). The adjustment factor of each
band is the median ratio over a hyperspectral library (:class:`~trishna_tbx.emulate.SBAF`),

.. math::
   :label: sbaf

   a_k = \operatorname{median}_i \left(\frac{R^{TRISHNA}_{k,i}}{R^{S2}_{k,i}}\right),
   \qquad R^{TRISHNA}_k \approx a_k\, R^{S2}_k,

fitted on the spectra above :math:`10^{-4}` sr\ :sup:`-1` (at least 5). The factors depend on the
water type (Blue: 1.28 in clear waters to 0.84 in very turbid waters), so
:func:`~trishna_tbx.emulate.owt_band_adjustment` fits them for each Bi et al. (2024) type and applies
to each pixel the factors of its type.

**Noise.** Gaussian white noise of standard deviation :eq:`noise` can be added to the TOA reflectance
(:func:`~trishna_tbx.emulate.add_noise`).

.. _methods-owt:

Optical water types
-------------------

The reference spectra (Bi et al., 2024: 10 types, 400-900 nm; Spyrakos et al., 2018: 13 types) are
**convolved** with the responses of each sensor :eq:`convolve`, rather than interpolated at the
band centres, which matters for the 60-70 nm wide TRISHNA bands. Each spectrum gets the type of
smallest spectral angle (:func:`~trishna_tbx.owt.sam`):

.. math::
   :label: sam

   \alpha_j = \arccos \frac{\sum_k x_k\, r_{j,k}}{\sqrt{\sum_k x_k^2 \sum_k r_{j,k}^2}},
   \qquad \mathrm{OWT} = \arg\min_j \alpha_j.

The separability of two types on a sensor is the spectral angle between their convolved reference
spectra (:func:`~trishna_tbx.owt.separability`). On the four TRISHNA water bands, the sub-types 3a/3b
(2.0°), 4a/4b (3.2°) and 1/2 (5.4°) are nearly collinear: they are merged into a 7-class TRISHNA
set (:data:`~trishna_tbx.owt.TRISHNA_MERGE`).

**Benchmark.** Synthetic members of each type,

.. math::
   :label: members

   R_{rs} = \overline{R_{rs}} + \sigma \left(z_1 + z_2\, \frac{\lambda - 600}{200}\right),
   \qquad z_i \sim \mathcal{N}(0, 0.5),

are classified on their hyperspectral spectrum (the reference) and on the bands of each sensor, with or
without noise (:func:`~trishna_tbx.owt.benchmark`). The spectral term :math:`z_2` changes the shape of
the spectrum, otherwise the spectral angle would be invariant.

.. _methods-algorithms:

Water-quality algorithms
------------------------

Each algorithm of :data:`~trishna_tbx.algorithms.REGISTRY` declares the wavelengths it needs. A sensor
covers an algorithm if it has a band within 15 nm of each (:func:`~trishna_tbx.algorithms.coverage`).
This avoids the silent nearest-band substitution (443 by 485 nm, 705 by 670 nm on TRISHNA). The
algorithms implemented, with the published coefficients:

.. list-table::
   :header-rows: 1
   :widths: 18 18 64

   * - Product
     - Algorithm
     - Equation
   * - Chl-a
     - OC2
     - :math:`\log_{10} \mathrm{Chl} = \sum_{i=0}^{4} c_i\, x^i`, :math:`x = \log_{10}(R_{rs}(490)/R_{rs}(560))`;
       valid 0.001-100 mg m\ :sup:`-3`
   * - :math:`K_d(490)`
     - KD2
     - :math:`K_d = K_w + 10^{\sum c_i x^i}`, same ratio, MERIS KD2E coefficients; valid up to 6.4 m\ :sup:`-1`
   * - :math:`K_d(\mathrm{PAR})`
     - Roy and Das 2022
     - :math:`3.09\, \exp(-90.17\,(R_{rs}(490) - R_{rs}(665)))`
   * - SPM
     - Nechad et al. 2010
     - :math:`A\, \rho_w / (1 - \rho_w/C)` at 665 nm, :math:`A = 355.85`, :math:`C = 0.1728`; valid for :math:`\rho_w/C < 0.9`
   * - Turbidity
     - Dogliotti et al. 2015
     - same form at 665 and 865 nm, blended linearly for :math:`0.05 < \rho_w(665) < 0.07`
   * - Indices
     - NDI, red/NIR, NDCI
     - :math:`(R_{560} - R_{665})/(R_{560} + R_{665})`, :math:`R_{865}/R_{665}`, :math:`(R_{705} - R_{665})/(R_{705} + R_{665})`

The coefficients of the band-ratio algorithms are those of the S2/MERIS/OCTS bands: for TRISHNA they
have to be re-fitted on simulated or in situ band values (:func:`~trishna_tbx.algorithms.retune_ocx`).

.. _methods-wv:

Water vapour
------------

The ratio of an absorption band to a window band (TRISHNA 910/860 nm, S2 945/865 nm) is modelled
without coupling between absorption and scattering:

.. math::
   :label: wvratio

   \frac{R_{abs}}{R_{win}} = c\, \frac{T_{abs}(W, M)}{T_{win}(W, M)},
   \qquad
   T_k(W, M) = \frac{\int e^{-\tau_{H_2O}(\lambda)\, W M}\, E_0 S_k\, d\lambda}{\int E_0 S_k\, d\lambda},

with :math:`\tau_{H_2O}` the water vapour optical thickness per kg m\ :sup:`-2` (REPTRAN fine, from
radcalnet_oc) and :math:`c` the continuum ratio of the surface (1 for a flat spectrum). The ratio is
tabulated in :math:`(W, M)` and inverted pixel by pixel at the pixel air mass
(:class:`~trishna_tbx.water_vapor.RatioRetrieval`). The error due to the noise alone is

.. math::
   :label: wvnoise

   \frac{\sigma_r}{r} = \sqrt{\left(\frac{\sigma_{abs}}{R_{abs}}\right)^2 + \left(\frac{\sigma_{win}}{R_{win}}\right)^2},
   \qquad \sigma_W = \frac{\sigma_r}{|\partial r / \partial W|}.

Over water the surface signal is too dark and the path radiance dominates: :eq:`wvratio` then
underestimates :math:`W`. The retrieval is meant for land pixels (then interpolated over water) or has
to be replaced by the coupled path-length model of hgrs and radcalnet_oc.

.. _methods-toa:

TOA simulations
---------------

In situ spectra (AERONET-OC reconstructed hyperspectrally by invRrs, HYPERNETS) or synthetic spectra are
propagated to the top of atmosphere with radcalnet_oc, with the atmosphere measured at the site
(AOT, Angström exponent, water vapour, ozone, NO\ :sub:`2`, pressure, wind;
:func:`~trishna_tbx.simulate.aeronet_oc_input`). The hyperspectral output (350-2500 nm) is convolved
with the responses of both sensors (:func:`~trishna_tbx.simulate.to_sensors`), which gives paired
TRISHNA and S2 band values for the band adjustment, the noise studies and the sensitivity to the
atmosphere (:func:`~trishna_tbx.simulate.sensitivity_grid`). See the
`radcalnet_oc methods <https://radcalnet-oc.readthedocs.io/en/latest/methods.html>`__.

.. _methods-inversion:

Information content of the inversion
------------------------------------

The semi-analytical model ``MultispectralSimple`` of invRrs (unknowns: Chl-a, :math:`a_{CDOM}(440)`,
SPM, offset) is integrated over the TRISHNA responses (:func:`~trishna_tbx.inversion.trishna_forward_model`).
The posterior covariance of the unknowns is linearized around a state :math:`x`, with a flat prior:

.. math::
   :label: infocontent

   \mathbf{S} = \left(\mathbf{J}^T \mathbf{S}_\epsilon^{-1} \mathbf{J}\right)^{-1},
   \qquad J_{kj} = \frac{\partial R_{rs,k}}{\partial x_j},

with :math:`\mathbf{S}_\epsilon` the diagonal covariance of the Rrs noise (5 % + 2·10\ :sup:`-4` sr\ :sup:`-1`
in the assessment). An unknown whose relative standard deviation exceeds 1, or whose correlation with
another is close to ±1, is not retrievable (:func:`~trishna_tbx.inversion.information_content`).

.. _methods-validation:

Validation with HYPERNETS
-------------------------

HYPERNETS gives hyperspectral :math:`R_{rs}` (400-900 nm) every 30 min, so the same measurement
validates both band sets once convolved :eq:`convolve`. The records of quality flag 0 and solar zenith
angle below 70° are kept, and matched with the Sentinel-2 acquisitions within 1 h. Each L1C subset is
cropped to 15 km around the site, averaged to 60 m and processed by grs with the TRISHNA-like and with
all the S2 bands (same version and settings); the median of the water pixels of a 3 × 3 window
(180 m) is compared with the in situ value. The archived standard GRS products (20 m) of the same
acquisitions are used as a reference of the processing.

Statistics (:func:`~trishna_tbx.validation.metrics`), with :math:`x` in situ and :math:`y` satellite:

.. math::
   :label: metrics

   \mathrm{bias} = \operatorname{median}(y - x),
   \qquad \mathrm{MdAPD} = 100\, \operatorname{median}\left|\frac{y - x}{x}\right|,
   \qquad \mathrm{RMSD} = \sqrt{\overline{(y - x)^2}},

plus the slope, intercept and correlation of the linear fit.
