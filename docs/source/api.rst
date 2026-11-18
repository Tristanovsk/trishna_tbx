.. _api:

API reference
=============

The sensors are importable from the package itself, the other tools from their modules:

.. code-block:: python

   from trishna_tbx import TRISHNA, make_s2
   from trishna_tbx import algorithms, emulate, owt, water_vapor, validation

   print(algorithms.coverage([make_s2(), TRISHNA]))

The modules that rely on a sibling package import it inside the functions: ``emulate`` (GRSdriver),
``owt`` (GRSl2bgen), ``simulate`` and ``water_vapor`` (radcalnet_oc), ``inversion`` (invRrs),
``validation`` (grs).

Sensors
-------

.. autosummary::
   :toctree: generated
   :template: autosummary/class.rst
   :nosignatures:

   trishna_tbx.sensor.Band
   trishna_tbx.sensor.Sensor

.. autosummary::
   :toctree: generated
   :nosignatures:

   trishna_tbx.sensor.make_trishna
   trishna_tbx.sensor.make_s2
   trishna_tbx.sensor.make_s2_trishna_like
   trishna_tbx.sensor.super_gaussian
   trishna_tbx.sensor.solar_irradiance
   trishna_tbx.sensor.TRISHNA
   trishna_tbx.sensor.TRISHNA_BANDS
   trishna_tbx.sensor.S2_BANDS
   trishna_tbx.sensor.S2_TO_TRISHNA
   trishna_tbx.sensor.WL

Emulation of TRISHNA from Sentinel-2
------------------------------------

.. autosummary::
   :toctree: generated
   :template: autosummary/class.rst
   :nosignatures:

   trishna_tbx.emulate.SBAF

.. autosummary::
   :toctree: generated
   :nosignatures:

   trishna_tbx.emulate.grs_band_set
   trishna_tbx.emulate.open_s2_l1c_as_trishna
   trishna_tbx.emulate.l1c_dataset_as_trishna
   trishna_tbx.emulate.restrict_bands
   trishna_tbx.emulate.owt_band_adjustment
   trishna_tbx.emulate.add_noise
   trishna_tbx.emulate.S2_BAND_ORDER
   trishna_tbx.emulate.S2_NOMINAL_WL

Optical water types
-------------------

.. autosummary::
   :toctree: generated
   :nosignatures:

   trishna_tbx.owt.load_reference
   trishna_tbx.owt.convolve_reference
   trishna_tbx.owt.sam
   trishna_tbx.owt.classify
   trishna_tbx.owt.separability
   trishna_tbx.owt.merge_classes
   trishna_tbx.owt.sample_members
   trishna_tbx.owt.benchmark
   trishna_tbx.owt.TRISHNA_MERGE
   trishna_tbx.owt.BI2024_LABELS

Water-quality algorithms
------------------------

.. autosummary::
   :toctree: generated
   :template: autosummary/class.rst
   :nosignatures:

   trishna_tbx.algorithms.Algo

.. autosummary::
   :toctree: generated
   :nosignatures:

   trishna_tbx.algorithms.coverage
   trishna_tbx.algorithms.run_all
   trishna_tbx.algorithms.retune_ocx
   trishna_tbx.algorithms.REGISTRY
   trishna_tbx.algorithms.ocx
   trishna_tbx.algorithms.oc2
   trishna_tbx.algorithms.kd490_kd2
   trishna_tbx.algorithms.kdpar_rd22
   trishna_tbx.algorithms.spm_nechad10
   trishna_tbx.algorithms.turb_dogliotti15
   trishna_tbx.algorithms.ndci
   trishna_tbx.algorithms.green_red_index
   trishna_tbx.algorithms.red_nir_ratio

Water vapour
------------

.. autosummary::
   :toctree: generated
   :template: autosummary/class.rst
   :nosignatures:

   trishna_tbx.water_vapor.RatioRetrieval

.. autosummary::
   :toctree: generated
   :nosignatures:

   trishna_tbx.water_vapor.band_transmittance
   trishna_tbx.water_vapor.h2o_optical_thickness
   trishna_tbx.water_vapor.compare_trishna_s2
   trishna_tbx.water_vapor.noise_error
   trishna_tbx.water_vapor.PAIRS

TOA simulations (radcalnet_oc)
------------------------------

.. autosummary::
   :toctree: generated
   :nosignatures:

   trishna_tbx.simulate.build_input
   trishna_tbx.simulate.aeronet_oc_input
   trishna_tbx.simulate.simulate_toa
   trishna_tbx.simulate.to_sensors
   trishna_tbx.simulate.sensitivity_grid

Inversion (invRrs)
------------------

.. autosummary::
   :toctree: generated
   :template: autosummary/class.rst
   :nosignatures:

   trishna_tbx.inversion.SensorRSR

.. autosummary::
   :toctree: generated
   :nosignatures:

   trishna_tbx.inversion.trishna_forward_model
   trishna_tbx.inversion.information_content

Validation with HYPERNETS
-------------------------

.. autosummary::
   :toctree: generated
   :nosignatures:

   trishna_tbx.validation.load_hypernets
   trishna_tbx.validation.convolve
   trishna_tbx.validation.find_matchups
   trishna_tbx.validation.process_matchup
   trishna_tbx.validation.find_l2a
   trishna_tbx.validation.extract_l2a_window
   trishna_tbx.validation.metrics

Readers
-------

.. autosummary::
   :toctree: generated
   :nosignatures:

   trishna_tbx.io.open_grs_l2a
   trishna_tbx.io.parse_window
   trishna_tbx.io.water_mask
   trishna_tbx.io.to_60m
   trishna_tbx.io.list_l2a
   trishna_tbx.io.cams_value
