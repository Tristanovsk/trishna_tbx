Usage
=====

``trishna_tbx`` is a Python toolbox: its functions are used from scripts and notebooks. Two scripts
reproduce the main numbers of the assessment, and the VRE notebooks (:doc:`tutorials`) show each
step on real data.

Sensors
-------

.. code-block:: python

   from trishna_tbx import TRISHNA, make_s2

   s2 = make_s2('S2A')                       # ESA spectral responses (GRSdriver)
   TRISHNA.srf                               # xarray.DataArray (band, wl), super-Gaussian responses
   TRISHNA.noise_reflectance(sza=30.)        # 1-sigma TOA reflectance noise of each band
   s2.noise_reflectance(30., n_aggregated=9) # S2 20 m bands averaged to 60 m

   rrs_bands = TRISHNA.convolve(rrs_hyper, 'wl', weight=solar_irradiance())

Measured TRISHNA spectral responses will replace the modelled ones: ``TRISHNA.srf`` is an ordinary
:class:`xarray.DataArray` that can be replaced, or a new sensor can be built with
:class:`~trishna_tbx.sensor.Sensor`.

TRISHNA-like images from Sentinel-2
-----------------------------------

Level 1 (for the atmospheric correction with grs):

.. code-block:: python

   import grs
   from trishna_tbx.emulate import open_s2_l1c_as_trishna

   prod = open_s2_l1c_as_trishna('S2B_MSIL1C_<...>.SAFE', resolution=60)
   process = grs.Process()
   process.execute(prod, odir='L2A', cams_file='cams_forecast_2025-05.nc')

For L1C subsets already saved as NetCDF (e.g. around in situ sites), use
:func:`~trishna_tbx.emulate.l1c_dataset_as_trishna`.

Level 2 (from an existing GRS L2A product):

.. code-block:: python

   from trishna_tbx.io import open_grs_l2a, water_mask, to_60m
   from trishna_tbx.emulate import owt_band_adjustment
   from trishna_tbx.sensor import make_s2

   ds = open_grs_l2a('S2B_MSIL2Agrs_<...>.nc', window='y=0:3600,x=0:960')
   rrs = to_60m(ds.Rrs.where(water_mask(ds)), native=20.)
   s2_tri = make_s2().subset(['B2', 'B3', 'B4', 'B8A'])
   rrs_tri, owt, factors = owt_band_adjustment(rrs.sel(wl=[490, 560, 665, 865]), s2_tri)

Products and optical water types
--------------------------------

.. code-block:: python

   from trishna_tbx import algorithms, owt

   algorithms.coverage([s2, TRISHNA])          # which algorithm each sensor can run
   products = algorithms.run_all(rrs_tri, TRISHNA)

   mean, sd = owt.load_reference('Bi2024')
   ref = owt.convolve_reference(mean, TRISHNA.subset(['Blue', 'Green', 'Red', 'NIR']))
   cls, angle = owt.classify(rrs_bands, ref, dim='band')
   cls7 = owt.merge_classes(cls)               # 7-class TRISHNA set

:meth:`Algo.run <trishna_tbx.algorithms.Algo.run>` refuses to run an algorithm whose bands the sensor
does not have (``strict=True``), instead of silently using the nearest band.

Water vapour, simulations and inversion
---------------------------------------

.. code-block:: python

   from trishna_tbx.water_vapor import RatioRetrieval, compare_trishna_s2, noise_error
   from trishna_tbx.simulate import aeronet_oc_input, simulate_toa, to_sensors
   from trishna_tbx.inversion import trishna_forward_model, information_content

   w = RatioRetrieval('TRISHNA').retrieve(r860, r910, sza, vza)   # kg m-2

   inp = aeronet_oc_input('Venise', '2021-06', hyp_dir='/DATA/AERONET/OCv3_invRrs/data')
   bands = to_sensors(simulate_toa(inp))                         # {'TRISHNA': ..., 'Sentinel-2 A': ...}

Validation with HYPERNETS
-------------------------

.. code-block:: python

   from trishna_tbx import validation as val

   hyp = val.load_hypernets('BEFR_hypernets.nc')
   mu = val.find_matchups(hyp.time.values, 'L1C/BEFR', max_dt='1h')
   res = [val.process_matchup(f, 'cams', 'cache', site_xy=(x0, y0)) for f in mu.l1c]
   val.metrics(insitu, satellite)

The GRS runs are cached: the first run of the 64 Berre matchups takes about 1.5 to 2 hours, the
next runs a few minutes.

Scripts
-------

``scripts/prelaunch_assessment.py`` reproduces the numbers of :doc:`assessment` that do not need an
image (algorithm coverage, OWT benchmark, band adjustment, water vapour; ``--inversion`` adds the
information content of the inversion, which needs invRrs):

.. code-block:: bash

   python scripts/prelaunch_assessment.py [--inversion]

``scripts/l2a_s2_vs_trishna.py`` compares the products and the OWT of a GRS L2A image with all
the S2 bands and with the TRISHNA-like bands, at 60 m:

.. code-block:: text

   usage: l2a_s2_vs_trishna.py [-h] [--native NATIVE] [--crop CROP] [--out OUT] l2a

   l2a              GRS L2A NetCDF file
   --native NATIVE  native pixel size (m), default 20
   --crop CROP      pixel window 'y0:y1,x0:x1'
   --out OUT        NetCDF file of the 60 m products and OWT

Notebooks
---------

The VRE notebooks are written by ``notebooks/build_notebooks.py`` (edit this file, not the
``.ipynb``):

.. code-block:: bash

   cd notebooks
   python build_notebooks.py                 # write all notebooks
   python build_notebooks.py --execute       # write and execute (keeps the outputs)
   python build_notebooks.py --execute owt   # only the notebooks whose path contains "owt"
   python build_notebooks.py --execute --allow-errors   # keep errors in the notebook (debugging)
   python check_nb.py                        # executed cells / first error of each notebook

The data paths are set in the *Data* cell of each notebook, or with environment variables:

.. list-table::
   :header-rows: 1
   :widths: 30 45 25

   * - Variable
     - Default data
     - Notebooks
   * - ``TRISHNA_TBX_S2_L1C``
     - S2B L1C 2025-05-26, tile T31TFJ (Berre lagoon, Gulf of Fos)
     - emulation, GRS, products
   * - ``TRISHNA_TBX_CAMS``
     - CAMS global forecast, 2025-05
     - emulation, GRS
   * - ``TRISHNA_TBX_L2A_EXAMPLE``
     - GRS L2A (Zarr pyramid, v3.1.1), tile T30TWT (Loire estuary), 2017-04-12, local S3 bucket
       ``grstbx-s3``
     - products
   * - ``TRISHNA_TBX_L2A_DIR``
     - GRS L2A, tile T33PVQ (Lake Chad), January-February 2022
     - OWT, datacube
   * - ``TRISHNA_TBX_AERONET_HYP``
     - invRrs hyperspectral AERONET-OC spectra
     - TOA simulations
   * - ``TRISHNA_TBX_HYPERNETS``
     - HYPERNETS BEFR (Berre lagoon) L2A, 2021-2023
     - validation
   * - ``TRISHNA_TBX_L1C_SUBSETS``
     - S2 L1C subsets (NetCDF) around BEFR
     - validation
   * - ``TRISHNA_TBX_CAMS_DIR``
     - monthly CAMS files ``cams_forecast_YYYY-MM.nc``
     - validation
   * - ``TRISHNA_TBX_OUTPUT``
     - ``./output`` (GRS results, datacube)
     - emulation, GRS, datacube
