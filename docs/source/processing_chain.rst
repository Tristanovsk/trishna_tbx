Processing Chain Description
=============================

This page documents trishna_tbx following the
`Processing Chain Documentation Template <https://processing-chain-guidelines.readthedocs.io/en/latest/Data_processing_chain_template/>`_.

Process description
--------------------

Description
~~~~~~~~~~~

trishna_tbx is a study toolbox, not an operational processor: it emulates TRISHNA VSWIR data from
Sentinel-2 and simulates them from in situ measurements, runs the existing water processing chains
(grs, GRSl2bgen algorithms, invRrs) on the TRISHNA and on the Sentinel-2 band sets, and compares the
results with each other and with in situ data. The equations are described in :doc:`methods`.

It is used through its Python API (:doc:`api`), two scripts (``scripts/prelaunch_assessment.py``,
``scripts/l2a_s2_vs_trishna.py``) and the VRE notebooks (:doc:`tutorials`).

Application Domain (Granule)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Depending on the function:

- **image**: one Sentinel-2 L1C (SAFE or NetCDF subset) or one GRS L2A product, at 60 m;
- **time series of images**: the GRS L2A products of one tile (datacube);
- **site**: the in situ time series of one AERONET-OC or HYPERNETS site and the matching Sentinel-2
  subsets;
- **spectral library**: reference or synthetic spectra, without geolocation (OWT benchmark, band
  adjustment, information content).

Scheduling and Triggers
~~~~~~~~~~~~~~~~~~~~~~~~

No scheduler: the functions are run interactively (notebooks) or in scripts. The GRS runs of the
validation are cached, so that the notebooks can be re-executed quickly.

Dataflow
--------

.. only:: html

   .. mermaid:: _diagrams/dataflow.mmd

.. only:: not html

   See the dataflow diagram in the
   `online documentation <https://trishna-tbx.readthedocs.io/en/latest/processing_chain.html#dataflow>`__.

Inputs
------

.. list-table::
   :header-rows: 1
   :widths: 25 15 60

   * - Data Type Name
     - Cardinality
     - Selection Criteria
   * - Sentinel-2 L1C
     - 0..n
     - ``.SAFE`` folder (read by GRSdriver) or NetCDF subset of GRSdriver; for the emulation of
       TRISHNA L1C and the atmospheric correction
   * - GRS L2A
     - 0..n
     - NetCDF or Zarr, stacked ``Rrs`` or ``Rrs_<wl>`` (beam profile), ``mask`` or ``ndwi``; for the
       Level-2 emulation, products, OWT and datacube
   * - CAMS
     - 0..n
     - global forecast files of grs (``cams_forecast_YYYY-MM.nc``), for grs
   * - HYPERNETS L2A
     - 0..1 per site
     - stacked NetCDF time series (``Rrs`` on ``time``, ``wl``, ``quality_flag``, ``solar_zenith_angle``)
   * - AERONET-OC
     - 0..1 per site
     - Level 2.0 of ``radcalnet_oc.aeronet_oc``, optionally with the hyperspectral invRrs
       reconstructions
   * - Auxiliary data of the sibling packages
     - 1..1
     - S2 responses (GRSdriver), OWT references (GRSl2bgen), solar spectrum and water vapour
       absorption (radcalnet_oc), look-up tables of grs and radcalnet_oc

Outputs
-------

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   * - Data Type Name
     - Cardinality
     - Description
   * - TRISHNA-like L1C
     - 0..n
     - :class:`xarray.Dataset` ready for grs (bands, angles and masks at 60 m, ``wl_to_process``)
   * - TRISHNA-like L2A and products
     - 0..n
     - :class:`xarray.Dataset` of Rrs, products (:func:`~trishna_tbx.algorithms.run_all`) and OWT;
       NetCDF with ``--out`` of ``l2a_s2_vs_trishna.py``
   * - Matchups
     - 0..n
     - cached NetCDF per acquisition (``<l1c name>_win3_crop15km.nc``): Rrs statistics of the window,
       number of water pixels, AOT, per band set
   * - Assessment tables
     - 1..n
     - :class:`pandas.DataFrame` (coverage, OWT agreement, SBAF, validation statistics) printed by the
       scripts and notebooks

Return Codes
------------

- Scripts: ``0`` when done; ``2`` for invalid command-line arguments (argparse); ``1`` for an error
  during the processing (Python traceback).
- API: the errors are raised as Python exceptions, e.g. :class:`ValueError` from
  :meth:`Algo.run <trishna_tbx.algorithms.Algo.run>` when a sensor lacks a band of the algorithm.
  :func:`~trishna_tbx.validation.process_matchup` catches the failures of grs (no water, cloudy subset),
  logs them and returns None.

Log Format
----------

trishna_tbx uses the Python :py:mod:`logging` module (warnings of
:func:`~trishna_tbx.validation.process_matchup`) and follows the logging configuration of the calling
application. :func:`~trishna_tbx.validation.process_matchup` sets the level of the root logger to
``WARNING`` to silence the progress messages of grs.

Required Resources
-------------------

.. list-table::
   :header-rows: 1
   :widths: 30 30 40

   * - Resource Type
     - Quantity
     - Notes
   * - CPU
     - 1 core for most tutorials
     - grs uses its own parallelism
   * - RAM
     - 4 GB; 8-16 GB for the GRS runs on a full tile
     - the validation crops the subsets to 15 km
   * - Execution time
     - seconds to minutes per tutorial; 1.5-2 h for the first GRS runs of the 64 Berre matchups
     - with the results cached, the validation notebook runs in about 3 minutes
   * - Disk — output
     - a few MB per matchup, a few hundred MB for a datacube
     - in ``<section>/output/`` of the notebooks

Data Types
----------

.. list-table::
   :header-rows: 1
   :widths: 18 32 15 35

   * - Data Name
     - Description
     - Granule
     - Nomenclature
   * - S2 L1C
     - Sentinel-2 Level-1C
     - tile or subset
     - ESA naming, ``S2B_MSIL1C_<YYYYMMDDTHHMMSS>_N<baseline>_R<orbit>_T<tile>_<...>``
   * - GRS L2A
     - water reflectance of GRS
     - tile or subset
     - ``S2B_MSIL2Agrs_<...>.nc``, or archive ``<root>/YYYY/MM/DD/<name>/<name>.nc``
   * - HYPERNETS
     - hyperspectral Rrs time series
     - one site
     - ``<SITE>_hypernets.nc`` (e.g. ``BEFR``)
   * - Matchup cache
     - window statistics, both band sets
     - one acquisition
     - ``<L1C name>_win<n>_crop<km>km.nc``
