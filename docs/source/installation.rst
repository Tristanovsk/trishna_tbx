Installation
============

Python environment
------------------

Python >= 3.10 is required (the package is developed with Python 3.12). With conda:

.. code-block:: bash

   conda create -n trishna_tbx python=3.12
   conda activate trishna_tbx

On the CNES JupyterHub, the ``sandbox_env`` kernel already has the sibling packages: only
``trishna_tbx`` has to be added.

Install the package
-------------------

From the root of a clone of the repository:

.. code-block:: bash

   pip install -e .

The core of the package (sensor models, band adjustment, OWT classification, algorithms) only
needs ``numpy``, ``xarray``, ``pandas``, ``scipy`` and ``openpyxl``.

Sibling packages
----------------

The processing chains come from the sibling packages, imported only by the functions that need
them:

.. list-table::
   :header-rows: 1
   :widths: 20 40 40

   * - Package
     - Used for
     - Used by
   * - `GRSdriver <https://github.com/Tristanovsk/GRSdriver>`__
     - reading of the Sentinel-2 L1C, ESA spectral responses of MSI
     - :func:`~trishna_tbx.sensor.make_s2` (Gaussian responses without it),
       :func:`~trishna_tbx.emulate.open_s2_l1c_as_trishna`
   * - `grs <https://grs.readthedocs.io>`__ (GRSprocessor)
     - atmospheric and sunglint correction
     - :func:`~trishna_tbx.validation.process_matchup`, notebooks of section 3;
       :func:`~trishna_tbx.water_vapor.h2o_optical_thickness` (without radcalnet_oc)
   * - `GRSl2bgen <https://github.com/Tristanovsk/GRSl2bgen>`__
     - reference spectra of the optical water types
     - :mod:`trishna_tbx.owt`
   * - `radcalnet_oc <https://radcalnet-oc.readthedocs.io>`__
     - TOA simulations, solar spectrum, water vapour absorption
     - :mod:`trishna_tbx.simulate`, :func:`~trishna_tbx.sensor.solar_irradiance` (black body
       without it), :mod:`trishna_tbx.water_vapor`
   * - `invRrs <https://github.com/Tristanovsk/invRrs>`__
     - semi-analytical inversion of Rrs
     - :mod:`trishna_tbx.inversion`
   * - `hgrs <https://github.com/Tristanovsk/hgrs>`__
     - coupled water vapour model (to be ported for water pixels)
     - --

Install them from their repositories (editable mode), or all at once with the ``full`` extra if
they are available in your package index:

.. code-block:: bash

   pip install -e ".[full]"

grs also needs its look-up tables and the CAMS files, see the
`grs documentation <https://grs.readthedocs.io/en/latest/installation.html>`__.

Run the tests
-------------

.. code-block:: bash

   pip install -e ".[dev]"
   pytest

The tests that need a sibling package (radcalnet_oc, GRSl2bgen) are skipped when it is not
installed.

Build this documentation
------------------------

.. code-block:: bash

   pip install -e ".[docs]"
   sphinx-build -b html docs/source docs/build/html

The tutorials are the notebooks of ``notebooks/``, copied to ``docs/source/tutorials`` at each
build. They are not executed: they are rendered with the outputs saved in them. To update them,
re-execute them with ``notebooks/build_notebooks.py --execute`` (see :doc:`usage`).

The sibling packages are mocked for the API reference: they are not needed to build the
documentation.
