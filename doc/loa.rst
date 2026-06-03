.. _loa vac:
.. _loa fastspec vac:

Loa fastspec VAC (DR2)
======================

.. contents:: Contents
    :depth: 3

Overview
--------

This page describes the ``Loa`` fastspec value-added catalog, which contains
full spectrophotometric fitting results (continuum + emission lines) and will
be publicly released as part of `DESI Data Release 2 (DESI/DR2)`_. For the
companion photometry-only catalog, see the :ref:`Loa fastphot VAC<loa fastphot vac>`.

Please refer to the :ref:`acknowledgments` section for the conditions for using
this VAC.

.. note::

   We only document the *latest* version of this VAC here; please see below for
   information regarding :ref:`previous versions<previous versions - loa>`.

Data Content & Access
---------------------

Data from the ``Loa`` fastspec VAC can be accessed at any of the following links:

============================ ==================================================================
Data url                     TBD
`NERSC`_ (for collaborators) TBD
============================ ==================================================================

For more information regarding the content and organization of the VAC, please
click on the following links:

* :ref:`Healpix Catalogs<healpix catalogs>`
* :ref:`Merged Catalogs<merged catalogs>`
* :ref:`Sample Selection<sample selection>`
* :ref:`Updated QSO Redshifts<qso redshifts>`

Summary Statistics
------------------

The following tables summarize the size and number of targets in each merged
catalog, organized by survey and program.

TBD

The following tables summarize the number of QSO targets whose redshift has been
updated using the procedure documented :ref:`here<qso redshifts>`.

TBD

Code & Template Versions
------------------------

The following tables document the code versions and environment variables used
to produce this VAC. For details regarding the revision history of
``FastSpecFit``, please see the `change log`_.

Note that the tagged dependencies can be retrieved from any FITS file with the
following bit of code::

  import fitsio
  from desiutil.depend import Dependencies
  codever = Dependencies(fitsio.read_header('/path/to/fastspecfit/file.fits', ext=0))
  for codename, version in codever.items():
      print(codename, version)

.. rst-class:: columns

================ ==========
Software Package Version(s)
================ ==========
python           TBD
numpy            TBD
scipy            TBD
astropy          TBD
yaml             TBD
matplotlib       TBD
fitsio           TBD
mpi4py           TBD
healpy           TBD
desiutil         TBD
desispec         TBD
desitarget       TBD
speclite         TBD
fastspecfit      TBD
================ ==========

.. rst-class:: columns

==================== ============================================================
Environment Variable Value
==================== ============================================================
DESI_ROOT            /global/cfs/cdirs/desi
DUST_DIR             TBD
FPHOTO_DIR           TBD
FTEMPLATES_DIR       /dvs_ro/cfs/cdirs/desi/public/external/templates/fastspecfit
FTEMPLATES_FILE      TBD (see `README.txt`_)
==================== ============================================================

.. _previous versions - loa:

Notes & Known Issues
--------------------

v1.0 (latest release)
~~~~~~~~~~~~~~~~~~~~~~

* Release date: TBD
* ``FastSpecFit`` version: TBD
* Templates: TBD (see `README.txt`_).
* Notes:

  * Initial release of the ``Loa`` fastspec VAC.

* Known issues:

  * None at this time.

.. _`DESI Data Release 2 (DESI/DR2)`: https://data.desi.lbl.gov/doc/releases/dr2
.. _`NERSC`: https://nersc.gov
.. _`open a ticket`: https://github.com/desihub/fastspecfit/issues
.. _`change log`: https://github.com/desihub/fastspecfit/blob/main/doc/changes.rst
.. _`README.txt`: https://data.desi.lbl.gov/public/external/templates/fastspecfit/README.txt
