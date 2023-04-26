.. _fuji vac:

Fuji VAC (EDR)
==============

.. contents:: Contents
    :depth: 3

Overview
--------

This page describes the ``Fuji`` value-added catalog, which will be publicly
released as part of the `DESI Early Data Release (DESI/EDR)`_ in
late-spring 2023.

Please refer to the :ref:`acknowledgments` section for the conditions for using
this VAC.

.. note::

   We only document the *v2.0* version of this VAC; the *v1.0* version is
   publicly accessible but has been deprecated and is not documented further.

Data Content & Access
---------------------

Data from the ``Fuji`` VAC can be accessed at any of the following links:

============================ ===================================================================
Data url                     https://data.desi.lbl.gov/public/edr/vac/fastspecfit/fuji/v2.0
Interactive web-app          https://fastspecfit.desi.lbl.gov
`NERSC`_ (for collaborators) ``/global/cfs/cdirs/desi/public/edr/vac/fastspecfit/fuji/v2.0``
============================ ===================================================================

For more information regarding the content and organization of the VAC, please
click on the following links:

* :ref:`Healpix Catalogs<healpix catalogs>`
* :ref:`Merged Catalogs<merged catalogs>`
* :ref:`Sample Selection<sample selection>`
* :ref:`Updated QSO Redshifts<qso redshifts>`

Summary Statistics
------------------
  
The following table summarizes the size and number of targets in each merged
catalog.

.. note::

   The last catalog listed in the table is a super-merge of all the preceding
   catalogs, i.e., a merge over all surveys and programs.

.. rst-class:: columns

=============================== ========= =================
File Name                       File Size Number of Targets
=============================== ========= =================
fastspec-fuji-cmx-other.fits    11.6 MB   2,771
fastspec-fuji-special-dark.fits 145 MB    35,647
fastspec-fuji-sv1-backup.fits   15.5 MB   3,683
fastspec-fuji-sv1-bright.fits   524 MB    126,677
fastspec-fuji-sv1-dark.fits     976 MB    235,881
fastspec-fuji-sv1-other.fits    141 MB    34,150
fastspec-fuji-sv2-backup.fits   602 KB    107
fastspec-fuji-sv2-bright.fits   193 MB    46,510
fastspec-fuji-sv2-dark.fits     218 MB    52,771
fastspec-fuji-sv3-backup.fits   6.62 MB   1,564
fastspec-fuji-sv3-bright.fits   1.08 GB   265,324
fastspec-fuji-sv3-dark.fits     2.41 GB   592,394
fastspec-fuji.fits              5.7 GB    1,397,479
=============================== ========= =================

The following table summarizes the number of QSO targets whose redshift has been
updated using the procedure documented :ref:`here<qso redshifts>`.

.. rst-class:: columns

=============================== ================= ===============================
Catalog                         Number of Targets Number with Corrected Redshifts
=============================== ================= ===============================
fastspec-fuji-cmx-other.fits    2,771             56
fastspec-fuji-special-dark.fits 35,647            311
fastspec-fuji-sv1-backup.fits   3,683             100
fastspec-fuji-sv1-bright.fits   126,677           64
fastspec-fuji-sv1-dark.fits     235,881           3,749
fastspec-fuji-sv1-other.fits    34,150            168
fastspec-fuji-sv2-backup.fits   107               0
fastspec-fuji-sv2-bright.fits   46,510            8
fastspec-fuji-sv2-dark.fits     52,771            1,019
fastspec-fuji-sv3-backup.fits   1,564             0
fastspec-fuji-sv3-bright.fits   265,324           132
fastspec-fuji-sv3-dark.fits     592,394           3,397
fastspec-fuji.fits              1,397,479         9,004
=============================== ================= ===============================

Known Issues
------------

This section documents any issues or problems which were identified with the VAC
after its final release.

* The tied amplitudes of [OII]7320,7330 doublet were reversed in the line fitting [`#119`_].

.. _`#119`: https://github.com/desihub/fastspecfit/issues/119

To report projects or to request new features please `open a ticket`_.

.. _`DESI Early Data Release (DESI/EDR)`: https://data.desi.lbl.gov/public/edr
.. _`NERSC`: https://nersc.gov
.. _`open a ticket`: https://github.com/desihub/fastspecfit/issues
