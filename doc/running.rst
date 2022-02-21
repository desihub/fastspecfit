.. _running_fastspecfit:

Running FastSpecFit
===================

.. note::
  Before running any of the examples below, we assume that you have successfully
  installed and set up ``FastSpecFit``. If not, please see the :ref:`install`
  instructions. 

Overview
--------

Running ``FastSpecFit`` is accomplished through a handful of high-level Python
scripts. The two primary, independent scripts which can be run on one (or a
small number) of `Redrock`_ redshift catalogs are:

  * :ref:`fastspec`, to model DESI spectrophotometry; and
  * :ref:`fastphot`, to model DESI broadband photometry.

Note that both scripts require (and assume) DESI-derived redshifts and the DESI
data model as inputs.

In addition, there are two key support routines, which we describe in more
detail below:

  * ``fastspecfit-qa``, to build quality assurance (QA) figures; and
  * ``mpi-fastspecfit``, to execute a variety of tasks (in parallel) on larger
    numbers of input files or catalogs.

.. _`fastspec`:

fastspec
--------

The arguments to ``fastspec`` can be inspected by invoking the script with the
`--help` option::

  fastspec --help
  usage: fastspec [-h] [-n NTARGETS] [--firsttarget FIRSTTARGET] [--targetids TARGETIDS] \
    [--mp MP] -o OUTFILE [--solve-vdisp] [--verbose] [--specprod SPECPROD] [redrockfiles ...]
  
  positional arguments:
    redrockfiles          Full path to input redrock file(s). (default: None)
  
  optional arguments:
    -h, --help            show this help message and exit
    -n NTARGETS, --ntargets NTARGETS
                          Number of targets to process in each file. (default: None)
    --firsttarget FIRSTTARGET
                          Index of first object to to process in each file (0-indexed). (default: 0)
    --targetids TARGETIDS
                          Comma-separated list of target IDs to process. (default: None)
    --mp MP               Number of multiprocessing processes per MPI rank or node. (default: 1)
    -o OUTFILE, --outfile OUTFILE
                          Full path to output filename. (default: None)
    --solve-vdisp         Solve for the velocity dispersion (only when using fastspec). (default: False)
    --verbose             Be verbose (for debugging purposes). (default: False)
    --specprod SPECPROD   Spectroscopic production. (default: everest)

Write me.

.. _`fastphot`:

fastphot
--------

Write me.

Other
-----

Next::

  fastphot -h
    usage: fastphot [-h] [-n NTARGETS] [--firsttarget FIRSTTARGET] [--targetids TARGETIDS] [--mp MP] -o OUTFILE [--exposures] [--solve-vdisp] [zbestfiles [zbestfiles ...]]
    
    positional arguments:
      zbestfiles            Full path to input zbest file(s). (default: None)
    
    optional arguments:
      -h, --help            show this help message and exit
      -n NTARGETS, --ntargets NTARGETS
                            Number of targets to process in each file. (default: None)
      --firsttarget FIRSTTARGET
                            Index of first object to to process in each file (0-indexed). (default: 0)
      --targetids TARGETIDS
                            Comma-separated list of target IDs to process. (default: None)
      --mp MP               Number of multiprocessing processes per MPI rank or node. (default: 1)
      -o OUTFILE, --outfile OUTFILE
                            Full path to output filename. (default: None)
      --exposures           Fit the individual exposures (not the coadds). (default: False)
      --solve-vdisp         Solve for the velocity disperion. (default: False)

Production Mode
---------------

Write me.

Examples
--------

.. _`RedRock`: https://github.com/desihub/redrock
