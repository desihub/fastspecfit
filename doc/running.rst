.. _running_fastspecfit:

Running FastSpecFit
===================

.. note::
  Before running any of the examples below, we assume that you have successfully
  :ref:`installed and set up<install>` ``FastSpecFit``.

.. _nersc installation:

Overview
--------

Running ``FastSpecFit`` is accomplished through a handful of high-level Python
scripts. The two primary, independent scripts which can be run on one (or a
small number) of `Redrock`_ redshift catalogs are:

  * :ref:`fastspec description`, to model DESI spectrophotometry; and
  * :ref:`fastphot description`, to model DESI broadband photometry.

Note that both scripts require (and assume) DESI-derived redshifts and the DESI
data model as inputs.

In addition, there are two key support routines, which we describe in more
detail below:

  * ``fastspecfit-qa``, to build quality assurance (QA) figures; and
  * ``mpi-fastspecfit``, to execute a variety of tasks (in parallel) on larger
    numbers of input files or catalogs.

.. _`fastspec description`:

fastspec
--------

To model the spectrum of a single object, we just need to provide ``fastspec``
an input Redrock catalog and an (arbitrary) output filename::

  $> fastspec $DESI_ROOT/spectro/redux/guadalupe/healpix/main/bright/300/30022/redrock-main-bright-30022.fits \
    --ntargets 1 --outfile fastspec-example.fits

See the :ref:`fastspec data model<fastspec>` for a full description of the
contents of the file which is written out. We can visualize the results by
generating a figure::

  $> fastspecfit-qa fastspec-example.fits







The arguments to ``fastspec`` can be inspected by invoking the script with the
`--help` option::

  $> fastspec --help
  
  usage: fastspec [-h] -o OUTFILE [--mp MP] [-n NTARGETS] [--firsttarget FIRSTTARGET] \
    [--targetids TARGETIDS] [--solve-vdisp] [--verbose] [redrockfiles ...]
  
  positional arguments:
    redrockfiles          Full path to input redrock file(s). (default: None)
  
  options:
    -h, --help            show this help message and exit
    -o OUTFILE, --outfile OUTFILE
                          Full path to output filename (required). (default: None)
    --mp MP               Number of multiprocessing threads per MPI rank. (default: 1)
    -n NTARGETS, --ntargets NTARGETS
                          Number of targets to process in each file. (default: None)
    --firsttarget FIRSTTARGET
                          Index of first object to to process in each file, zero-indexed. (default: 0)
    --targetids TARGETIDS
                          Comma-separated list of TARGETIDs to process. (default: None)
    --solve-vdisp         Solve for the velocity dispersion (only when using fastspec). (default: False)
    --verbose             Be verbose (for debugging purposes). (default: False)

For example, do::

  fastspec ${DESI_ROOT}/spectro/redux/fuji/tiles/cumulative/80613/20210324/redrock-4-80613-thru20210324.fits -o fastspec.fits --targetids 39633345008634465


Write me.

.. _`fastphot description`:

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
