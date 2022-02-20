.. _install:

Installation
============

``FastSpecFit`` and all its dependencies are installed inside a Docker
container, making it easy to run either on a personal laptop (if you have the
necessary data) or at NERSC using *shifter*.

Development Mode
----------------

To run **FastSpecFit** in development mode requires a local check-out of the
*fastspecfit* repository and a simple one-time setup procedure::

  cd /path/to/fastspecfit
  git clone https://github.com/desihub/fastspecfit.git

Then, edit the ``/path/to/fastspecfit/bin/fastspecfit-setup.sh`` file with the
path to your *fastspecfit* check-out as well as the *$DESI_ROOT*,
*$FASTSPECFIT_DATA*, and *$FASTSPECFIT_TEMPLATES* environment variables. (This
setup file is not used anywhere else, so you may choose to copy it somewhere
more convenient, like your home directory.)

Next, simply load the shifter image, export the necessary environment variables,
and you are ready to go::

  /path/to/fastspecfit/bin/fastspecfit-setup.sh shifter
  source /path/to/fastspecfit/bin/fastspecfit-setup.sh env
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

