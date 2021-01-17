.. _install:

Installation
============

**FastSpecFit** and all its dependencies are installed inside a Docker
container, making it easy to run either on a personal laptop (if you have the
necessary data) or at NERSC using *shifter*.

Development Mode
----------------

To run **FastSpecFit** in development mode requires a local check-out of the
​*fastspecfit* repository and a simple one-time setup procedure::

  git clone git@github.com:desihub/fastspecfit.git
  wget -O fastspecfit-setup.sh https://gist.github.com/moustakas/27bfdc6ea031d470f9600efa183aa5a1/raw

Once you have downloaded the *fastspecfit-setup.sh* script, edit it with the
path to your *fastspecfit* check-out as well as the *$DESI_ROOT*,
*$FASTSPECFIT_DATA*, and *$FASTSPECFIT_TEMPLATES* environment variables.

Next, simply load the shifter image, export the necessary environment variables,
and you are ready to go::

  ./fastspecfit-setup.sh shifter
  source fastspecfit-setup.sh env
  fastspecfit -h
    usage: fastspecfit [-h] [-n NTARGETS] [--firsttarget FIRSTTARGET] [--targetids TARGETIDS] [--mp MP] -o OUTFILE [--exposures] [--qa] [--photfit]
                   [--solve-vdisp]
                   [zbestfiles [zbestfiles ...]]


Production Mode
---------------

Write me.

Examples
--------

