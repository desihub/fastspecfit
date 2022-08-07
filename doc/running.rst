.. _running_fastspecfit:

Running FastSpecFit
===================

.. note::
   
   Before running any of the examples below, we assume that you have
   successfully :ref:`installed and set up<install>` ``FastSpecFit``.

.. contents:: Contents
    :depth: 3

Overview
--------

Running ``FastSpecFit`` is accomplished through a handful of high-level Python
scripts. The two primary, independent scripts which can be run on one (or a
small number) of `Redrock`_ redshift catalogs are:

  * :ref:`fastspec example`, to model DESI spectrophotometry; and
  * :ref:`fastphot example`, to model DESI broadband photometry.

Note that both scripts require (and assume) DESI spectra and redshits as
inputs. In addition, there are two key support routines, which we describe in
more detail below:

  * ``fastspecfit-qa``, to build quality assurance (QA) figures; and
  * ``mpi-fastspecfit``, to execute a variety of tasks (in parallel) on larger
    numbers of input files or catalogs.

.. _`RedRock`: https://github.com/desihub/redrock

.. _`fastspec example`:

One fastspec Example
--------------------

To model the spectrum of a single object, we just need to provide ``fastspec``
an input Redrock catalog and an (arbitrary) output filename::


  $> fastspec /global/cfs/cdirs/desi/spectro/redux/fuji/healpix/sv1/bright/71/7108/redrock-sv1-bright-7108.fits \
    --targetids 39633345008634465 --outfile fastspec-example.fits
    
  INFO:continuum.py:137:__init__: Reading /global/cfs/cdirs/desi/science/gqp/templates/SSP-CKC14z/v1.0/SSP_Padova_CKC14z_Kroupa_Z0.0190.fits
  INFO:continuum.py:137:__init__: Reading /global/cfs/cdirs/desi/science/gqp/templates/SSP-CKC14z/v1.0/SSP_Padova_CKC14z_Kroupa_Z0.0190.fits
  INFO:fastspecfit.py:165:fastspec: Initializing the classes took: 1.19 sec
  INFO:io.py:296:select: Reading and parsing 1 unique redrockfile(s)
  INFO:io.py:348:select: specprod=fuji, coadd_type=healpix, survey=sv1, program=bright, healpix=7108
  INFO:io.py:443:select: Updating QSO redshifts using a QN threshold of 0.95.
  INFO:io.py:569:select: Gathered photometric metadata for 1 objects in 0.32 sec
  INFO:io.py:658:read_and_unpack: Reading 1 spectra from /global/cfs/cdirs/desi/spectro/redux/fuji/healpix/sv1/bright/71/7108/coadd-sv1-bright-7108.fits
  INFO:spectra.py:291:read_spectra: iotime 0.355 sec to read coadd-sv1-bright-7108.fits at 2022-08-07T04:00:52.389206
  INFO:continuum.py:836:get_linesigma: Forbidden masking sigma=98.374 km/s and S/N=73.893
  INFO:continuum.py:836:get_linesigma: Balmer masking sigma=88.367 km/s and S/N=84.792
  INFO:continuum.py:836:get_linesigma: UV/Broad masking sigma=0.000 km/s and S/N=0.000
  INFO:io.py:855:read_and_unpack: Read data for 1 objects in 1.24 sec
  INFO:fastspecfit.py:181:fastspec: Reading and unpacking the 1 spectra to be fitted took: 3.43 sec
  INFO:continuum.py:1992:continuum_specfit: Preparing the models took 0.55 sec
  INFO:continuum.py:2001:continuum_specfit: Fitting for the reddening took: 0.02 sec
  INFO:continuum.py:2007:continuum_specfit: Best-fitting spectroscopic A(V)=0.8598+/-0.0663
  INFO:continuum.py:2016:continuum_specfit: Solving for velocity dispersion: S/N_B=3.21, S/N_R=6.00, redshift=0.369
  INFO:continuum.py:2099:continuum_specfit: Fitting for the velocity dispersion took: 0.31 sec
  INFO:continuum.py:2109:continuum_specfit: Finding vdisp failed; adopting vdisp=150.00 km/s
  INFO:continuum.py:2169:continuum_specfit: Spectroscopic DN(4000)=0.950+/-0.028, Age=0.31 Gyr
  INFO:continuum.py:2230:continuum_specfit: Smooth continuum correction: b=-2.865%, r=0.754%, z=-1.058%
  INFO:fastspecfit.py:53:fastspec_one: Continuum-fitting object 0 [targetid 39633345008634465] took 1.11 sec
  INFO:emlines.py:1233:fit: Initial line-fitting with 28 free parameters took 0.79 sec (niter=119) with chi2=1.542
  INFO:emlines.py:1277:fit: Second (broad) line-fitting with 39 free parameters took 2.46 sec (niter=405) with chi2=1.633
  INFO:emlines.py:1286:fit: Chi2 with broad lines = 20.29041 and without broad lines = 13.94330
  INFO:emlines.py:1296:fit: All broad lines have been dropped, using narrow-line only model.
  INFO:fastspecfit.py:61:fastspec_one: Line-fitting object 0 [targetid=39633345008634465] took 3.51 sec
  INFO:fastspecfit.py:205:fastspec: Fitting everything took: 4.64 sec
  INFO:io.py:1019:write_fastspecfit: Writing results for 1 objects to fastspec.fits
  INFO:io.py:1074:write_fastspecfit: Writing out took 1.51 sec

See the :ref:`fastspec data model<fastspec datamodel>` for a full description of
the contents of the ``fastspec-example.fits`` file which is written out. We can
visualize the results by invoking the following command::

  $> fastspecfit-qa fastspec-example.fits
  
  INFO:io.py:984:read_fastspecfit: Read 1 object(s) from fastspec.fits
  INFO:continuum.py:137:__init__: Reading /global/cfs/cdirs/desi/science/gqp/templates/SSP-CKC14z/v1.0/SSP_Padova_CKC14z_Kroupa_Z0.0190.fits
  INFO:continuum.py:137:__init__: Reading /global/cfs/cdirs/desi/science/gqp/templates/SSP-CKC14z/v1.0/SSP_Padova_CKC14z_Kroupa_Z0.0190.fits
  INFO:fastspecfit-qa:102:main: Initializing the classes took: 1.43 sec
  INFO:io.py:296:select: Reading and parsing 1 unique redrockfile(s)
  INFO:io.py:348:select: specprod=fuji, coadd_type=healpix, survey=sv1, program=bright, healpix=7108
  INFO:io.py:443:select: Updating QSO redshifts using a QN threshold of 0.95.
  INFO:io.py:569:select: Gathered photometric metadata for 1 objects in 0.07 sec
  INFO:io.py:658:read_and_unpack: Reading 1 spectra from /global/cfs/cdirs/desi/spectro/redux/fuji/healpix/sv1/bright/71/7108/coadd-sv1-bright-7108.fits
  INFO:spectra.py:291:read_spectra: iotime 0.161 sec to read coadd-sv1-bright-7108.fits at 2022-08-07T04:01:18.246061
  INFO:continuum.py:836:get_linesigma: Forbidden masking sigma=98.374 km/s and S/N=73.893
  INFO:continuum.py:836:get_linesigma: Balmer masking sigma=88.367 km/s and S/N=84.792
  INFO:continuum.py:836:get_linesigma: UV/Broad masking sigma=0.000 km/s and S/N=0.000
  INFO:io.py:855:read_and_unpack: Read data for 1 objects in 0.62 sec
  INFO:emlines.py:2111:qa_fastspec: Writing ./fastspec-sv1-bright-7108-39633345008634465.png
  INFO:fastspecfit-qa:186:main: QA for everything took: 4.77 sec

.. image:: _static/fastspec-sv1-bright-7108-39633345008634465.png

This figure shows, from top to bottom: (top) the fit to the underlying stellar
continuum for each of the three DESI cameras (blue, green and red); (middle) the
fit to the (residual) emission-line spectrum after subtracting from the data the
best-fitting stellar continuum model and the smooth continuum correction (shown
as a light gray curve in the top panel); and (bottom) panels which zoom into all
the individual lines modeled by ``FastSpecFit``.

In some cases it may be convenient to generate your own figure of the data and
the best-fitting models, which you can do by reading the data yourself and using
the spectra stored in the ``MODELS`` FITS extension:

.. code-block:: python

  import numpy as np
  import fitsio 
  from astropy.table import Table
  import matplotlib.pyplot as plt
  
  from desiutil.dust import dust_transmission
  from desispec.io import read_spectra
  from desispec.coaddition import coadd_cameras
  
  specfile = '/global/cfs/cdirs/desi/spectro/redux/fuji/healpix/sv1/bright/71/7108/coadd-sv1-bright-7108.fits'
  fastfile = 'fastspec-example.fits'

  meta = Table(fitsio.read(fastfile, 'METADATA'))
  fast = Table(fitsio.read(fastfile, 'FASTSPEC'))
  
  models, hdr = fitsio.read(fastfile, 'MODELS', header=True)
  modelwave = hdr['CRVAL1'] + np.arange(hdr['NAXIS1']) * hdr['CDELT1']
  
  spec = read_spectra(specfile).select(targets=meta['TARGETID'])
  coadd_spec = coadd_cameras(spec)
  bands = coadd_spec.bands[0]
  
  mw_transmission_spec = dust_transmission(coadd_spec.wave[bands], meta['EBV'])
  
  fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
  ax1.plot(coadd_spec.wave[bands], coadd_spec.flux[bands].flatten() / mw_transmission_spec,
           color='gray', alpha=0.7, label='Data')
  ax1.plot(modelwave, models[0, 0, :], label='Stellar Continuum Model', ls='-', color='blue')
  ax1.plot(modelwave, models[0, 1, :], label='Smooth Continuum Correction', ls='--', color='k')
  ax1.set_ylim(-2.5, 7.5)
  ax1.legend(fontsize=8, loc='upper right')

  ax2.plot(coadd_spec.wave[bands], coadd_spec.flux[bands].flatten() / mw_transmission_spec,
           color='gray', alpha=0.7, label='Data')
  ax2.plot(modelwave, np.sum(models, axis=1).flatten(), label='Final Model', ls='-', color='red')
  ax2.legend(fontsize=8, loc='upper left')
  ax2.set_xlabel(r'Observed-frame Wavelength ($\AA$)')

  fig.subplots_adjust(hspace=0.05, top=0.95, right=0.95)
  fig.text(0.05, 0.5, r'Flux Density ($10^{-17}~{\rm erg}~{\rm s}^{-1}~{\rm cm}^{-2}~\AA^{-1}$)',
            ha='center', va='center', rotation='vertical')

  fig.savefig('fastspec-example.png')

.. image:: _static/fastspec-example.png

.. note::
   
   All the quantities and models returned by ``FastSpecFit`` are measured from
   DESI spectra which have been corrected for Galactic extinction, so the data
   have to be extinction-corrected when generating the figure above.

.. _`fastphot example`:

One fastphot Example
--------------------

``FastSpecFit`` can also model the broadband photometry (at the given DESI
redshift) using ``fastphot``. Using the same example object as above, we have::

  $> fastphot /global/cfs/cdirs/desi/spectro/redux/fuji/healpix/sv1/bright/71/7108/redrock-sv1-bright-7108.fits \
    --targetids 39633345008634465 --outfile fastphot-example.fits
    
  INFO:fastspecfit.py:123:parse: /global/homes/i/ioannis/code/desihub/fastspecfit/bin/fastphot /global/cfs/cdirs/desi/spectro/redux/fuji/healpix/sv1/bright/71/7108/redrock-sv1-bright-7108.fits --targetids 39633345008634465 --outfile fastphot-example.fits
  INFO:continuum.py:137:__init__: Reading /global/cfs/cdirs/desi/science/gqp/templates/SSP-CKC14z/v1.0/SSP_Padova_CKC14z_Kroupa_Z0.0190.fits
  INFO:fastspecfit.py:249:fastphot: Initializing the classes took: 1.28 sec
  INFO:io.py:296:select: Reading and parsing 1 unique redrockfile(s)
  INFO:io.py:348:select: specprod=fuji, coadd_type=healpix, survey=sv1, program=bright, healpix=7108
  INFO:io.py:443:select: Updating QSO redshifts using a QN threshold of 0.95.
  INFO:io.py:569:select: Gathered photometric metadata for 1 objects in 0.14 sec
  INFO:io.py:658:read_and_unpack: Reading 1 spectra from /global/cfs/cdirs/desi/spectro/redux/fuji/healpix/sv1/bright/71/7108/coadd-sv1-bright-7108.fits
  INFO:io.py:855:read_and_unpack: Read data for 1 objects in 0.08 sec
  INFO:fastspecfit.py:260:fastphot: Reading and unpacking the 1 spectra to be fitted took: 1.75 sec
  INFO:continuum.py:1815:continuum_fastphot: Preparing the models took 0.22 sec
  INFO:continuum.py:1843:continuum_fastphot: Fitting the photometry took: 0.05 sec
  INFO:continuum.py:1853:continuum_fastphot: Finding photometric A(V) failed; adopting A(V)=0.0000
  INFO:continuum.py:1896:continuum_fastphot: Photometric DN(4000)=1.170, Age=1.25 Gyr, Mr=-20.61 mag, Mstar=6.089e+09
  INFO:fastspecfit.py:83:fastphot_one: Continuum-fitting object 0 [targetid 39633345008634465] took 0.40 sec
  INFO:fastspecfit.py:276:fastphot: Fitting everything took: 0.41 sec
  INFO:io.py:1019:write_fastspecfit: Writing results for 1 objects to fastphot-example.fits
  INFO:io.py:1074:write_fastspecfit: Writing out took 0.11 sec

  $> fastspecfit-qa fastphot-example.fits
  
  INFO:fastspecfit-qa:44:parse: /global/homes/i/ioannis/code/desihub/fastspecfit/bin/fastspecfit-qa fastphot-example.fits
  INFO:io.py:984:read_fastspecfit: Read 1 object(s) from fastphot-example.fits
  INFO:continuum.py:137:__init__: Reading /global/cfs/cdirs/desi/science/gqp/templates/SSP-CKC14z/v1.0/SSP_Padova_CKC14z_Kroupa_Z0.0190.fits
  INFO:continuum.py:137:__init__: Reading /global/cfs/cdirs/desi/science/gqp/templates/SSP-CKC14z/v1.0/SSP_Padova_CKC14z_Kroupa_Z0.0190.fits
  INFO:fastspecfit-qa:102:main: Initializing the classes took: 1.95 sec
  INFO:io.py:296:select: Reading and parsing 1 unique redrockfile(s)
  INFO:io.py:348:select: specprod=fuji, coadd_type=healpix, survey=sv1, program=bright, healpix=7108
  INFO:io.py:443:select: Updating QSO redshifts using a QN threshold of 0.95.
  INFO:io.py:569:select: Gathered photometric metadata for 1 objects in 0.11 sec
  INFO:io.py:658:read_and_unpack: Reading 1 spectra from /global/cfs/cdirs/desi/spectro/redux/fuji/healpix/sv1/bright/71/7108/coadd-sv1-bright-7108.fits
  INFO:io.py:855:read_and_unpack: Read data for 1 objects in 0.07 sec
  INFO:continuum.py:2489:qa_fastphot: Writing ./fastphot-sv1-bright-7108-39633345008634465.png
  INFO:fastspecfit-qa:186:main: QA for everything took: 1.71 sec

.. image:: _static/fastphot-sv1-bright-7108-39633345008634465.png

Once again, please refer to the :ref:`fastphot data model<fastphot datamodel>`
for a full description of the contents of the ``fastphot-example.fits`` file.


.. note::
   
   The current version of ``FastSpecFit`` only models the rest-frame optical
   spectra of galaxies; there is no re-radiated dust emission. Consequently, in
   the figure above the *grzW1* photometric points which are used in the fit are
   shown using filled symbols while the open symbols (representing *W2*, *W3*,
   and *W4*) are not used in the fit.

.. _`production`:

More Examples
-------------

In the examples above, we selected one specific object using the ``--targetids``
optional input, which can also be a comma-separated list. For example::

  $> fastspec /global/cfs/cdirs/desi/spectro/redux/fuji/healpix/sv1/bright/71/7108/redrock-sv1-bright-7108.fits \
    --targetids 39633345008634465,39633334917139798,39633348330522913 \
    --outfile fastspec-example2.fits

Alternatively, you may want to fit a subset of the targets on this healpixel,
say the first 20 objects, in which case you would use the ``--ntargets`` keyword::

  $> fastspec /global/cfs/cdirs/desi/spectro/redux/fuji/healpix/sv1/bright/71/7108/redrock-sv1-bright-7108.fits \
    --ntargets 20 --outfile fastspec-example3.fits

If you don't want to start at the zeroth object, you can offset by an integer
number of targets using the ``--firsttarget`` option, which in this example
would fit objects 50 through 70::

  $> fastspec /global/cfs/cdirs/desi/spectro/redux/fuji/healpix/sv1/bright/71/7108/redrock-sv1-bright-7108.fits \
    --firsttarget 50 --ntargets 20 --outfile fastspec-example4.fits

Finally, when fitting more than one object, you probably want to use
multiprocessing, so that multiple objects are fit simultaneously. We can use
parallelism (assuming you're on a machine with more than one core) using the
``--mp`` input::

  $> fastspec /global/cfs/cdirs/desi/spectro/redux/fuji/healpix/sv1/bright/71/7108/redrock-sv1-bright-7108.fits \
    --firsttarget 50 --ntargets 20 --mp 20 --outfile fastspec-example5.fits

You can see all the options by calling either ``fastspec`` or ``fastphot`` with
the ``--help`` option, although most users will only invoke the options
documented above::

  $> fastspec --help
  usage: fastspec [-h] -o OUTFILE [--mp MP] [-n NTARGETS] [--firsttarget FIRSTTARGET] [--targetids TARGETIDS] [--solve-vdisp] [--ssptemplates SSPTEMPLATES]
                  [--mapdir MAPDIR] [--dr9dir DR9DIR] [--verbose]
                  [redrockfiles ...]
  
  positional arguments:
    redrockfiles          Full path to input redrock file(s). (default: None)
  
  optional arguments:
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
    --ssptemplates SSPTEMPLATES
                          Optional name of the SSP templates. (default: None)
    --mapdir MAPDIR       Optional directory name for the dust maps. (default: None)
    --dr9dir DR9DIR       Optional directory name for the DR9 photometry. (default: None)
    --verbose             Be verbose (for debugging purposes). (default: False)

What if you want to fit a particular survey, program, or healpixel. Do you
really need to specify the full path to each individual Redrock file? No!
``FastSpecFit`` knows how the DESI data are organized, but to access this
information we need to use the higher-level ``mpi-fastspecfit`` script. For
example, to fit all the objects in the *Fuji* spectroscopic production from
``survey=sv``, ``program=bright`` and ``healpix=7108``, we would do (here, on a
single interactive Perlmutter node)::

  $> salloc -N 1 -C cpu -A desi -t 00:10:00 --qos interactive -L cfs
  $> source /global/cfs/cdirs/desi/software/desi_environment.sh main
  $> module load fastspecfit/main
  $> export FASTSPECFIT_TEMPLATES=$DESI_ROOT/science/gqp/templates/SSP-CKC14z
  $> time mpi-fastspecfit --specprod fuji --survey sv1 --program bright \
    --healpix 7108 --mp 128 --outdir-data .
  $> ls -l ./fuji/healpix/sv1/bright/71/7108
  
  INFO:mpi.py:223:_findfiles: Building file list for survey=sv1 and program=bright
  INFO:mpi.py:309:plan: Found 1/1 redrockfiles (left) to do.
  INFO:mpi-fastspecfit:46:run_fastspecfit: Planning took 0.16 sec
  INFO:mpi-fastspecfit:96:run_fastspecfit: Rank 0, ntargets=264: fastspec /global/cfs/cdirs/desi/spectro/redux/fuji/healpix/sv1/bright/71/7108/redrock-sv1-bright-7108.fits -o ./fuji/healpix/sv1/bright/71/7108/fastspec-sv1-bright-7108.fits.gz --mp 128
  INFO:mpi-fastspecfit:119:run_fastspecfit:   rank 0 done in 113.58 sec
  INFO:mpi-fastspecfit:140:run_fastspecfit: All done at Sun Aug  7 06:17:02 2022
  
  real	1m55.770s
  user	14m38.856s
  sys	1m16.424s  

  total 12092
  -rw-rw-r-- 1 ioannis ioannis 12007670 Aug  7 06:17 fastspec-sv1-bright-7108.fits.gz
  -rw-rw-r-- 1 ioannis ioannis   370424 Aug  7 06:17 fastspec-sv1-bright-7108.log

Since fitting can be relatively expensive (in this case, it took about two
minutes to fit 264 targets with 128 cores), you may want to see what's going to
happen before fitting large numbers of objects, which we can do using the
``--plan`` and/or ``--dry-run`` options::

  $> mpi-fastspecfit --specprod fuji --survey sv1 --program bright \
    --healpix 7108 --outdir-data . --plan
    
  INFO:mpi.py:223:_findfiles: Building file list for survey=sv1 and program=bright
  INFO:mpi.py:309:plan: Found 1/1 redrockfiles (left) to do.

  $> mpi-fastspecfit --specprod fuji --survey sv1 --program bright \
    --healpix 7108 --outdir-data . --dry-run
    
  INFO:mpi.py:223:_findfiles: Building file list for survey=sv1 and program=bright
  INFO:mpi.py:309:plan: Found 1/1 redrockfiles (left) to do.
  INFO:mpi-fastspecfit:46:run_fastspecfit: Planning took 0.01 sec
  INFO:mpi-fastspecfit:96:run_fastspecfit: Rank 0, ntargets=264: fastspec /global/cfs/cdirs/desi/spectro/redux/fuji/healpix/sv1/bright/71/7108/redrock-sv1-bright-7108.fits -o ./fuji/healpix/sv1/bright/71/7108/fastspec-sv1-bright-7108.fits.gz --mp 128

If you leave off any combination of the ``--survey``, ``--program``, and/or
``--healpix`` options, the code will assume that you want all the possible
values of these keywords. For example, to see how many SV3 Redrock files would
need to be fit (not recommended without MPI parallelism!), one would do::

  $> mpi-fastspecfit --specprod fuji --survey sv3 --outdir-data . --plan
  INFO:mpi.py:223:_findfiles: Building file list for survey=sv3 and program=bright
  INFO:mpi.py:223:_findfiles: Building file list for survey=sv3 and program=dark
  INFO:mpi.py:223:_findfiles: Building file list for survey=sv3 and program=other
  INFO:mpi.py:223:_findfiles: Building file list for survey=sv3 and program=backup
  INFO:mpi.py:309:plan: Found 1023/1023 redrockfiles (left) to do.
  INFO:mpi.py:326:plan: Skipping 70 files with no targets.

.. note::  

  One must always specify the spectroscopic production when calling
  ``mpi-fastspecfit``, in this case ``--specprod fuji``. Also, to fit the
  broadband photometry instead of the DESI spectroscopy, simply call any of the
  examples in this section with the ``--fastphot`` option::

    $> mpi-fastspecfit --specprod fuji --survey sv3 --outdir-data . --plan --fastphot
    INFO:mpi.py:223:_findfiles: Building file list for survey=sv3 and program=bright
    INFO:mpi.py:223:_findfiles: Building file list for survey=sv3 and program=dark
    INFO:mpi.py:223:_findfiles: Building file list for survey=sv3 and program=other
    INFO:mpi.py:223:_findfiles: Building file list for survey=sv3 and program=backup
    INFO:mpi.py:309:plan: Found 1023/1023 redrockfiles (left) to do.
    INFO:mpi.py:326:plan: Skipping 70 files with no targets.

Finally, ``mpi-fastspecfit`` also knows about the tile-based *cumulative*,
*per-night*, and *per-exposure* coadds via the ``--coadd-type`` optional
input. For example::

  $> mpi-fastspecfit --specprod fuji --coadd-type cumulative --tile 80613 --outdir-data . --plan
  INFO:mpi.py:309:plan: Found 10/10 redrockfiles (left) to do.

  $> mpi-fastspecfit --specprod fuji --coadd-type pernight --tile 80613 --outdir-data . --plan
  INFO:mpi.py:309:plan: Found 57/57 redrockfiles (left) to do.
  
  $> mpi-fastspecfit --specprod fuji --coadd-type perexp --tile 80613 --outdir-data . --plan
  INFO:mpi.py:309:plan: Found 283/283 redrockfiles (left) to do.
