#!/usr/bin/env python
"""
fastspecfit.fastspecfit
=======================

FastSpec wrapper. Call with, e.g.,

  # nice BGS example
  fastspec /global/cfs/cdirs/desi/spectro/redux/everest/tiles/cumulative/80613/20210324/redrock-4-80613-thru20210324.fits -o fastspec.fits --targetids 39633345008634465 --specprod everest
  fastspec /global/cfs/cdirs/desi/spectro/redux/everest/healpix/sv1/bright/70/7022/redrock-sv1-bright-7022.fits -o fastspec2.fits --ntargets 1 --specprod everest

  # redrock is wrong!
  fastspec /global/cfs/cdirs/desi/spectro/redux/everest/tiles/cumulative/80605/redrock-0-80605-deep.fits -o fastspec.fits --targetids 39627652595714901

  # good test of needing smoothing continuum residuals before line-fitting
  fastspec /global/cfs/cdirs/desi/spectro/redux/everest/tiles/cumulative/80605/redrock-9-80605-deep.fits -o fastspec.fits --targetids 39627658622930703

  fastspec /global/cfs/cdirs/desi/spectro/redux/everest/tiles/cumulative/80613/redrock-0-80613-deep.fits -o fastspec.fits --targetids 39633314155332057
  fastspec /global/cfs/cdirs/desi/spectro/redux/everest/tiles/cumulative/80613/redrock-0-80606-deep.fits -o fastspec.fits --ntargets 2

  # nice QSO with broad lines
  fastspec /global/cfs/cdirs/desi/spectro/redux/everest/tiles/cumulative/80607/20201219/redrock-9-80607-thru20201219.fits -o fastspec2.fits --targetids 39633331528141827 --specprod everest
  fastspec /global/cfs/cdirs/desi/spectro/redux/everest/tiles/cumulative/80607/20201219/redrock-9-80607-thru20201219.fits -o fastspec2.fits --targetids 39633321176600909 --specprod everest
  fastspecfit-qa --fastspecfile ./fastspec2.fits -o cosmo-www/tmp/


ff = Table.read('/global/cfs/cdirs/desi/spectro/redux/everest/tiles/cumulative/80607/20201219/redrock-9-80607-thru20201219.fits', 'FIBERMAP')
tt = Table.read('/global/cfs/cdirs/desi/spectro/redux/everest/tiles/cumulative/80607/20201219/redrock-9-80607-thru20201219.fits', 'REDSHIFTS')

In [16]: tt[(tt['Z'] > 2) * (ff['COADD_FIBERSTATUS'] == 0)]
Out[16]:
<Table length=19>
     TARGETID             CHI2                  COEFF [10]                  Z                   ZERR          ZWARN NPIXELS SPECTYPE SUBTYPE NCOEFF      DELTACHI2
      int64             float64                  float64                 float64              float64         int64  int64   bytes6  bytes20 int64        float64
------------------ ------------------ ----------------------------- ------------------ ---------------------- ----- ------- -------- ------- ------ -------------------
 39633321180792100 208156.92278671265 0.00036161256944187563 .. 0.0  4.771305349553081   2.67536442758405e-05     0    7930      QSO              4  26860.353954315186
 39633321176600909 16523.136709213257  8.831772363315664e-05 .. 0.0 2.2657412125056293  0.0001064761302343068     0    7930      QSO              4   10004.88240519166
 39633321180791247  12710.44403065741 2.2275734055524796e-05 .. 0.0   2.79282043898751  0.0003073365407085948     4    7930      QSO              4 0.31629960238933563
 39633321180791630 123169.91798686981  0.0016733355408052341 .. 0.0  2.625127767580395 4.8572220298588665e-05     0    7930      QSO              4  177285.64080905914
 39633328097199918 15715.593243002892  4.814717482059453e-05 .. 0.0  2.263183148835255  0.0001787750120986952     0    7929      QSO              4  1742.2603384945542
616094076948709456  9683.685187131166 1.6850354871158607e-07 .. 0.0  2.018166042597203  0.0002754084591784837     5    7930      QSO              4  6.7541270181536674
 39633331528143479 170210.14914035797  0.0005605767003717567 .. 0.0 2.9467324291846437  3.677543414159961e-05     0    7930      QSO              4  309962.47056770325
 39633324653676573  12008.00667237863  7.001341136996581e-05 .. 0.0   2.03741457380998  0.0003714630610027555     0    7929      QSO              4   784.9251936562359
 39633324653676217  11427.49386062473 3.5152059216539954e-05 .. 0.0 2.0325382432226378  0.0003267446617409548     0    7929      QSO              4   467.0210649073124
 39633321176598475 14995.089814335108  7.537447147364333e-05 .. 0.0 2.4066543500367796 0.00014704600656576563     0    7930      QSO              4   4943.020963191986
 39633331532336872 14732.979373201728 0.00023425889404080574 .. 0.0 2.0153782188148104 0.00023455017770337706     0    7929      QSO              4  7651.7109602838755
 39633328097202474 19191.189551770687  0.0003443592319551414 .. 0.0 2.2344579631752475 0.00011693256745068771     0    7930      QSO              4    38755.5963024199
616094080404816345   9997.95639693737  6.634271448721053e-06 .. 0.0 2.0026654995655147 0.00019580992192198883     5    7929      QSO              4   8.437787368893623
616094080400622381  9627.908837988973 1.8039854586058742e-06 .. 0.0 2.0322822167365344 0.00023962923329210812     5    7929      QSO              4  3.9985851272940636
 39633321189179675 12645.148697257042  6.765447744035747e-05 .. 0.0 2.8952281257102097  0.0001885014403313145     0    7929      QSO              4   4907.485580228269
 39633321176598312  61853.88678389788  0.0006109744669765684 .. 0.0   2.53071304738548 0.00013073979474828362     0    7929      QSO              4   48060.86926621199
 39633328097199634 29596.591431498528 0.00020488152277961015 .. 0.0 2.0451104198087444  8.248902899025006e-05     0    7930      QSO              4   22535.68895442784
 39633321172405812  13977.65833118558  7.064589556932252e-05 .. 0.0 2.2761846709951663 0.00015072825430395817     0    7930      QSO              4    4329.20822699368
 39633321172405061 61947.770018577576  0.0013760717996423705 .. 0.0 2.0597296941700125  4.091995977008691e-05     0    7922      QSO              4  131703.95274591446  


Fastphot wrapper. Call with, e.g.,
  fastphot /global/cfs/cdirs/desi/spectro/redux/everest/tiles/perexp/80607/00068028/redrock-0-80607-exp00068028.fits -o fastphot-perexp.fits --ntargets 2 --specprod everest
  fastphot /global/cfs/cdirs/desi/spectro/redux/everest/tiles/pernight/80607/20201218/redrock-0-80607-20201218.fits -o fastphot-pernight.fits --ntargets 2 --specprod everest
  fastphot /global/cfs/cdirs/desi/spectro/redux/everest/tiles/cumulative/80607/20201219/redrock-0-80607-thru20201219.fits -o fastphot-thrunight.fits --ntargets 2 --specprod everest
  fastphot /global/cfs/cdirs/desi/spectro/redux/everest/healpix/sv1/bright/70/7022/redrock-sv1-bright-7022.fits -o fastphot-hpx.fits --ntargets 2 --specprod everest

"""
import pdb # for debugging

import os, time
import numpy as np

from desiutil.log import get_logger
log = get_logger()

## ridiculousness! - this seems to come from healpy, blarg
#import tempfile
#os.environ['MPLCONFIGDIR'] = tempfile.mkdtemp()

def _fastspec_one(args):
    """Multiprocessing wrapper."""
    return fastspec_one(*args)

def fastspec_one(iobj, data, out, meta, CFit, EMFit, solve_vdisp=False):
    """Fit one spectrum."""
    #log.info('Continuum-fitting object {}'.format(iobj))
    t0 = time.time()

    cfit, continuummodel, smooth_continuum = CFit.continuum_specfit(data, solve_vdisp=solve_vdisp)
    for col in cfit.colnames:
        out[col] = cfit[col]

    ## Copy over the reddening-corrected fluxes -- messy!
    #for iband, band in enumerate(CFit.fiber_bands):
    #    meta['FIBERTOTFLUX_{}'.format(band.upper())] = data['fiberphot']['nanomaggies'][iband]
    #    #result['FIBERTOTFLUX_IVAR_{}'.format(band.upper())] = data['fiberphot']['nanomaggies_ivar'][iband]
    #for iband, band in enumerate(CFit.bands):
    #    meta['FLUX_{}'.format(band.upper())] = data['phot']['nanomaggies'][iband]
    #    meta['FLUX_IVAR_{}'.format(band.upper())] = data['phot']['nanomaggies_ivar'][iband]
        
    log.info('Continuum-fitting object {} [targetid {}] took {:.2f} sec'.format(
        iobj, meta['TARGETID'], time.time()-t0))

    # Fit the emission-line spectrum.
    t0 = time.time()
    emfit = EMFit.fit(data, continuummodel, smooth_continuum)
    for col in emfit.colnames:
        out[col] = emfit[col]
    log.info('Line-fitting object {} (targetid={}) took {:.2f} sec'.format(
        iobj, meta['TARGETID'], time.time()-t0))

    return out, meta

def _fastphot_one(args):
    """Multiprocessing wrapper."""
    return fastphot_one(*args)

def fastphot_one(iobj, data, out, meta, CFit):
    """Fit one spectrum."""
    #log.info('Continuum-fitting object {}'.format(iobj))
    t0 = time.time()
    cfit, _ = CFit.continuum_fastphot(data)
    for col in cfit.colnames:
        out[col] = cfit[col]

    # Copy over the reddening-corrected fluxes -- messy!
    for iband, band in enumerate(CFit.fiber_bands):
        meta['FIBERTOTFLUX_{}'.format(band.upper())] = data['fiberphot']['nanomaggies'][iband]
        #result['FIBERTOTFLUX_IVAR_{}'.format(band.upper())] = data['fiberphot']['nanomaggies_ivar'][iband]
    for iband, band in enumerate(CFit.bands):
        meta['FLUX_{}'.format(band.upper())] = data['phot']['nanomaggies'][iband]
        meta['FLUX_IVAR_{}'.format(band.upper())] = data['phot']['nanomaggies_ivar'][iband]
        
    log.info('Continuum-fitting object {} took {:.2f} sec'.format(iobj, time.time()-t0))
    
    return out, meta

def parse(options=None):
    """Parse input arguments.

    """
    import argparse, sys

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-n', '--ntargets', type=int, help='Number of targets to process in each file.')
    parser.add_argument('--firsttarget', type=int, default=0, help='Index of first object to to process in each file (0-indexed).')
    parser.add_argument('--targetids', type=str, default=None, help='Comma-separated list of target IDs to process.')
    parser.add_argument('--mp', type=int, default=1, help='Number of multiprocessing processes per MPI rank or node.')
    parser.add_argument('-o', '--outfile', type=str, required=True, help='Full path to output filename.')

    parser.add_argument('--solve-vdisp', action='store_true', help='Solve for the velocity dispersion (only when using fastspec).')

    # make specprod required until this ticket is addressed--
    # https://github.com/desihub/desispec/issues/1077
    parser.add_argument('--specprod', type=str, default='everest', choices=['everest', 'denali', 'daily'],
                        required=True, help='Spectroscopic production.')
    #parser.add_argument('--coadd-type', type=str, default='healpix', choices=['healpix', 'cumulative', 'pernight', 'perexp'],
    #                    help='Type of spectral coadds corresponding to the input redrockfiles.')

    parser.add_argument('redrockfiles', nargs='*', help='Full path to input redrock file(s).')

    if options is None:
        args = parser.parse_args()
        log.info(' '.join(sys.argv))
    else:
        args = parser.parse_args(options)
        log.info('fastspec {}'.format(' '.join(options)))

    return args

def fastspec(args=None, comm=None):
    """Main fastspec module.

    """
    from astropy.table import Table
    from fastspecfit.continuum import ContinuumFit
    from fastspecfit.emlines import EMLineFit
    from fastspecfit.io import DESISpectra, write_fastspecfit

    if isinstance(args, (list, tuple, type(None))):
        args = parse(args)

    if args.targetids:
        targetids = [int(x) for x in args.targetids.split(',')]
    else:
        targetids = args.targetids

    # Initialize the continuum- and emission-line fitting classes.
    t0 = time.time()
    CFit = ContinuumFit()
    EMFit = EMLineFit()
    Spec = DESISpectra(specprod=args.specprod)
    log.info('Initializing the classes took: {:.2f} sec'.format(time.time()-t0))

    # Read the data.
    t0 = time.time()

    Spec.find_specfiles(args.redrockfiles, firsttarget=args.firsttarget,
                        targetids=targetids, ntargets=args.ntargets)
    if len(Spec.specfiles) == 0:
        return
    data = Spec.read_and_unpack(CFit, fastphot=False, synthphot=True, remember_coadd=True)

    out, meta = Spec.init_output(CFit=CFit, EMFit=EMFit, fastphot=False)
    log.info('Reading and unpacking the {} spectra to be fitted took: {:.2f} sec'.format(
        Spec.ntargets, time.time()-t0))

    # Fit in parallel
    t0 = time.time()
    fitargs = [(iobj, data[iobj], out[iobj], meta[iobj], CFit, EMFit, args.solve_vdisp)
               for iobj in np.arange(Spec.ntargets)]
    if args.mp > 1:
        import multiprocessing
        with multiprocessing.Pool(args.mp) as P:
            _out = P.map(_fastspec_one, fitargs)
    else:
        _out = [fastspec_one(*_fitargs) for _fitargs in fitargs]
    _out = list(zip(*_out))
    out = Table(np.hstack(_out[0]))
    meta = Table(np.hstack(_out[1]))
    log.info('Fitting everything took: {:.2f} sec'.format(time.time()-t0))

    # Write out.
    write_fastspecfit(out, meta, outfile=args.outfile, specprod=Spec.specprod,
                      coadd_type=Spec.coadd_type, fastphot=False)

def fastphot(args=None, comm=None):
    """Main fastphot module.

    """
    from astropy.table import Table
    from fastspecfit.continuum import ContinuumFit
    from fastspecfit.io import DESISpectra, write_fastspecfit

    if isinstance(args, (list, tuple, type(None))):
        args = parse(args)

    if args.targetids:
        targetids = [int(x) for x in args.targetids.split(',')]
    else:
        targetids = args.targetids

    # Initialize the continuum-fitting classes.
    t0 = time.time()
    CFit = ContinuumFit()
    Spec = DESISpectra(specprod=args.specprod)
    log.info('Initializing the classes took: {:.2f} sec'.format(time.time()-t0))

    # Read the data.
    t0 = time.time()
    Spec.find_specfiles(args.redrockfiles, firsttarget=args.firsttarget,
                        targetids=targetids, ntargets=args.ntargets)
    if len(Spec.specfiles) == 0:
        return
    data = Spec.read_and_unpack(CFit, fastphot=True, synthphot=False)
    
    out, meta = Spec.init_output(CFit=CFit, fastphot=True)
    log.info('Reading and unpacking the {} spectra to be fitted took: {:.2f} sec'.format(
        Spec.ntargets, time.time()-t0))

    # Fit in parallel
    t0 = time.time()
    fitargs = [(iobj, data[iobj], out[iobj], meta[iobj], CFit)
               for iobj in np.arange(Spec.ntargets)]
    if args.mp > 1:
        import multiprocessing
        with multiprocessing.Pool(args.mp) as P:
            _out = P.map(_fastphot_one, fitargs)
    else:
        _out = [fastphot_one(*_fitargs) for _fitargs in fitargs]
    _out = list(zip(*_out))
    out = Table(np.hstack(_out[0]))
    meta = Table(np.hstack(_out[1]))
    log.info('Fitting everything took: {:.2f} sec'.format(time.time()-t0))

    # Write out.
    write_fastspecfit(out, meta, outfile=args.outfile, specprod=Spec.specprod,
                      coadd_type=Spec.coadd_type, fastphot=True)
