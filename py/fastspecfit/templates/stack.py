"""
fastspecfit.templates.stack
===========================

Tools for generating stacked spectra.

"""
import pdb # for debugging

import os, time
import numpy as np

from desiutil.log import get_logger
log = get_logger()

def _stack_onebin(args):
    """Multiprocessing wrapper."""
    return deredshift_one(*args)

def stack_onebin(flux2d, ivar2d, cflux2d, templatewave, continuumwave):
    """Build one stack."""

    templatewave1, templateflux1, templateivar1, _, templatepix = iterative_stack(
        templatewave, flux2d, ivar2d, constant_ivar=False, verbose=False)
    
    _, continuumflux1, _, _, templatecpix = iterative_stack(continuumwave, cflux2d, verbose=False)

    return templateflux1, templateivar1, templatepix, continuumflux1, templatecpix

def write_binned_stacks(metadata, wave, flux, ivar, cwave, cflux, outfile):
    """Write out the stacked spectra.
    
    """
    from astropy.io import fits
    
    hduflux = fits.PrimaryHDU(flux.astype('f4'))
    hduflux.header['EXTNAME'] = 'FLUX'

    hduivar = fits.ImageHDU(ivar.astype('f4'))
    hduivar.header['EXTNAME'] = 'IVAR'

    hduwave = fits.ImageHDU(wave.astype('f8'))
    hduwave.header['EXTNAME'] = 'WAVE'
    hduwave.header['BUNIT'] = 'Angstrom'
    hduwave.header['AIRORVAC'] = ('vac', 'vacuum wavelengths')

    hducflux = fits.ImageHDU(cflux.astype('f4'))
    hducflux.header['EXTNAME'] = 'CFLUX'

    hducwave = fits.ImageHDU(cwave.astype('f8'))
    hducwave.header['EXTNAME'] = 'CWAVE'
    hducwave.header['BUNIT'] = 'Angstrom'
    hducwave.header['AIRORVAC'] = ('vac', 'vacuum wavelengths')
    
    hdumeta = fits.convenience.table_to_hdu(metadata)
    hdumeta.header['EXTNAME'] = 'METADATA'
    #hdumeta.header['SPECPROD'] = (specprod, 'spectroscopic production name')

    hx = fits.HDUList([hduflux, hduivar, hduwave, hducflux, hducwave, hdumeta])
    hx.writeto(outfile, overwrite=True, checksum=True)
    log.info('Writing {} stacked spectra to {}'.format(len(metadata), outfile))

def quick_stack(wave, flux2d, ivar2d=None, constant_ivar=False):
    """Simple inverse-variance-weighted stack.
    
    """
    if ivar2d is None:
        ivar2d = np.ones_like(flux2d)
    
    if constant_ivar:
        _ivar2d = (ivar2d > 0) * 1.0 # * (flux2d > 0)
    else:
        _ivar2d = ivar2d #* (flux2d > 0)
        
    ivar = np.sum(_ivar2d, axis=0)
    #nperpix2 = np.sum((_ivar2d > 0) * (flux2d > 0), axis=0).astype(int)
    nperpix = np.sum((_ivar2d > 0), axis=0).astype(int)

    good = np.where((ivar > 0) * (nperpix > 0.9*np.max(nperpix)))[0]
    flux = np.sum(_ivar2d[:, good] * flux2d[:, good], axis=0) / ivar[good]
    
    #pos = np.where(flux > 0)[0]
    #good = good[pos]
    #flux = flux[pos]
    #ivar = ivar[pos]

    return wave[good], flux, ivar, nperpix[good], good

def iterative_stack(wave, flux2d, ivar2d=None, constant_ivar=False, maxdiff=0.01, 
                    maxiter=500, normwave=4500, smooth=None, debug=False, verbose=True):
    """Iterative stacking algorithm taken from Bovy, Hogg, & Moustakas 2008.
       https://arxiv.org/pdf/0805.1200.pdf

    """
    from scipy.ndimage import median_filter
    
    ngal, npix = flux2d.shape
    if ivar2d is None:
        ivar2d = np.ones_like(flux2d)

    # initial template, no relative scaling
    templatewave, templateflux, templateivar, nperpix, goodpix = quick_stack(
        wave, flux2d, ivar2d, constant_ivar=constant_ivar)

    _flux2d = flux2d[:, goodpix]
    _ivar2d = ivar2d[:, goodpix]
    
    for ii in np.arange(maxiter):
        # compute the maximum likelihood scale factor.
        scale = np.sum(templateflux[np.newaxis, :] * _flux2d, axis=1) / np.sum(_flux2d * _flux2d, axis=1)
        #print(scale)
        
        _flux2d *= scale[:, np.newaxis]
        _ivar2d /= scale[:, np.newaxis]**2

        # now update the template
        _templatewave, _templateflux, _templateivar, _, _ = quick_stack(
            wave[goodpix], _flux2d, _ivar2d, constant_ivar=constant_ivar)
    
        diff = _templateflux - templateflux
    
        if ii % 2 == 0:
            #print(ii, np.median(diff), np.max(np.abs(diff)))
            if debug:
                if smooth:
                    for jj in np.arange(ngal):
                        #plt.plot(wave[goodpix], flux2d[jj, goodpix], alpha=0.5, color='gray')
                        plt.plot(wave[goodpix], median_filter(_flux2d[jj, :], smooth), alpha=0.5, color='gray')
                    plt.plot(templatewave, median_filter(templateflux, smooth), color='k', lw=2)
                    plt.plot(_templatewave, median_filter(_templateflux, smooth), color='orange', alpha=0.5)
                else:
                    for jj in np.arange(ngal):
                        #plt.plot(wave[goodpix], flux2d[jj, goodpix], alpha=0.5, color='gray')
                        plt.plot(wave[goodpix], _flux2d[jj, :], alpha=0.5, color='gray')
                    plt.plot(templatewave, templateflux, color='k', lw=2)
                    plt.plot(_templatewave, _templateflux, color='orange', alpha=0.5)
                    
                #plt.xlim(6400, 6800)
                plt.xlim(3900, 6500)
                #plt.ylim(np.min(templateflux), np.max(templateflux))
                plt.show()
        
        templateflux = _templateflux
        templateivar = _templateivar
        
        _maxdiff = np.max(np.abs(diff))
        if _maxdiff < maxdiff or ii == maxiter - 1:
            if ii == maxiter - 1:
                msg = 'Did not converge'
            else:
                msg = 'Converged'
            if verbose:
                print('{} after {}/{} iterations with a maximum difference of {:.4f}.'.format(
                    msg, ii, maxiter, _maxdiff))
            
            # normalize
            normflux = np.median(templateflux[(templatewave > (normwave-10)) * (templatewave < (normwave+10))])
            #normflux = np.interp(normwave, templatewave, templateflux)
            templateflux /= normflux
            templateivar *= normflux**2
            
            break

    return templatewave, templateflux, templateivar, nperpix, goodpix

def stack_in_bins(sample, data, templatewave, mp=1, continuumwave=None, outfile=None):
    """Stack spectra in bins given the output of a spectra_in_bins run.

    """
    npix, ntemplate = len(templatewave), len(sample)
    templateflux = np.zeros((ntemplate, npix), dtype='f4')
    templateivar = np.zeros((ntemplate, npix), dtype='f4')

    if continuumwave is not None:
        ncpix = len(continuumwave)
        continuumflux = np.zeros((ntemplate, ncpix), dtype='f4')

    # parallelize the stack-building step
    mpargs = []
    for samp in sample:
        ibin = samp['ibin'].astype(str)
        mpargs.append([data[ibin]['flux'], data[ibin]['ivar'], data[ibin]['cflux'], templatewave, continuumwave])

    t0 = time.time()
    if mp > 1:
        import multiprocessing
        with multiprocessing.Pool(args.mp) as P:
            results = P.map(_stack_onebin, mpargs)
    else:
        results = [stack_onebin(*_mpargs) for _mpargs in mpargs]
    log.info('Building all stacks took: {:.2f} min'.format((time.time()-t0) / 60))

    # unpack the results and (optionally) write out
    results = list(zip(*results))
    for itemp in np.arange(ntemplate):
        templatepix = results[2][itemp]
        templatecpix = results[4][itemp]
        
        templateflux[itemp, templatepix] = results[0][itemp]
        templateivar[itemp, templatepix] = results[1][itemp]
        continuumflux[itemp, templatecpix] = results[3][itemp]
    
    if outfile:
        write_binned_stacks(sample, templatewave, templateflux, templateivar, 
                            continuumwave, continuumflux, outfile)

