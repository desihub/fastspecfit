"""
fastspecfit.templates.stack
===========================

Tools for generating stacked spectra.

"""
import pdb # for debugging

import os, time
import numpy as np
from astropy.table import Table, Column

from desiutil.log import get_logger
log = get_logger()

def _fastspec_onestack(args):
    """Multiprocessing wrapper."""
    return fastspec_onestack(*args)

def fastspec_onestack(fastmeta, fastspec, flux, ivar, meta,
                      wave, CFit, EMFit, qadir, qaprefix):
    """Fit one stacked spectrum.

    """
    from scipy.sparse import identity
    from desispec.resolution import Resolution

    good = np.where(ivar > 0)[0]
    npix = len(good)
    if npix == 0:
        log.warning('No good pixels!')
        raise ValueError

    # mock up a fastspec-compatible dictionary of DESI data
    data = {}
    data['zredrock'] = fastspec['CONTINUUM_Z']
    data['cameras'] = ['all']
    data['photsys'] = 'S'
    data['wave'] = [wave[good] * (1+data['zredrock'])]
    data['flux'] = [flux[good]]
    data['ivar'] = [ivar[good]]
    data['res'] = [Resolution(identity(n=npix))] # hack!
    data['snr'] = [np.median(flux[good] * np.sqrt(ivar[good]))]

    linemask = np.ones(npix, bool)
    if False:
        for line in CFit.linetable:
            zwave = line['restwave'] * (1 + data['zredrock'])
            if line['isbroad']:
                sigma = CFit.linemask_sigma_broad
            else:
                sigma = CFit.linemask_sigma
            I = np.where((data['wave'][0] >= (zwave - 1.5*sigma * zwave / C_LIGHT)) *
                         (data['wave'][0] <= (zwave + 1.5*sigma * zwave / C_LIGHT)))[0]
            if len(I) > 0:
                linemask[I] = False  # False = affected by line
    data['linemask'] = [linemask]

    # adapt the emission-line fitting range for each object.
    minwave, maxwave = data['wave'][0].min()-1.0, data['wave'][0].max()+1.0
    EMFit.log10wave = np.arange(np.log10(minwave), np.log10(maxwave), EMFit.dlogwave)

    # fit the continuum and the (residual) emission-line spectrum
    cfit, continuummodel, smooth_continuum = CFit.continuum_specfit(data)    
    emfit = EMFit.fit(data, continuummodel, smooth_continuum, synthphot=False)

    # pack up the results and write out
    for col in cfit.colnames:
        if col in fastspec.colnames:
            fastspec[col] = cfit[col]
    for col in emfit.colnames:
        fastspec[col] = emfit[col]

    # optional?
    #EMFit.qa_fastspec(data, fastspec, fastmeta, wavelims=(minwave, maxwave),
    #                  outdir=qadir, outprefix=qaprefix)

    return fastspec, fastmeta
    
def _stack_onebin(args):
    """Multiprocessing wrapper."""
    return stack_onebin(*args)

def stack_onebin(flux2d, ivar2d, templatewave, cflux2d, continuumwave):
    """Build one stack."""

    templatewave1, templateflux1, templateivar1, _, templatepix = iterative_stack(
        templatewave, flux2d, ivar2d, constant_ivar=False, verbose=False)
    
    if continuumwave is None:
        return templateflux1, templateivar1, templatepix
    else:
        _, continuumflux1, _, _, templatecpix = iterative_stack(continuumwave, cflux2d, verbose=False)
        return templateflux1, templateivar1, templatepix, continuumflux1, templatecpix

def write_binned_stacks(outfile, wave, flux, ivar, metadata=None, cwave=None,
                        cflux=None, fastspec=None, fastmeta=None):
    """Write out the stacked spectra.
    
    """
    from astropy.io import fits

    nobj, _ = flux.shape
    
    hduflux = fits.PrimaryHDU(flux.astype('f4'))
    hduflux.header['EXTNAME'] = 'FLUX'

    hduivar = fits.ImageHDU(ivar.astype('f4'))
    hduivar.header['EXTNAME'] = 'IVAR'

    hduwave = fits.ImageHDU(wave.astype('f8'))
    hduwave.header['EXTNAME'] = 'WAVE'
    hduwave.header['BUNIT'] = 'Angstrom'
    hduwave.header['AIRORVAC'] = ('vac', 'vacuum wavelengths')

    hdulist = [hduflux, hduivar, hduwave]

    if cflux is not None and cwave is not None:
        hducflux = fits.ImageHDU(cflux.astype('f4'))
        hducflux.header['EXTNAME'] = 'CFLUX'

        hducwave = fits.ImageHDU(cwave.astype('f8'))
        hducwave.header['EXTNAME'] = 'CWAVE'
        hducwave.header['BUNIT'] = 'Angstrom'
        hducwave.header['AIRORVAC'] = ('vac', 'vacuum wavelengths')

        hdulist = hdulist + [hducflux, hducwave]
        
    if metadata is not None and fastmeta is None:
        hdumetadata = fits.convenience.table_to_hdu(metadata)
        hdumetadata.header['EXTNAME'] = 'METADATA'
        #hdumetadata.header['SPECPROD'] = (specprod, 'spectroscopic production name')
        hdulist = hdulist + [hdumetadata]

    if metadata is None and fastmeta is not None:
        hdufastmeta = fits.convenience.table_to_hdu(fastmeta)
        hdufastmeta.header['EXTNAME'] = 'METADATA'
        hdulist = hdulist + [hdufastmeta]

    if fastspec is not None:
        hdufast = fits.convenience.table_to_hdu(fastspec)
        hdufast.header['EXTNAME'] = 'FASTSPEC'
        hdulist = hdulist + [hdufast]
    
    hx = fits.HDUList(hdulist)

    hx.writeto(outfile, overwrite=True, checksum=True)
    log.info('Writing {} spectra to {}'.format(nobj, outfile))

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
    
    nobj, npix = flux2d.shape
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
                    for jj in np.arange(nobj):
                        #plt.plot(wave[goodpix], flux2d[jj, goodpix], alpha=0.5, color='gray')
                        plt.plot(wave[goodpix], median_filter(_flux2d[jj, :], smooth), alpha=0.5, color='gray')
                    plt.plot(templatewave, median_filter(templateflux, smooth), color='k', lw=2)
                    plt.plot(_templatewave, median_filter(_templateflux, smooth), color='orange', alpha=0.5)
                else:
                    for jj in np.arange(nobj):
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

def stack_in_bins(sample, data, templatewave, mp=1, continuumwave=None, stackfile=None):
    """Stack spectra in bins given the output of a spectra_in_bins run.

    See bin/desi-templates for how to generate the required input.

    """
    npix, ntemplate = len(templatewave), len(sample)
    templateflux = np.zeros((ntemplate, npix), dtype='f4')
    templateivar = np.zeros((ntemplate, npix), dtype='f4')

    if continuumwave is not None:
        continuumflux = np.zeros((ntemplate, len(continuumwave)), dtype='f4')

    # parallelize the stack-building step
    mpargs = []
    for samp in sample:
        ibin = samp['IBIN'].astype(str)
        if continuumwave is not None:
            mpargs.append([data[ibin]['flux'], data[ibin]['ivar'], templatewave, data[ibin]['cflux'], continuumwave])
        else:
            mpargs.append([data[ibin]['flux'], data[ibin]['ivar'], templatewave])

    t0 = time.time()
    if mp > 1:
        import multiprocessing
        with multiprocessing.Pool(mp) as P:
            results = P.map(_stack_onebin, mpargs)
    else:
        results = [stack_onebin(*_mpargs) for _mpargs in mpargs]
    log.info('Building all stacks took: {:.2f} min'.format((time.time()-t0) / 60))

    # unpack the results and (optionally) write out
    results = list(zip(*results))
    for itemp in np.arange(ntemplate):
        templatepix = results[2][itemp]
        templateflux[itemp, templatepix] = results[0][itemp]
        templateivar[itemp, templatepix] = results[1][itemp]

        if continuumwave is not None:
            templatecpix = results[4][itemp]
            continuumflux[itemp, templatecpix] = results[3][itemp]
    
    if stackfile:
        write_binned_stacks(stackfile, templatewave, templateflux, templateivar,
                            metadata=sample, cwave=continuumwave, cflux=continuumflux)

def fastspecfit_stacks(stackfile, mp=1, qadir=None, qaprefix=None, fastspecfile=None):
    """Model stacked spectra using fastspecfit.

    Called by bin/desi-templates.

    stackfile is the output of stack_in_bins

    """
    import fitsio
    from fastspecfit.continuum import ContinuumFit
    from fastspecfit.emlines import EMLineFit
    from fastspecfit.util import C_LIGHT    
    #rom fastspecfit.io import write_fastspecfit

    if not os.path.isfile(stackfile):
        log.warning('Stack file {} not found!'.format(stackfile))
        raise IOError
    
    meta = Table(fitsio.read(stackfile, ext='METADATA'))
    wave = fitsio.read(stackfile, ext='WAVE')
    flux = fitsio.read(stackfile, ext='FLUX')
    ivar = fitsio.read(stackfile, ext='IVAR')
    #cwave = fitsio.read(stackfile, ext='CWAVE')
    #cflux = fitsio.read(stackfile, ext='CFLUX')

    nobj = len(meta)
    log.info('Read {} stacked spectra from {}'.format(nobj, stackfile))

    # initialize the continuum-fitting class and the output tables
    CFit = ContinuumFit()
    EMFit = EMLineFit()

    fastmeta = meta.copy()
    for col in fastmeta.colnames:
        fastmeta.rename_column(col, col.upper())
    fastmeta.add_column(Column(name='FIBER', data=np.arange(nobj, dtype=np.int32)), index=0)
    fastmeta.add_column(Column(name='TILEID', data=np.ones(nobj, dtype=np.int32)), index=0)
    fastmeta.add_column(Column(name='TARGETID', data=np.arange(nobj, dtype=np.int64)), index=0)

    from astropy.table import hstack
    fastspec = hstack((CFit.init_spec_output(nobj=nobj), EMFit.init_output(CFit.linetable, nobj=nobj)))
    fastspec.rename_column('CONTINUUM_SNR_B', 'CONTINUUM_SNR_ALL')
    fastspec.remove_columns(['CONTINUUM_SNR_R', 'CONTINUUM_SNR_Z'])
    fastspec.add_column(Column(name='TARGETID', data=np.arange(nobj, dtype=np.int64)), index=0)

    fastspec['CONTINUUM_Z'] = meta['ZOBJ']

    # Fit in parallel
    t0 = time.time()
    fitargs = [(fastmeta[iobj], fastspec[iobj], flux[iobj, :], ivar[iobj, :], meta[iobj],
                wave, CFit, EMFit, qadir, qaprefix) for iobj in np.arange(nobj)]
    if mp > 1:
        import multiprocessing
        with multiprocessing.Pool(mp) as P:
            _out = P.map(_fastspec_onestack, fitargs)
    else:
        _out = [fastspec_onestack(*_fitargs) for _fitargs in fitargs]
    _out = list(zip(*_out))
    fastspec = Table(np.hstack(_out[0])) # overwrite
    fastmeta = Table(np.hstack(_out[1]))
    log.info('Fitting everything took: {:.2f} sec'.format(time.time()-t0))

    if fastspecfile:
        write_binned_stacks(fastspecfile, wave, flux, ivar,
                            fastspec=fastspec, fastmeta=fastmeta)
