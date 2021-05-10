"""
fastspecfit.templates.templates
===============================

Tools for generating stacked spectra and building templates.

"""
import pdb # for debugging

import os, time
import multiprocessing
import numpy as np
import fitsio
from astropy.table import Table, Column

from fastspecfit.util import C_LIGHT

from desiutil.log import get_logger
log = get_logger()

def _rebuild_fastspec_spectrum(args):
    """Multiprocessing wrapper."""
    return rebuild_fastspec_spectrum(*args)

def rebuild_fastspec_spectrum(fastspec, wave, flux, ivar, CFit, EMFit,
                              full_resolution=False, emlines_only=False):
    """Rebuilding a single fastspec model continuum and emission-line spectrum given
    a fastspec astropy.table.Table.

    full_resolution - returns a full-resolution spectrum in the rest-frame!

    """
    from fastspecfit.emlines import EMLineModel
    
    if full_resolution:
        modelwave = CFit.sspwave # rest-frame

        # set the redshift equal to zero here but...
        redshift = 0.0 # fastspec['CONTINUUM_Z']

        # steal the code we need from emlines.EMLineFit.emlinemodel_bestfit
        EMLine = EMLineModel(redshift=redshift, npixpercamera=len(modelwave), log10wave=np.log10(modelwave))
        lineargs = [fastspec[linename.upper()] for linename in EMLine.param_names]

        emlinemodel = EMLine._emline_spectrum(*lineargs)

        # ...multiply by 1+z because the amplitudes were measured from the
        # redshifted spectra (so when we evaluate the models in
        # EMLine._emline_spectrum, the lines need to be brighter by
        # 1+z). Shoould really check this by re-measuring...

        # oh, wait, we do that below, ugh.
        
        if emlines_only:
            continuum = np.zeros_like(emlinemodel)
        else:
            continuum, _ = CFit.SSP2data(CFit.sspflux, CFit.sspwave,
                                         redshift=redshift,
                                         AV=fastspec['CONTINUUM_AV'],
                                         vdisp=fastspec['CONTINUUM_VDISP'],
                                         coeff=fastspec['CONTINUUM_COEFF'],
                                         synthphot=False)
        
        #smooth_continuum = np.zeros_like(continuum, dtype='f4')
        modelflux = continuum + emlinemodel
        
        modelflux *= (1 + redshift) # rest frame

        #import matplotlib.pyplot as plt
        #plt.plot(modelwave, modelflux) ; plt.xlim(3200, 10000) ; plt.savefig('junk.png')
        #pdb.set_trace()

        return modelwave, modelflux

    else:
        icam = 0 # one (merged) "camera"

        # mock up a fastspec-compatible dictionary of DESI data
        data = fastspec_to_desidata(fastspec, wave, flux, ivar)

        continuum, _ = CFit.SSP2data(CFit.sspflux, CFit.sspwave,
                                     redshift=data['zredrock'], 
                                     specwave=data['wave'],
                                     specres=data['res'],
                                     cameras=data['cameras'],
                                     AV=fastspec['CONTINUUM_AV'],
                                     vdisp=fastspec['CONTINUUM_VDISP'],
                                     coeff=fastspec['CONTINUUM_COEFF'],
                                     synthphot=False)
    
        continuum = continuum[icam]
        smooth_continuum = CFit.smooth_residuals(
            continuum, data['wave'][icam], data['flux'][icam],
            data['ivar'][icam], data['linemask'][icam])

        modelwave = data['wave'][icam]

        # Adjust the emission-line model range for each object otherwise the
        # trapezoidal rebinning borks..
        minwave, maxwave = modelwave.min()-1.0, modelwave.max()+1.0
        EMFit.log10wave = np.arange(np.log10(minwave), np.log10(maxwave), EMFit.dlogwave)

        emlinemodel = EMFit.emlinemodel_bestfit(data['wave'], data['res'], fastspec)[icam]

        return modelwave, continuum, smooth_continuum, emlinemodel, data

def write_templates(outfile, wave, flux, metadata, weights=None, empca=False):
    """Write out the final templates
    
    """
    from astropy.io import fits

    nobj, _ = flux.shape

    if empca:
        hduname_wave = 'WAVELENGTH'
        hduname_weights = 'WEIGHTS'
    else:
        hduname_wave = 'WAVE'
        
    hduflux = fits.PrimaryHDU(flux.astype('f4'))
    hduflux.header['EXTNAME'] = 'FLUX'
    hdulist = [hduflux]

    if empca and weights is not None:
        hduweights = fits.ImageHDU(weights.astype('f4'))
        hduweights.header['EXTNAME'] = 'WEIGHTS'
        hdulist = hdulist + [hduweights]        
        
    hduwave = fits.ImageHDU(wave.astype('f8'))
    hduwave.header['EXTNAME'] = hduname_wave
    hduwave.header['BUNIT'] = 'Angstrom'
    hduwave.header['AIRORVAC'] = ('vac', 'vacuum wavelengths')
    hdulist = hdulist + [hduwave]

    hdumetadata = fits.convenience.table_to_hdu(metadata)
    hdumetadata.header['EXTNAME'] = 'METADATA'
    hdulist = hdulist + [hdumetadata]

    hx = fits.HDUList(hdulist)

    hx.writeto(outfile, overwrite=True, checksum=True)
    log.info('Writing {} templates to {}'.format(nobj, outfile))

def _fastspec_onestack(args):
    """Multiprocessing wrapper."""
    return fastspec_onestack(*args)

def fastspec_onestack(fastmeta, fastspec, flux, ivar, meta,
                      wave, CFit, EMFit, qadir, qaprefix):
    """Fit one stacked spectrum.

    """
    # mock up a fastspec-compatible dictionary of DESI data
    data = fastspec_to_desidata(fastspec, wave, flux, ivar)

    # Adjust the emission-line fitting range for each object (otherwise the
    # trapezoidal rebinning barfs with the default wavelength array).
    minwave, maxwave = data['wave'][0].min()-3.0, data['wave'][0].max()+3.0
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

def stack_onebin(flux2d, ivar2d, templatewave, normwave, cflux2d, continuumwave):
    """Build one stack."""

    templatewave1, templateflux1, templateivar1, _, templatepix = iterative_stack(
        templatewave, flux2d, ivar2d, constant_ivar=False, normwave=normwave)

    if continuumwave is None:
        return templateflux1, templateivar1, templatepix
    else:
        _, continuumflux1, _, _, templatecpix = iterative_stack(
            continuumwave, cflux2d, normwave=normwave)
        return templateflux1, templateivar1, templatepix, continuumflux1, templatecpix

def fastspec_to_desidata(fastspec, wave, flux, ivar):
    """Convenience wrapper to turn spectra in a fastspecfit-compatible dictionary of
    DESI data.

    """
    from scipy.sparse import identity
    from desispec.resolution import Resolution

    good = np.where(ivar > 0)[0]
    npix = len(good)
    if npix == 0:
        log.warning('No good pixels in the spectrum!')
        raise ValueError

    data = {}
    data['zredrock'] = fastspec['CONTINUUM_Z']
    data['cameras'] = ['all']
    data['photsys'] = 'S'
    data['wave'] = [wave[good] * (1+data['zredrock'])]
    data['flux'] = [flux[good]]
    data['ivar'] = [ivar[good]]
    data['res'] = [Resolution(identity(n=npix))] # hack!
    data['snr'] = [np.median(flux[good] * np.sqrt(ivar[good]))]
    data['good'] = good

    # line-masking is doing some odd things around MgII so skip for now

    linemask = np.ones(npix, bool)
    #for line in CFit.linetable:
    #    zwave = line['restwave'] * (1 + data['zredrock'])
    #    if line['isbroad']:
    #        sigma = CFit.linemask_sigma_broad
    #    else:
    #        sigma = CFit.linemask_sigma
    #    I = np.where((data['wave'][0] >= (zwave - 1.5*sigma * zwave / C_LIGHT)) *
    #                 (data['wave'][0] <= (zwave + 1.5*sigma * zwave / C_LIGHT)))[0]
    #    if len(I) > 0:
    #        linemask[I] = False  # False = affected by line
    data['linemask'] = [linemask]
    
    return data
    
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
    if len(good) == 0:
        log.warning('No good pixels!')
        raise ValueError
    flux = np.sum(_ivar2d[:, good] * flux2d[:, good], axis=0) / ivar[good]
    
    #pos = np.where(flux > 0)[0]
    #good = good[pos]
    #flux = flux[pos]
    #ivar = ivar[pos]

    return wave[good], flux, ivar, nperpix[good], good

def iterative_stack(wave, flux2d, ivar2d=None, constant_ivar=False, maxdiff=0.01, 
                    maxiter=500, normwave=4500, smooth=None, debug=False, verbose=False):
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
            if normflux <= 0:
                log.warning('Normalization flux is negative or zero!')
                normflux = 1.0
                #pdb.set_trace()
                #raise ValueError

            templateflux /= normflux
            templateivar *= normflux**2
            
            break

    return templatewave, templateflux, templateivar, nperpix, goodpix

def stack_in_bins(sample, data, templatewave, mp=1, normwave=4500.0,
                  continuumwave=None, stackfile=None):
    """Stack spectra in bins given the output of a spectra_in_bins run.

    See bin/desi-templates for how to generate the required input.

    """
    npix, ntemplate = len(templatewave), len(sample)
    templateflux = np.zeros((ntemplate, npix), dtype='f4')
    templateivar = np.zeros((ntemplate, npix), dtype='f4')

    if continuumwave is not None:
        continuumflux = np.zeros((ntemplate, len(continuumwave)), dtype='f4')
    else:
        continuumflux = None

    # parallelize the stack-building step
    mpargs = []
    for samp in sample:
        ibin = samp['IBIN'].astype(str)
        if continuumwave is not None:
            mpargs.append([data[ibin]['flux'], data[ibin]['ivar'], templatewave,
                           normwave, data[ibin]['cflux'], continuumwave])
        else:
            mpargs.append([data[ibin]['flux'], data[ibin]['ivar'], templatewave,
                           normwave, None, None])

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

    stackfile is the output of stack_in_bins.

    """
    import fitsio
    from fastspecfit.continuum import ContinuumFit
    from fastspecfit.emlines import EMLineFit

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
    fastspec.rename_column('CONTINUUM_SNR_B', 'CONTINUUM_SNR_ALL') # new column!
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

def remove_undetected_lines(fastspec, linetable=None, snrmin=3.0, oiidoublet=0.73,
                            siidoublet=1.3, niidoublet=2.936, oiiidoublet=2.875,
                            fastmeta=None, devshift=False):
    """replace weak or undetected emission lines with their upper limits

    devshift - remove individual velocity shifts
    oiidoublet - [OII] 3726/3729
    siidoublet - [SII] 6716/6731
    niidoublet - [NII] 6584/6548
    oiiidoublet - [OIII] 5007/4959
    
    """
    if linetable is None:
        from fastspecfit.emlines import read_emlines
        linetable = read_emlines()
    
    linenames = linetable['name']
    redshift = fastspec['CONTINUUM_Z']
    
    for linename in linenames:
        # Fix a quick bug that affected Denali
        # https://github.com/desihub/fastspecfit/issues/21
        fix = np.where((fastspec['{}_AMP'.format(linename.upper())] != 0) *
                       (fastspec['{}_AMP_IVAR'.format(linename.upper())] == 0))[0]
        if len(fix) > 0:
            fastspec['{}_AMP'.format(linename.upper())][fix] = 0.0 

        amp = fastspec['{}_AMP'.format(linename.upper())].data
        amp_ivar = fastspec['{}_AMP_IVAR'.format(linename.upper())].data
        cont_ivar = fastspec['{}_CONT_IVAR'.format(linename.upper())].data

        if devshift:
            fastspec['{}_VSHIFT'.format(linename.upper())] = 0.0

        # we deal with the tied doublets below
        if 'oiii_4959' in linename or 'nii_6548' in linename: 
            #fix = np.where((fastspec['OIII_4959_AMP_IVAR'] > 0) * (fastspec['OIII_5007_AMP_IVAR'] == 0))[0]
            #if len(fix) > 0:
            #    fastspec['OIII_4959_AMP'][fix] = 0.0
            #    fastspec['OIII_4959_AMP_IVAR'][fix] = 0.0
            continue

        snr = amp * np.sqrt(amp_ivar)
        losnr = np.where((snr < snrmin) * (amp_ivar > 0))[0]
        #losnr = np.where(np.logical_or((amp_ivar > 0) * (snr < snrmin), (amp_ivar == 0) * (cont_ivar > 0)))[0]

        if len(losnr) > 0:
            # optionally set a minimum line-amplitude
            if False:
                csig = 1 / np.sqrt(cont_ivar[losnr])
                if np.any(csig) < 0:
                    raise ValueError
                fastspec['{}_AMP'.format(linename.upper())][losnr] = snrmin * csig
                if linename == 'oiii_5007':
                    fastspec['OIII_4959_AMP'][losnr] = fastspec['OIII_5007_AMP'][losnr].data / oiiidoublet
                if linename == 'nii_6584':
                    fastspec['NII_6548_AMP'][losnr] = fastspec['NII_6584_AMP'][losnr].data / niidoublet
            else:
                # remove the line
                fastspec['{}_AMP'.format(linename.upper())][losnr] = 0.0
                #fastspec['{}_AMP_IVAR'.format(linename.upper())][losnr] = 0.0
                if linename == 'oiii_5007':
                    fastspec['OIII_4959_AMP'][losnr] = 0.0
                if linename == 'nii_6584':
                    fastspec['NII_6548_AMP'][losnr] = 0.0

    # corner case of when the stronger doublet is on the edge of the wavelength range
    fix = np.where((fastspec['OIII_5007_AMP'] == 0) * (fastspec['OIII_4959_AMP'] != 0))[0]
    if len(fix) > 0:
        fastspec['OIII_4959_AMP'][fix] = 0.0
        fastspec['OIII_4959_AMP_IVAR'][fix] = 0.0
    fix = np.where((fastspec['NII_6584_AMP'] == 0) * (fastspec['NII_6548_AMP'] != 0))[0]
    if len(fix) > 0:
        fastspec['NII_6548_AMP'][fix] = 0.0
        fastspec['NII_6548_AMP_IVAR'][fix] = 0.0

    # double-check for negative lines
    for linename in linenames:
        neg = np.where(fastspec['{}_AMP'.format(linename.upper())] < 0)[0]
        if len(neg) > 0:
            print(neg)
            raise ValueError('Negative {}!'.format(linename))

    # restore any missing components of the [OII] or [SII] doublets 
    fix_oii3726 = np.where((fastspec['OII_3726_AMP'] == 0) * (fastspec['OII_3729_AMP'] > 0))[0]
    if len(fix_oii3726) > 0:
        fastspec['OII_3726_AMP'][fix_oii3726] = fastspec['OII_3729_AMP'][fix_oii3726] * oiidoublet

    fix_oii3729 = np.where((fastspec['OII_3729_AMP'] == 0) * (fastspec['OII_3726_AMP'] > 0))[0]
    if len(fix_oii3726) > 0:
        fastspec['OII_3729_AMP'][fix_oii3729] = fastspec['OII_3726_AMP'][fix_oii3729] / oiidoublet

    fix_sii6716 = np.where((fastspec['SII_6731_AMP'] > 0) * (fastspec['SII_6716_AMP'] == 0))[0]
    if len(fix_sii6716) > 0:
        fastspec['SII_6716_AMP'][fix_sii6716] = fastspec['SII_6731_AMP'][fix_sii6716] * siidoublet

    fix_sii6731 = np.where((fastspec['SII_6731_AMP'] == 0) * (fastspec['SII_6716_AMP'] > 0))[0]
    if len(fix_sii6731) > 0:
        fastspec['SII_6731_AMP'][fix_sii6731] = fastspec['SII_6716_AMP'][fix_sii6731] / siidoublet

    # higher-order Balmer lines require H-beta and/or H-alpha
    fix_hgamma = np.where((fastspec['HGAMMA_AMP'] > 0) * (fastspec['HBETA_AMP'] == 0) * (fastspec['HBETA_AMP_IVAR'] > 0))[0]
    if len(fix_hgamma) > 0:
        fastspec['HGAMMA_AMP'][fix_hgamma] = 0.0 # 16, 21

    fix_hdelta = np.where((fastspec['HDELTA_AMP'] > 0) * (fastspec['HBETA_AMP'] == 0) * (fastspec['HBETA_AMP_IVAR'] > 0))[0]
    if len(fix_hdelta) > 0:
        fastspec['HDELTA_AMP'][fix_hdelta] = 0.0

    fix_hepsilon = np.where((fastspec['HEPSILON_AMP'] > 0) * (fastspec['HBETA_AMP'] == 0) * (fastspec['HBETA_AMP_IVAR'] > 0))[0]
    if len(fix_hepsilon) > 0:
        fastspec['HEPSILON_AMP'][fix_hepsilon] = 0.0

    #fastspec['HGAMMA_AMP', 'HBETA_AMP'][fix_hgamma]
    #I = (fastmeta['ZOBJ'] == 0.95) * (fastmeta['MR'] == -22.75) * (fastmeta['RW1'] == -0.625)
    #fastspec[I]['OII_3726_AMP', 'OII_3729_AMP']
    #pdb.set_trace()

    # go back through and update FLUX and EW based on the new line-amplitudes
    for oneline in linetable:
        linename = oneline['name'].upper()

        amp = fastspec['{}_AMP'.format(linename.upper())].data
        amp_ivar = fastspec['{}_AMP_IVAR'.format(linename.upper())].data
        snr = amp * np.sqrt(amp_ivar)
        hisnr = np.where(snr >= snrmin)[0]
        if len(hisnr) == 0:
            continue

        linez = redshift[hisnr] + fastspec['{}_VSHIFT'.format(linename)][hisnr].data / C_LIGHT
        linezwave = oneline['restwave'] * (1 + linez)

        linesigma = fastspec['{}_SIGMA'.format(linename)][hisnr].data # [km/s]
        log10sigma = linesigma / C_LIGHT / np.log(10)     # line-width [log-10 Angstrom]
            
        # get the emission-line flux
        linesigma_ang = linezwave * linesigma / C_LIGHT # [observed-frame Angstrom]
        linenorm = np.sqrt(2.0 * np.pi) * linesigma_ang

        fastspec['{}_FLUX'.format(linename)][hisnr] = fastspec['{}_AMP'.format(linename)][hisnr].data * linenorm

        cpos = np.where(fastspec['{}_CONT'.format(linename)][hisnr] > 0.0)[0]
        if len(cpos) > 0:
            factor = (1 + redshift[hisnr][cpos]) / fastspec['{}_CONT'.format(linename)][hisnr][cpos] # --> rest frame
            fastspec['{}_EW'.format(linename)][hisnr][cpos] = fastspec['{}_FLUX'.format(linename)][hisnr][cpos] * factor   # rest frame [A]

    return fastspec

def read_stacked_fastspec(fastspecfile, read_spectra=True):
    fastmeta = Table(fitsio.read(fastspecfile, ext='METADATA'))
    fastspec = Table(fitsio.read(fastspecfile, ext='FASTSPEC'))
    if read_spectra:
        wave = fitsio.read(fastspecfile, ext='WAVE')
        flux = fitsio.read(fastspecfile, ext='FLUX')
        ivar = fitsio.read(fastspecfile, ext='IVAR')
        return wave, flux, ivar, fastmeta, fastspec
    else:
        return fastmeta, fastspec

def read_templates(templatefile, empca=False):
    wave = fitsio.read(templatefile, ext='WAVE')
    flux = fitsio.read(templatefile, ext='FLUX')
    meta = Table(fitsio.read(templatefile, ext='METADATA'))
    return wave, flux, meta

def build_templates(fastspecfile, mp=1, minwave=None, maxwave=None,
                    templatefile=None, empca=False):
    """Build the final templates.

    Called by bin/desi-templates.

    fastspecfile is the output of fastspecfit_stacks

    """
    import fitsio
    from astropy.table import join
    from fastspecfit.continuum import ContinuumFit
    from fastspecfit.emlines import EMLineFit

    if not os.path.isfile(fastspecfile):
        log.warning('fastspecfile {} not found!'.format(fastspecfile))
        raise IOError

    CFit = ContinuumFit(minwave=minwave, maxwave=maxwave)
    EMFit = EMLineFit()

    wave, flux, ivar, fastmeta, fastspec = read_stacked_fastspec(fastspecfile, read_spectra=True)

    nobj = len(fastmeta)
    log.info('Read {} fastspec model fits from {}'.format(nobj, fastspecfile))

    # Setting full_resolution here because we want to ensure the templates are
    # at full spectral resolution and in the rest-frame.
    full_resolution = True

    # Remove low-significance lines.
    fastspec = remove_undetected_lines(fastspec, EMFit.linetable, fastmeta=fastmeta)
    
    # Add lines that could not be measured from the higher-redshift stacked
    # spectra.
    if False:
        t0 = time.time()

        emlines_only = True
        mpargs = [(fastspec[iobj], wave, flux[iobj, :], ivar[iobj, :], CFit, EMFit, 
                   full_resolution, emlines_only) for iobj in np.arange(nobj)]
        if mp > 1:
            with multiprocessing.Pool(mp) as P:
                _out = P.map(_rebuild_fastspec_spectrum, mpargs)
        else:
            _out = [rebuild_fastspec_spectrum(*_mpargs) for _mpargs in mpargs]
        _out = list(zip(*_out))

        modelwave = _out[0][0] # all are identical
        emlineflux = np.vstack(_out[1])

        Mr_norm = ((fastmeta['MR'] - np.min(fastmeta['MR'])) / (np.max(fastmeta['MR']) - np.min(fastmeta['MR']))).data
        rW1_norm = ((fastmeta['RW1'] - np.min(fastmeta['RW1'])) / (np.max(fastmeta['RW1']) - np.min(fastmeta['RW1']))).data
        gi_norm = ((fastmeta['GI'] - np.min(fastmeta['GI'])) / (np.max(fastmeta['GI']) - np.min(fastmeta['GI']))).data

        maxwave = 6540.0
        zcut = 9800 / maxwave - 1
        #waverange = np.where(modelwave < maxwave)[0]

        hiz = np.where((fastmeta['ZOBJ'] > zcut) * (np.sum(emlineflux > 0, axis=1) > 0))[0] # need
        loz = np.where((fastmeta['ZOBJ'] < zcut) * (np.sum(emlineflux > 0, axis=1) > 0))[0] # have

        loz_emlineflux = emlineflux[loz, :]

        #iewhb = np.where(fastspec['HBETA_EW'] > 0)[0]
        #ewhb_norm = (fastspec['HBETA_EW'][iewhb] - np.min(fastspec['HBETA_EW'][iewhb])) /
        #    (np.max(fastspec['HBETA_EW'][iewhb]) - np.min(fastspec['HBETA_EW'][iewhb]))
        #ewhb_norm = ((fastspec['HBETA_EW'] - np.min(fastspec['HBETA_EW'][iewhb])) / (np.max(fastspec['HBETA_EW'][iewhb]) - np.min(fastspec['HBETA_EW'][iewhb]))).data
        #I = np.where((fastmeta['ZOBJ'] > 9800 / 6540-1) * (fastspec['HBETA_AMP'] > 0) * (fastspec['HALPHA_AMP'] == 0))[0]
        #J = np.where((fastmeta['ZOBJ'] < 9800 / 6540-1) * (fastspec['HBETA_AMP'] > 0))[0]

        for ihiz in hiz:
            ## Get the scale factor between this spectrum and all the other spectra.
            #good = emlineflux[loz, :]
            #hiz_emlineflux = emlineflux[ihiz, :]
            #scale = np.sum(hiz_emlineflux[np.newaxis, :] * emlineflux[loz, :], axis=1) / np.sum(emlineflux[loz, :] * emlineflux[loz, :], axis=1)

            hiz_emlineflux = emlineflux[ihiz, :]
            good = np.where(hiz_emlineflux > 0)[0]

            weight = (hiz_emlineflux > 0) * 1

            emdist = np.sqrt(np.sum((hiz_emlineflux[np.newaxis, :] - loz_emlineflux)**2, axis=1))
            #emdist = np.sum(weight * (hiz_emlineflux[np.newaxis, :] - loz_emlineflux)**2, axis=1) / np.sum(weight)
            nemdist = (emdist - np.min(emdist)) / (np.max(emdist) - np.min(emdist))
            dist = np.sqrt(nemdist**2 + (Mr_norm[loz] - Mr_norm[ihiz])**2 + (rW1_norm[loz] - rW1_norm[ihiz])**2 + (gi_norm[loz] - gi_norm[ihiz])**2)

            #dist = np.sqrt((ewhb_norm[I] - ewhb_norm[jj])**2 + (Mr_norm[I] - Mr_norm[jj])**2 + (rW1_norm[I] - rW1_norm[jj])**2)
            M = loz[np.argmin(emdist)]
            #M = loz[np.argmin(dist)]
            print(fastmeta['ZOBJ'][ihiz], fastmeta['ZOBJ'][M], fastmeta['MR'][ihiz], fastmeta['MR'][M],
                  fastmeta['RW1'][ihiz], fastmeta['RW1'][M], fastmeta['GI'][ihiz], fastmeta['GI'][M],
                  fastspec['HBETA_AMP'][ihiz], fastspec['HBETA_AMP'][M], fastspec['HALPHA_AMP'][ihiz],
                  fastspec['HALPHA_AMP'][M])#, dist)
            #for linename in ['HALPHA', 

            import matplotlib.pyplot as plt
            xlim = (3200, 7000)
            ww = np.where((modelwave > xlim[0]) * (modelwave < xlim[1]))[0]
            ylim = (0, np.max(emlineflux[M, ww]))

            #plt.clf() ; plt.plot(modelwave, emlineflux[ihiz, :]) ;  plt.xlim(3200, 7000) ; plt.savefig('junk.png')
            plt.clf() ; plt.plot(modelwave, emlineflux[ihiz, :]) ; plt.plot(modelwave, emlineflux[M, :], color='orange', alpha=0.6) ; plt.xlim(xlim) ; plt.ylim(ylim) ; plt.savefig('junk.png')
            pdb.set_trace()

        I = np.where((fastspec['HALPHA_AMP_IVAR'] == 0) * (np.sum(emlineflux > 0, axis=1) > 0))[0]

        log.info('Updating the emission-line parameters took: {:.2f} sec'.format(time.time()-t0))

        pdb.set_trace()

    # Now rebuid the full-spectrum models again with the updated emission-line
    # parameters.
    t0 = time.time()
    emlines_only = False
    mpargs = [(fastspec[iobj], wave, flux[iobj, :], ivar[iobj, :], CFit, EMFit, 
               full_resolution, emlines_only) for iobj in np.arange(nobj)]
    if mp > 1:
        with multiprocessing.Pool(mp) as P:
            _out = P.map(_rebuild_fastspec_spectrum, mpargs)
    else:
        _out = [rebuild_fastspec_spectrum(*_mpargs) for _mpargs in mpargs]
    _out = list(zip(*_out))

    # rest-frame spectra
    modelwave = _out[0][0] # all are identical
    modelflux = np.vstack(_out[1])

    log.info('Rebuilding the final model spectra took: {:.2f} sec'.format(time.time()-t0))

    #speccols = ['CONTINUUM_SNR_ALL',
    #            'BALMER_Z', 'FORBIDDEN_Z', 'BROAD_Z',
    #            'BALMER_SIGMA', 'FORBIDDEN_SIGMA', 'BROAD_SIGMA',
    #            'OII_3726_EW', 'OII_3726_EW_IVAR',
    #            'OII_3729_EW', 'OII_3729_EW_IVAR',
    #            'OIII_5007_EW', 'OIII_5007_EW_IVAR',
    #            'HBETA_EW', 'HBETA_EW_IVAR',
    #            'HALPHA_EW', 'HALPHA_EW_IVAR']
    metadata = join(fastmeta, fastspec, keys='TARGETID')
    
    if empca:
        weights = np.zeros_like(modelflux) + fastmeta['NOBJ'].data[:, np.newaxis]
    else:
        weights = None

    if templatefile:
        write_templates(templatefile, modelwave, modelflux, metadata,
                        weights=weights, empca=empca)

