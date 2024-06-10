"""
fastspecfit.sandbox
===================

Sandbox code.

"""
import pdb
import numpy as np

from desiutil.log import get_logger
log = get_logger()


PIXKMS = 25.


def build_continuum_model(wave, templateflux, ebv, vdisp, coeff, redshift=0.):
    """Build the continuum model, given free parameters.

    """
    from scipy.ndimage import gaussian_filter1d

    power = -0.7
    AV = ebv * (wave / (1. + redshift) / 5500.)**(-power) # [mag]

    continuum = gaussian_filter1d(templateflux, sigma=vdisp/PIXKMS, axis=0).dot(coeff) * 10**(-0.4 * AV)

    return continuum


def objective(params, templateflux, wave, flux, weights, redshift):
    """Objective function.

    """
    ebv  = params[0]
    vdisp = params[1]
    coeff = params[2:]

    modelflux = build_continuum_model(wave, templateflux, ebv, vdisp, coeff, redshift=redshift)
    
    return weights * (flux - modelflux)
    


def fit_continuum(ztemplateflux, specwave, specflux, specivar, redshift=0., 
                  ebv_guess=0.1, vdisp_guess=125.):
    """Fit the stellar continuum using bounded non-linear least-squares.

    templateflux - 
    
    """
    from scipy.optimize import least_squares

    farg = (ztemplateflux, specwave, specflux, np.sqrt(specivar), redshift)

    npix, ntemplate = ztemplateflux.shape
    initial_guesses = np.hstack(([ebv_guess, vdisp_guess], np.ones(ntemplate)))
    bounds = [(0., 3.), (50., 500.), ] + [(0., np.inf)] * ntemplate

    fit_info = least_squares(objective, initial_guesses, args=farg, bounds=tuple(zip(*bounds)))
    bestparams = fit_info.x

    ebv, vdisp = bestparams[:2]
    coeff = bestparams[2:]
    print(ebv, vdisp, coeff)

    continuum = build_continuum_model(specwave, ztemplateflux, ebv, vdisp, coeff, redshift=redshift)

    import matplotlib.pyplot as plt
    plt.clf()
    plt.plot(specwave, specflux)
    plt.plot(specwave, continuum, color='red')
    plt.ylim(0, 5)
    plt.savefig('junk.png')

    pdb.set_trace()

    

    


