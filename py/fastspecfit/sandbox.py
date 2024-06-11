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


def build_continuum_model(wave, templatewave, templateflux, ebv, vdisp, coeff, redshift, dluminosity, massnorm=1e10):
    """Build the continuum model, given free parameters.

    """
    from scipy.ndimage import gaussian_filter1d
    from fastspecfit.io import FLUXNORM

    npix, ntemplate = templateflux.shape

    if redshift > 0:
        ztemplatewave = templatewave * (1. + redshift)
        #T = self.full_IGM(redshift, ztemplatewave) # ToDo: IGM attenuation
        T = FLUXNORM * self.massnorm * (10. / (1e6 * dluminosity))**2 / (1. + redshift)
        ztemplateflux = templateflux * T[:, np.newaxis]
    else:
        errmsg = 'Input redshift not defined, zero, or negative!'
        log.warning(errmsg)
        ztemplatewave = templatewave#.copy() # ???
        ztemplateflux = FLUXNORM * massnorm * templateflux
        #ztemplateflux = FLUXNORM * self.massnorm * templateflux

    # ToDo: synthesize photometry
    
    


    pdb.set_trace()

    power = -0.7
    AV = ebv * (wave / (1. + redshift) / 5500.)**power # [mag]

    continuum = gaussian_filter1d(templateflux, sigma=vdisp/PIXKMS, axis=0).dot(coeff) * 10**(-0.4 * AV)

    return continuum


def objective(params, templatewave, templateflux, wave, flux, weights, redshift, dluminosity):
    """Objective function.

    """
    ebv  = params[0]
    vdisp = params[1]
    coeff = params[2:]

    modelflux = build_continuum_model(wave, templatewave, templateflux, ebv, vdisp, coeff, 
                                      redshift, dluminosity)

    return weights * (flux - modelflux)


def fit_continuum(templatewave, templateflux, specwave, specflux, specivar, specres, 
                  camerapix=None, redshift=0., dluminosity=0.,
                  ebv_guess=0.05, vdisp_guess=125., log=None):
    """Fit the stellar continuum using bounded non-linear least-squares.

    templateflux - 
    
    """
    from scipy.optimize import least_squares

    if log is None:
        from desiutil.log import get_logger
        log = get_logger()

    npix, ntemplate = templateflux.shape
    weights = np.sqrt(specivar)

    farg = (templatewave, templateflux, specwave, specflux, weights, redshift, dluminosity)
    initial_guesses = np.hstack(([ebv_guess, vdisp_guess], 1e-3 * np.ones(ntemplate)))
    bounds = [(0., 3.), (50., 500.), ] + [(0., 1e4)] * ntemplate

    fit_info = least_squares(objective, initial_guesses, args=farg, bounds=tuple(zip(*bounds)),
                             tr_solver='lsmr', tr_options={'regularize': True}, method='trf',
                             max_nfev=5000, xtol=1e-10, verbose=2)
    bestparams = fit_info.x

    ebv, vdisp = bestparams[:2]
    coeff = bestparams[2:]
    print(ebv, vdisp, coeff)

    continuum = build_continuum_model(specwave, ztemplateflux, ebv, vdisp, coeff, redshift, dluminosity)

    import matplotlib.pyplot as plt
    plt.clf()
    plt.plot(specwave, specflux)
    plt.plot(specwave, continuum, color='red')
    plt.ylim(0, 5)
    plt.xlim(5500, 6500)
    plt.savefig('junk.png')

    pdb.set_trace()

    

    


