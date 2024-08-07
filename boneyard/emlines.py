        # Do a fast box-car integration to get the initial line-amplitudes and
        # line-widths...actually these initial guesses are not really working
        # but keep the code here.
        initguess = False
        if initguess:
            sigma_cont = 200.0
            init_linesigmas = []
            for pp in self.EMLineModel.param_names:
                if getattr(self.EMLineModel, pp).tied:
                    #print('Skipping tied parameter {}'.format(pp))
                    continue

                if 'amp' in pp:
                    pinfo = getattr(self.EMLineModel, pp)
                    linename = pinfo.name.replace('_amp', '')

                    iline = np.where(self.linetable['name'] == linename)[0]
                    if len(iline) != 1:
                        log.warning('No matching line found!')
                        raise ValueError

                    oneline = self.linetable[iline][0]
                    linezwave = oneline['restwave'] * (1 + redshift)
                    lineindx = np.where((emlinewave > (linezwave - 5*sigma_cont * linezwave / C_LIGHT)) *
                                        (emlinewave < (linezwave + 5.*sigma_cont * linezwave / C_LIGHT)) *
                                        #(emlineflux * np.sqrt(emlineivar) > 1)
                                        (emlineivar > 0))[0]

                    linesigma = getattr(self.EMLineModel, '{}_sigma'.format(linename)).default # initial guess
                    #print(linename, linesigma, len(lineindx))

                    lineamp = pinfo.default # backup
                    if len(lineindx) > 10:
                        linesigma_ang = linezwave * sigma_cont / C_LIGHT # [observed-frame Angstrom]
                        linenorm = np.sqrt(2.0 * np.pi) * linesigma_ang

                        lineflux = np.sum(emlineflux[lineindx])
                        lineamp = np.abs(lineflux / linenorm)

                        # estimate the velocity width from potentially strong, isolated lines; fragile!
                        if lineflux > 0 and linename in ['mgii_2803', 'oiii_5007', 'hbeta', 'halpha']:
                            linevar = np.sum(emlineflux[lineindx] * (emlinewave[lineindx] - linezwave)**2) / np.sum(emlineflux[lineindx]) / linezwave * C_LIGHT # [km/s]
                            if linevar > 0:
                                linesigma = np.sqrt(linevar)
                                init_linesigmas.append(linesigma)

                    if not pinfo.tied:# and False:
                        setattr(self.EMLineModel, pp, lineamp)
                    #print(pinfo.name, len(lineindx), lineflux, lineamp, linesigma)

            # update the initial veolocity widths
            if len(init_linesigmas) >= 3:
                init_linesigma = np.median(init_linesigmas)
                if init_linesigma > 0 and init_linesigma < 300:
                    for pp in self.EMLineModel.param_names:
                        if 'sigma' in pp:
                            setattr(self.EMLineModel, pp, init_linesigma)
