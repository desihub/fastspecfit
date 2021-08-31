"""
fastspecfit.html
================

Code to support building HTML visualizations.

"""
import pdb # for debugging
import os
import numpy as np

from desiutil.log import get_logger
log = get_logger()

cfsroot = '/global/cfs/cdirs/'
httpsroot = 'https://data.desi.lbl.gov/'

def _build_htmlpage_one(args):
    """Wrapper function for the multiprocessing."""
    return build_htmlpage_one(*args)

def build_htmlpage_one(htmlhome, meta, httpfile, pngfile, targiddatafile, nexttargetid,
                       prevtargetid, nexttargidhtmlfile, prevtargidhtmlfile):
    """Build the web page for a single object.

    """
    #log.info('  Building {}'.format(targiddatafile))
    targetid = meta['TARGETID']
    
    with open(targiddatafile, 'w') as html:
        html.write('<html><body>\n')
        html.write('<style type="text/css">\n')
        html.write('table, td, th {padding: 5px; text-align: center; border: 1px solid black}\n')
        html.write('</style>\n')

        #_properties(html, targetid)

        html.write('<br /><br />\n')
        html.write('<a href="{}">Home</a>\n'.format(htmlhome))
        html.write('<br />\n')
        html.write('<a href="{}">Next ({})</a>\n'.format(nexttargidhtmlfile, nexttargetid))
        html.write('<br />\n')
        html.write('<a href="{}">Previous ({})</a>\n'.format(prevtargidhtmlfile, prevtargetid))
        html.write('<br />\n')

        html.write('<h2>TargetID {}</h2>'.format(targetid))
        #html.write('<br /><br />\n')
        
        html.write('<table>\n')
        html.write('<tr width="90%">\n')
        html.write('<td width="50%">FastSpec QA</td><td width="50%">FastPhot QA</td>\n')
        html.write('</tr>\n')

        html.write('<tr width="90%">\n')
        for prefix in ('fastspec', 'fastphot'):
            _httpfile = httpfile.replace('prefix', prefix)
            _pngfile = pngfile.replace('prefix', prefix)
            if os.path.isfile(_pngfile):
                html.write('<td><a href="{0}"><img src="{0}" height="auto" width="512px"></a></td>\n'.format(_httpfile))
            else:
                html.write('<td>Not Available</td>\n')
        html.write('</tr>\n')
                
        #html.write('<br /><b><i>Last updated {}</b></i>\n'.format(js))
        html.write('<br />\n')
        html.write('</html></body>\n')

def build_htmlhome(htmldir, fastfit, metadata, tileinfo, coadd_type='cumulative',
                   specprod='denali', htmlhome='index.html', mp=1):
    """Build the home (index.html) page.

    """
    htmldir_https = htmldir.replace(cfsroot, httpsroot)

    htmlhomefile = os.path.join(htmldir, htmlhome)
    htmlhomefile_https = os.path.join(htmldir_https, htmlhome)
    log.info('Building {}'.format(htmlhomefile))

    targetclasses = ('ELG', 'LRG', 'QSO', 'BGS_ANY', 'MWS_ANY')
    #surveys = set(tileinfo['SURVEY'])
    surveys = ['SV1', 'SV2', 'SV3', 'Main']

    tilesurveys = np.array([ss.upper() for ss in np.atleast_1d(tileinfo['SURVEY'])])

    #pdb.set_trace()
    with open(htmlhomefile, 'w') as html:
        html.write('<html><body>\n')
        html.write('<style type="text/css">\n')
        html.write('table, td, th {padding: 5px; text-align: center; border: 1px solid black;}\n')
        html.write('p {display: inline-block;;}\n')
        html.write('</style>\n')

        html.write('<br />\n')
        html.write('<h1><a href="https://fastspecfit.readthedocs.io/en/latest/index.html">FastSpecFit</a> Visualizations</h1>\n')
        html.write('<h3><a href="https://data.desi.lbl.gov/desi/spectro/redux/{0}">{0} Data Release</a></h3>\n'.format(specprod.capitalize()))

        #html.write('<p style="width: 75%">\n')
        #html.write("""Write me.</p>\n""")

        html.write('<br />\n')
        html.write('<table>\n')
        html.write('<tr>\n')
        html.write('<th></th>\n')
        html.write('<th></th>\n')
        html.write('<th></th>\n')
        html.write('<th></th>\n')
        html.write('<th></th>\n')
        html.write('<th colspan="2">Effective Time (minutes)</th>\n')
        html.write('</tr>\n')
        
        html.write('<tr>\n')
        html.write('<th>Tile</th>\n')
        html.write('<th>Survey</th>\n')
        html.write('<th>Program</th>\n')
        html.write('<th>Nexp</th>\n')
        html.write('<th>Last Night</th>\n')
        html.write('<th>ELG, dark</th>\n')
        html.write('<th>BGS, bright</th>\n')
        html.write('</tr>\n')

        for survey in surveys:
            these = np.where(tilesurveys == survey.upper())[0]
            if len(these) == 0:
                continue
            html.write('<h2>{} Survey</h2>\n'.format(survey))
            for tinfo in tileinfo[these]:
                tile = tinfo['TILEID']
                tiledatafile = os.path.join(htmldir, survey.lower(), 'tiles', coadd_type, str(tile), '{}-{}-{}.html'.format(
                    survey.lower(), tile, coadd_type))
                tilehtmlfile = os.path.join(htmldir_https, survey.lower(), 'tiles', coadd_type, str(tile), '{}-{}-{}.html'.format(
                    survey.lower(), tile, coadd_type))

                html.write('<tr>\n')
                html.write('<td><a href="{}">{}</a></td>\n'.format(tilehtmlfile, tile))
                html.write('<td>{}</td>\n'.format(tinfo['SURVEY'].upper()))
                html.write('<td>{}</td>\n'.format(tinfo['FAPRGRM'].upper()))
                html.write('<td>{}</td>\n'.format(tinfo['NEXP']))
                html.write('<td>{}</td>\n'.format(tinfo['LASTNIGHT']))
                html.write('<td>{:.0f}</td>\n'.format(tinfo['ELG_EFFTIME_DARK']/60))
                html.write('<td>{:.0f}</td>\n'.format(tinfo['BGS_EFFTIME_BRIGHT']/60))
                html.write('</tr>\n')
        html.write('</table>\n')

        # close up shop
        html.write('<br /><br />\n')
        html.write('</html></body>\n')

    # Now the individual tile pages are more finely subdivided by target classes.
    for survey in surveys:
        these = np.where(tilesurveys == survey.upper())[0]
        if len(these) == 0:
            continue

        if survey.upper() == 'SV1':
            from desitarget.sv1.sv1_targetmask import scnd_mask, desi_mask#, bgs_mask, mws_mask
            desibit = 'SV1_DESI_TARGET'
            bgsbit = 'SV1_BGS_TARGET'
            mwsbit = 'SV1_MWS_TARGET'
            scndbit = 'SV1_SCND_TARGET'
        elif survey.upper() == 'SV2':
            from desitarget.sv2.sv2_targetmask import scnd_mask, desi_mask#, bgs_mask, mws_mask
            desibit = 'SV2_DESI_TARGET'
            bgsbit = 'SV2_BGS_TARGET'
            mwsbit = 'SV2_MWS_TARGET'
            scndbit = 'SV2_SCND_TARGET'
        elif survey.upper() == 'SV3':
            from desitarget.sv3.sv3_targetmask import scnd_mask, desi_mask#, bgs_mask, mws_mask
            desibit = 'SV3_DESI_TARGET'
            bgsbit = 'SV3_BGS_TARGET'
            mwsbit = 'SV3_MWS_TARGET'
            scndbit = 'SV3_SCND_TARGET'
        elif survey.upper() == 'MAIN':
            from desitarget.targetmask import scnd_mask, desi_mask#, bgs_mask, mws_mask
            desibit = 'DESI_TARGET'
            bgsbit = 'BGS_TARGET'
            mwsbit = 'MWS_TARGET'
            scndbit = 'SCND_TARGET'
        else:
            NotImplementedError
        
        # need to account for the join suffix when we have both fastspec and
        # fastphot output
        if 'Z_SPEC' in metadata.colnames:
            zcolumn = 'Z_SPEC'
        elif 'Z_PHOT' in metadata.colnames:
            zcolumn = 'Z_PHOT'
        else:
            zcolumn = 'Z'
        
        if 'DESI_TARGET_SPEC' in metadata.colnames:
            desibit = desibit+'_SPEC'
            bgsbit = bgsbit+'_SPEC'
            mwsbit = mwsbit+'_SPEC'
                
        for tinfo in tileinfo[these]:
            tile = tinfo['TILEID']
            tiledatadir = os.path.join(htmldir, survey.lower(), 'tiles', coadd_type, str(tile))
            if not os.path.isdir(tiledatadir):
                os.makedirs(tiledatadir, exist_ok=True)
            tiledatafile = os.path.join(tiledatadir, '{}-{}-{}.html'.format(
                survey.lower(), tile, coadd_type))
            log.info('Building {}'.format(tiledatafile))

            with open(tiledatafile, 'w') as html:
                html.write('<html><body>\n')
                html.write('<style type="text/css">\n')
                html.write('table, td, th {padding: 5px; text-align: center; border: 1px solid black;}\n')
                html.write('p {width: "75%";}\n')
                html.write('</style>\n')

                html.write('<h3><a href="{}">Home</a></h3>\n'.format(htmlhomefile_https))

                html.write('<h2>Tile {}</h2>\n'.format(tile))
                html.write('<br />\n')
                html.write('<table>\n')
                html.write('<tr>\n')
                html.write('<th></th>\n')
                html.write('<th></th>\n')
                html.write('<th></th>\n')
                html.write('<th colspan="2">Effective Time (minutes)</th>\n')
                html.write('</tr>\n')

                html.write('<tr>\n')
                html.write('<th>Survey</th>\n')
                html.write('<th>Program</th>\n')
                html.write('<th>Nexp</th>\n')
                html.write('<th>ELG, dark</th>\n')
                html.write('<th>BGS, bright</th>\n')
                html.write('</tr>\n')
                html.write('<tr>\n')
                html.write('<td>{}</td>\n'.format(tinfo['SURVEY'].upper()))
                html.write('<td>{}</td>\n'.format(tinfo['FAPRGRM'].upper()))
                html.write('<td>{}</td>\n'.format(tinfo['NEXP']))
                html.write('<td>{:.0f}</td>\n'.format(tinfo['ELG_EFFTIME_DARK']/60))
                html.write('<td>{:.0f}</td>\n'.format(tinfo['BGS_EFFTIME_BRIGHT']/60))
                html.write('</tr>\n')
                html.write('</table>\n')

                #html.write('<br /><br />\n')

                html.write('<h3>Primary Targets</h3>\n')
                html.write('<table>\n')
                html.write('<tr>\n')
                html.write('<th>Target Class</th>\n')
                html.write('<th>N</th>\n')
                html.write('</tr>\n')

                for targetclass in targetclasses:
                    targintile = np.where((tile == metadata['TILEID']) *
                                          metadata[desibit] & desi_mask.mask(targetclass) != 0)[0]
                    nobj = len(targintile)
                    
                    targhtmlfile = os.path.join(htmldir_https, survey.lower(), 'tiles', coadd_type, str(tile), '{}-{}-{}-{}.html'.format(
                        survey.lower(), targetclass.lower(), tile, coadd_type))
                    targdatafile = os.path.join(htmldir, survey.lower(), 'tiles', coadd_type, str(tile), '{}-{}-{}-{}.html'.format(
                        survey.lower(), targetclass.lower(), tile, coadd_type))

                    html.write('<tr>\n')
                    html.write('<td><a href="{}">{}</a></td>\n'.format(targhtmlfile, targetclass))
                    html.write('<td>{}</td>\n'.format(nobj))
                    html.write('</tr>\n')

                    with open(targdatafile, 'w') as targhtml:
                        targhtml.write('<html><body>\n')
                        targhtml.write('<style type="text/css">\n')
                        targhtml.write('table, td, th {padding: 5px; text-align: center; border: 1px solid black;}\n')
                        targhtml.write('p {width: "75%";}\n')
                        targhtml.write('</style>\n')
                        targhtml.write('<h3><a href="{}">Home</a></h3>\n'.format(htmlhomefile_https))
                        targhtml.write('<h3>{} Survey - Tile {} - {} Target Class (N={})</h3>\n'.format(
                            survey, tile, targetclass, nobj))

                        targhtml.write('<table>\n')
                        targhtml.write('<tr>\n')
                        targhtml.write('<th>TargetID</th>\n')
                        targhtml.write('<th>Redshift</th>\n')
                        targhtml.write('</tr>\n')
                        for meta in metadata[targintile]:
                            targetid = meta['TARGETID']
                            targidhtmlfile = os.path.join(htmldir_https, survey.lower(), 'tiles', coadd_type, str(tile), targetclass.lower(),
                                                          '{}-{}-{}-{}-{}.html'.format(survey.lower(), targetclass.lower(), tile, targetid, coadd_type))
                            targiddatadir = os.path.join(htmldir, survey.lower(), 'tiles', coadd_type, str(tile), targetclass.lower())
                            targiddatafile = os.path.join(targiddatadir, '{}-{}-{}-{}-{}.html'.format(
                                survey.lower(), targetclass.lower(), tile, targetid, coadd_type))
                            if not os.path.isdir(targiddatadir):
                                os.makedirs(targiddatadir, exist_ok=True)
                            
                            targhtml.write('<tr>\n')
                            targhtml.write('<td><a href="{}">{}</a></td>\n'.format(targidhtmlfile, targetid))
                            targhtml.write('<td>{:.5f}</td>\n'.format(meta[zcolumn]))
                            targhtml.write('</tr>\n')
                        targhtml.write('</table>\n')
                        targhtml.write('<br /><br />\n')
                        targhtml.write('</html></body>\n')

                #html.write('<h3>Secondary Targets</h3>\n')
                #html.write('<table>\n')
                #html.write('<tr>\n')
                #html.write('<th>Secondary Program</th>\n')
                #html.write('<th>N</th>\n')
                #html.write('</tr>\n')

                html.write('</table>\n')
                html.write('<br /><br />\n')
                html.write('</html></body>\n')

                # Finally build the actually individual files.
                for targetclass in targetclasses:
                    targintile = np.where((tile == metadata['TILEID']) *
                                          metadata[desibit] & desi_mask.mask(targetclass) != 0)[0]
                    if len(targintile) == 0:
                        continue

                    targetids = metadata[targintile]['TARGETID'].data
                    targidhtmlfiles, targiddatafiles, httpfiles, pngfiles = [], [], [], []
                    for targetid in targetids:
                        targidhtmlfiles.append(os.path.join(htmldir_https, survey.lower(), 'tiles', coadd_type, str(tile), targetclass.lower(),
                                                            '{}-{}-{}-{}-{}.html'.format(survey.lower(), targetclass.lower(), tile, targetid, coadd_type)))
                        targiddatadir = os.path.join(htmldir, survey.lower(), 'tiles', coadd_type, str(tile), targetclass.lower())
                        targiddatafiles.append(os.path.join(targiddatadir, '{}-{}-{}-{}-{}.html'.format(
                            survey.lower(), targetclass.lower(), tile, targetid, coadd_type)))

                        httpfiles.append(os.path.join(htmldir_https, 'tiles', coadd_type, str(tile), 'prefix-{}-{}-{}.png'.format(
                            tile, coadd_type, targetid)))
                        pngfiles.append(os.path.join(htmldir, 'tiles', coadd_type, str(tile), 'prefix-{}-{}-{}.png'.format(
                            tile, coadd_type, targetid)))

                    nexttargetids = np.roll(np.atleast_1d(targetids), -1)
                    prevtargetids = np.roll(np.atleast_1d(targetids), 1)
                    nexttargidhtmlfiles = np.roll(np.atleast_1d(targidhtmlfiles), -1)
                    prevtargidhtmlfiles = np.roll(np.atleast_1d(targidhtmlfiles), 1)

                    htmlargs = []
                    for meta, targiddatafile, httpfile, pngfile, nexttargetid, prevtargetid, nexttargidhtmlfile, prevtargidhtmlfile in zip(
                            metadata, targiddatafiles, httpfiles, pngfiles, nexttargetids, prevtargetids, nexttargidhtmlfiles, prevtargidhtmlfiles):
                        htmlargs.append([htmlhomefile_https, meta, httpfile, pngfile, targiddatafile, nexttargetid,
                                         prevtargetid, nexttargidhtmlfile, prevtargidhtmlfile])

                    if mp > 1:
                        import multiprocessing
                        with multiprocessing.Pool(mp) as P:
                            P.map(_build_htmlpage_one, htmlargs)
                    else:
                        [build_htmlpage_one(*_htmlargs) for _htmlargs in htmlargs]
