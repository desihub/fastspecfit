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

def build_htmlhome(htmldir, fastfit, metadata, tileinfo, coadd_type='cumulative',
                   specprod='denali', htmlhome='index.html'):
    """Build the home (index.html) page.

    """
    htmldir_https = htmldir.replace(cfsroot, httpsroot)

    htmlhomefile = os.path.join(htmldir, htmlhome)
    htmlhomefile_https = os.path.join(htmldir_https, htmlhome)
    print('Building {}'.format(htmlhomefile))

    targetclasses = ('ELG', 'LRG', 'QSO', 'BGS_ANY', 'MWS_ANY')

    with open(htmlhomefile, 'w') as html:
        html.write('<html><body>\n')
        html.write('<style type="text/css">\n')
        html.write('table, td, th {padding: 5px; text-align: center; border: 1px solid black;}\n')
        html.write('p {display: inline-block;;}\n')
        html.write('</style>\n')

        html.write('<br />\n')
        html.write('<h1><a href="https://fastspecfit.readthedocs.io/en/latest/index.html">FastSpecFit</a> Visualizations</h1>\n')
        html.write('<h3><a href="https://data.desi.lbl.gov/desi/spectro/redux/{0}">{0} Data Release</a></h3>\n'.format(specprod))

        #html.write('<p style="width: 75%">\n')
        #html.write("""Write me.</p>\n""")

        html.write('<br />\n')
        html.write('<table>\n')
        html.write('<tr>\n')
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
        html.write('<th>ELG, dark</th>\n')
        html.write('<th>BGS, bright</th>\n')
        html.write('</tr>\n')

        for tinfo in tileinfo:
            tile = tinfo['TILEID']
            tiledatafile = os.path.join(htmldir, 'tiles', coadd_type, str(tile), '{}-{}.html'.format(tile, coadd_type))
            tilehtmlfile = os.path.join(htmldir_https, 'tiles', coadd_type, str(tile), '{}-{}.html'.format(tile, coadd_type))
            
            html.write('<tr>\n')
            html.write('<td><a href="{}">{}</a></td>\n'.format(tilehtmlfile, tile))
            html.write('<td>{}</td>\n'.format(tinfo['SURVEY'].upper()))
            html.write('<td>{}</td>\n'.format(tinfo['FAPRGRM'].upper()))
            html.write('<td>{}</td>\n'.format(tinfo['NEXP']))
            html.write('<td>{:.0f}</td>\n'.format(tinfo['ELG_EFFTIME_DARK']/60))
            html.write('<td>{:.0f}</td>\n'.format(tinfo['BGS_EFFTIME_BRIGHT']/60))
            html.write('</tr>\n')
        html.write('</table>\n')

        # close up shop
        html.write('<br /><br />\n')
        html.write('</html></body>\n')

    # Now the individual tile pages are more finely subdivided by target classes.
    for tinfo in tileinfo:
        tile = tinfo['TILEID']
        tiledatadir = os.path.join(htmldir, 'tiles', coadd_type, str(tile))
        if not os.path.isdir(tiledatadir):
            os.makedirs(tiledatadir, exist_ok=True)
        tiledatafile = os.path.join(tiledatadir, '{}-{}.html'.format(tile, coadd_type))
        log.info('Building {}'.format(tiledatafile))

        if tinfo['SURVEY'].upper() == 'SV1':
            from desitarget.sv1.sv1_targetmask import desi_mask#, bgs_mask, mws_mask
            desibit = 'SV1_DESI_TARGET'
            bgsbit = 'SV1_BGS_TARGET'
            mwsbit = 'SV1_MWS_TARGET'
        elif tinfo['SURVEY'].upper() == 'SV2':
            from desitarget.sv2.sv2_targetmask import desi_mask#, bgs_mask, mws_mask
            desibit = 'SV2_DESI_TARGET'
            bgsbit = 'SV2_BGS_TARGET'
            mwsbit = 'SV2_MWS_TARGET'
        else:
            NotImplementedError

        # need to account for the join suffix when we have both fastspec and
        # fastphot output
        if 'DESI_TARGET_SPEC' in metadata.colnames:
            desibit = desibit+'_SPEC'
            bgsbit = bgsbit+'_SPEC'
            mwsbit = mwsbit+'_SPEC'

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
            
            html.write('<br /><br />\n')
            html.write('<table>\n')
            html.write('<tr>\n')
            html.write('<th>Target Class</th>\n')
            html.write('<th>N</th>\n')
            html.write('</tr>\n')
                
            for targetclass in targetclasses:
                targhtmlfile = os.path.join(htmldir_https, 'tiles', coadd_type, str(tile), '{}-{}-{}.html'.format(tile, targetclass, coadd_type))
                targdatafile = os.path.join(htmldir, 'tiles', coadd_type, str(tile), '{}-{}-{}.html'.format(tile, targetclass, coadd_type))
                log.info('  Building {}'.format(targdatafile))

                targintile = np.where((tile == metadata['TILEID']) *
                                      metadata[desibit] & desi_mask.mask(targetclass) != 0)[0]
                nobj = len(targintile)

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
                    targhtml.write('<h3>Tile {} - {} (N={})</h3>\n'.format(tile, targetclass, nobj))
                    targhtml.write('<table>\n')
                    targhtml.write('<tr>\n')
                    targhtml.write('<th>TargetID</th>\n')
                    targhtml.write('<th>fastspec</th>\n')
                    targhtml.write('<th>fastphot</th>\n')
                    targhtml.write('</tr>\n')
                
                    for targetid in metadata['TARGETID'][targintile]:
                        targhtml.write('<tr>\n')
                        targhtml.write('<td>{}</td>\n'.format(targetid))

                        for prefix in ('fastspec', 'fastphot'):
                            httpfile = os.path.join(htmldir_https, 'tiles', coadd_type, str(tile), '{}-{}-{}-{}.png'.format(
                                prefix, tile, coadd_type, targetid))
                            pngfile = os.path.join(htmldir, 'tiles', coadd_type, str(tile), '{}-{}-{}-{}.png'.format(
                                prefix, tile, coadd_type, targetid))
                            if os.path.isfile(pngfile):
                                targhtml.write('<td><a href="{0}"><img src="{0}" height="auto" width="512px"></a></td>\n'.format(httpfile))
                            else:
                                targhtml.write('<td>Not Available</td>\n')
                    targhtml.write('</table>\n')
                    targhtml.write('<br /><br />\n')
                    targhtml.write('</html></body>\n')

            html.write('</table>\n')
            html.write('<br /><br />\n')
            html.write('</html></body>\n')

def _build_htmlpage_one(args):
    """Wrapper function for the multiprocessing."""
    return build_htmlpage_one(*args)

def build_htmlpage_one(ii, targetid, htmltargetid, specfile, nexttargetid, prevtargetid, nexthtmltargetid, prevhtmltargetid):
    """Build the web page for a single object.

    """
    import fitsio
    from glob import glob

    with open(htmltargetid, 'w') as html:
        html.write('<html><body>\n')
        html.write('<style type="text/css">\n')
        html.write('table, td, th {padding: 5px; text-align: center; border: 1px solid black}\n')
        html.write('</style>\n')

        #_html_group_properties(html, gal)

        html.write('<br /><br />\n')
        html.write('<a href="{}">Home</a>\n'.format(htmlhome))
        html.write('<br />\n')
        html.write('<a href="{}">Next ({})</a>\n'.format(nexthtmltargetid, nexttargetid))
        html.write('<br />\n')
        html.write('<a href="{}">Previous ({})</a>\n'.format(prevhtmltargetid, prevtargetid))
        html.write('<br />\n')

        html.write('<h1>{}</h1>'.format(targetid))
        html.write('<br /><br />\n')
        
        html.write('<table>\n')
        html.write('<tr>\n')
        #html.write('<td>{}</td>\n'.format(targetid))
        html.write('<td><a href="{}">{}</a></td>\n'.format(htmlspecfile, targetid))
        #html.write('<td><a href="{0}"><img src="{0}" height="256px"></a></td>\n'.format(htmlspecfile))
        #html.write('<td>{}</td>\n'.format(gal['INDEX']))
        html.write('</tr>\n')
        
        html.write('<br /><b><i>Last updated {}</b></i>\n'.format(js))
        html.write('<br />\n')
        html.write('</html></body>\n')

    ## Now the individual pages.
    #htmltargetids = [specfile.replace(cfsroot, httpsroot).replace('.png', '.html') for specfile in fastspecfiles]
    #
    #nexttargetids = np.roll(np.atleast_1d(targetids), -1)
    #prevtargetids = np.roll(np.atleast_1d(targetids), 1)
    #nexthtmltargetids = np.roll(np.atleast_1d(htmltargetids), -1)
    #prevhtmltargetids = np.roll(np.atleast_1d(htmltargetids), 1)
    #
    #htmlargs = []
    #for ii, (targetid, htmltargetid, specfile, nexttargetid, prevtargetid, nexthtmltargetid, prevhtmltargetid) in enumerate(
    #        zip(targetids, htmltargetids, fastspecfiles, nexttargetids, prevtargetids, nexthtmltargetids, prevhtmltargetids)):
    #    htmlargs.append([ii, targetid, htmltargetid, specfile, nexttargetid, prevtargetid, nexthtmltargetid, prevhtmltargetid])
    #
    #if args.mp > 1:
    #    import multiprocessing
    #    with multiprocessing.Pool(args.mp) as P:
    #        P.map(_build_htmlpage_one, htmlargs)
    #else:
    #    [build_htmlpage_one(*_htmlargs) for _htmlargs in htmlargs]

