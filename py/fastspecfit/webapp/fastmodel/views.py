#!/usr/bin/env python

"""Holds the functions that send http responses to the browser, including
rendering the html pages index.html, explore.html, and fastmodel.html, or sending a
download file.

All logic that must be done before the browser renders the html occurs here,
including sessions, serialization, querying database, applying filters, and
pagination.

"""
import os, pickle, tempfile
import numpy as np

import astropy.io.fits
from astropy.table import Table, Column

if __name__ == '__main__':
    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "fastspecfit.webapp.settings")
    import django
    django.setup()

from django.shortcuts import render
from django.core.paginator import EmptyPage, PageNotAnInteger, Paginator
from django.http import HttpResponse
from django.urls import reverse
from django.db.models import Case, When

from fastspecfit.webapp.fastmodel.filters import FastModelFilter
from fastspecfit.webapp.fastmodel.models import FastModel

def explore(req):
    """Returns the explore.html file, or renders the explore.html page after it
    applies the filter, stores result to session, and sets up pagination.

    Args:
        req: the http request

    Returns:
        File stream if user clicked download, otherwise render for explore.html

    """
    import fitsio
    from fastspecfit.webapp import settings

    # If download button was pressed return the selected subset of the FITS table.
    if req.method == 'POST':
        from fastspecfit.webapp.load import DATADIR, fastspecfile

        print('download: req.GET:', req.GET)
        query = pickle.loads(req.session['fastmodel_query'])
        print('Query:', query)
        qs = FastModel.objects.all()
        qs.query = query
        inds = qs.values_list('row_index')
        datafile = os.path.join(DATADIR, fastspecfile)
        inds = np.array(inds)
        inds = inds[:,0]
        print('Query indices:', inds.shape)
        import fitsio
        fin = fitsio.FITS(datafile)
        hdu = fin['METADATA'] # fragile!
        t = hdu.read(rows=inds)
        hdr = hdu.read_header()
        phdr = fin[0].read_header()
        hdu2 = fin['FASTSPEC'] # fragile!
        t2 = hdu2.read(rows=inds)
        hdr2 = hdu2.read_header()
        #print('Read', t)
        fits = fitsio.FITS('mem://', 'rw')
        fits.write(None, header=phdr)
        fits.write(t, header=hdr, extname='METADATA')
        fits.write(t2, header=hdr2, extname='FASTSPEC')
        rawdata = fits.read_raw()
        fits.close()
        filename = 'fastspecfit-query.fits'
        res = HttpResponse(rawdata, content_type='application/fits')
        res['Content-Disposition'] = 'attachment; filename="%s"' % filename
        return res

    # Render the page based on new filter. Automatically sort by targetid if no
    # other sort value given.
    sort = None
    if "sort" in req.GET:
        sort = req.GET.get('sort')

    queryset = None

    if req.GET.get('catalog', ''):
        dirnm = settings.USER_QUERY_DIR
        catname = os.path.join(dirnm, req.GET.get('catalog')+'.fits')
        F = fitsio.FITS(catname)
        cols = F[1].get_colnames()
        upper = True
        reqcols = ['SURVEY', 'PROGRAM', 'TARGETID', 'HEALPIX']
        if not np.all(np.isin(reqcols, cols)):
            upper = False
            reqcols = ['survey', 'program', 'targetid', 'healpix']
            if not np.all(np.isin(reqcols, cols)):
                raise ValueError('One or more required columns are missing from uploaded catalog.')

        T = Table(F[1].read(columns=reqcols, upper=True))
        target_name = ['{}-{}-{}-{}'.format(survey.lower(), program.lower(), healpix, targetid)
                       for survey, program, healpix, targetid in zip(T['SURVEY'], T['PROGRAM'], T['HEALPIX'], T['TARGETID'])]
        # query and preserve order
        # https://stackoverflow.com/questions/4916851/django-get-a-queryset-from-array-of-ids-in-specific-order
        inorder = Case(*[When(target_name__iexact=targname, then=pos) for pos, targname in enumerate(target_name)])
        queryset = FastModel.objects.filter(target_name__in=target_name).order_by(inorder)
    else:
        cone_ra  = req.GET.get('conera','')
        cone_dec = req.GET.get('conedec','')
        cone_rad = req.GET.get('coneradius','')

        # save for form default
        cone_rad_arcmin = cone_rad
        if len(cone_ra) and len(cone_dec) and len(cone_rad):
            try:
                from django.db.models import F
                cone_ra = float(cone_ra)
                cone_dec = float(cone_dec)
                cone_rad = float(cone_rad) / 60.
                dd = np.deg2rad(cone_dec)
                rr = np.deg2rad(cone_ra)
                cosd = np.cos(dd)
                x, y, z = cosd * np.cos(rr), cosd * np.sin(rr), np.sin(dd)
                r2 = np.deg2rad(cone_rad)**2
                queryset = FastModel.objects.all().annotate(
                    r2=((F('ux')-x)*(F('ux')-x) +
                        (F('uy')-y)*(F('uy')-y) +
                        (F('uz')-z)*(F('uz')-z)))
                queryset = queryset.filter(r2__lt=r2)
                if sort is None:
                    sort='r2'
                #queryset = fastmodel_near_radec(cone_ra, cone_dec, cone_rad).order_by(sort)
            except ValueError:
                pass

    if queryset is None:
        queryset = FastModel.objects.all()

    if req.GET.get('catalog', None) is None:
        if sort is None:
            sort = 'targetid'
        queryset = queryset.order_by(sort)

    #apply filter to FastModel model, then store in queryset.
    fastmodel_filter = FastModelFilter(req.GET, queryset)
    fastmodel_filtered = fastmodel_filter.qs

    #use pickle to serialize queryset (just the SQL query), and store in session
    req.session['fastmodel_query'] = pickle.dumps(fastmodel_filtered.query)

    #use django pagination functionality
    paginator = Paginator(fastmodel_filtered, 50)
    page_num = req.GET.get('page')
    page = paginator.get_page(page_num)

    # Include pagination values we will use in html page in the return
    # statement.
    for sam in page:
        print(sam)

    if req.GET.get('catalog', ''):
        return render(req, 'explore.html', {'page': page, 'paginator': paginator})
    else:
        return render(req, 'explore.html', {'page': page, 'paginator': paginator,
                                            'cone_ra':cone_ra, 'cone_dec':cone_dec,
                                            'cone_rad':cone_rad_arcmin})

def target_test(req):

    return render(req, 'target-test.html')

def target(req, target_name):

    # grab this one (unique) target
    target = FastModel.objects.all().filter(target_name=target_name)
    target.order_by('targetid')

    result_index = req.GET.get('index', '-1')
    try:
        result_index = int(result_index, 10)
    except:
        result_index = -1

    has_next = has_prev = False
    if result_index > -1:
        i_next,_ = get_next_target(req, result_index)
        i_prev,_ = get_next_target(req, result_index, direction=-1)
        has_next = i_next is not None
        has_prev = i_prev is not None

    return render(req, 'target.html', {'target': target[0],
                                       'target_name': target_name,
                                       'result_index': result_index,
                                       'has_next': has_next,
                                       'has_prev': has_prev,})

def get_next_target(req, index, qs=None, direction=1):
    # "index" is actually 1-indexed...
    index -= 1
    if qs is None:
        query = pickle.loads(req.session['fastmodel_query'])
        qs = FastModel.objects.all()
        qs.query = query
    N = qs.count()
    if index >= N or index < 0:
        return None,None
    # Find the next target.
    obj = qs[index]
    targ = obj.target_name
    while True:
        index += direction
        if index >= N or index < 0:
            return None,None
        if qs[index].target_name != targ:
            return index+1, qs[index].target_name

def target_prev(req, index):
    from django.shortcuts import redirect
    index = int(index,10)
    nextindex, nexttarget = get_next_target(req, index, direction=-1)
    if nextindex is None:
        return HttpResponse('bad index')
    return redirect(reverse(target, args=(nexttarget,)) + '?index=%i' % (nextindex))

def target_next(req, index):
    from django.shortcuts import redirect
    from django.urls import reverse
    index = int(index,10)
    nextindex, nexttarget = get_next_target(req, index)
    if nextindex is None:
        return HttpResponse('bad index')
    return redirect(reverse(target, args=(nexttarget,)) + '?index=%i' % (nextindex))

def index(req):
    """
    Renders the homepage from index.html

    Args:
        req: the http request

    Returns:
        Render for index.html

    """
    return render(req, 'index.html')

def send_file(fn, content_type, unlink=False, modsince=None, expires=3600, filename=None):
    """Creates a streaminghttpresponse to send download file to browser

    Taken from unwise.views.py.

    """
    import datetime
    from django.http import HttpResponseNotModified, StreamingHttpResponse

    '''
    modsince: If-Modified-Since header string from the client.
    '''
    st = os.stat(fn)
    f = open(fn, 'rb')
    if unlink:
        os.unlink(fn)
    # file was last modified.
    lastmod = datetime.datetime.fromtimestamp(st.st_mtime)

    if modsince:
        #print('If-modified-since:', modsince #Sat, 22 Nov 2014 01:12:39 GMT)
        ifmod = datetime.datetime.strptime(modsince, '%a, %d %b %Y %H:%M:%S %Z')
        #print('Parsed:', ifmod)
        #print('Last mod:', lastmod)
        dt = (lastmod - ifmod).total_seconds()
        if dt < 1:
            return HttpResponseNotModified()

    res = StreamingHttpResponse(f, content_type=content_type)
    # res['Cache-Control'] = 'public, max-age=31536000'
    res['Content-Length'] = st.st_size
    if filename is not None:
        res['Content-Disposition'] = 'attachment; filename="%s"' % filename
    # expires in an hour?
    now = datetime.datetime.utcnow()
    then = now + datetime.timedelta(0, expires, 0)
    timefmt = '%a, %d %b %Y %H:%M:%S GMT'
    res['Expires'] = then.strftime(timefmt)
    res['Last-Modified'] = lastmod.strftime(timefmt)
    return res

def fastmodel_near_radec(ra, dec, rad, tablename='fastmodel',
                         extra_where='', clazz=FastModel):
    #from astrometry.util.starutil import deg2distsq
    dec = np.deg2rad(dec)
    ra = np.deg2rad(ra)
    cosd = np.cos(dec)
    x,y,z = cosd * np.cos(ra), cosd * np.sin(ra), np.sin(dec)
    radius = rad + np.sqrt(2.)/2. * 2048 * 2.75 / 3600. * 1.01

    ## FIXME
    r2 = np.deg2rad(radius)**2
    #r2 = deg2distsq(radius)
    sample = clazz.objects.raw(
        ('SELECT *, ((ux-(%g))*(ux-(%g))+(uy-(%g))*(uy-(%g))+(uz-(%g))*(uz-(%g))) as r2'
         + ' FROM %s where r2 <= %g %s ORDER BY r2') %
        (x,x,y,y,z,z, tablename, r2, extra_where))

    return sample

def index(req, **kwargs):
    print('Host is', req.META.get('HTTP_HOST', None))
    if is_decaps(req):
        return decaps(req)
    if is_m33(req):
        return m33(req)
    return _index(req, **kwargs)

def upload_catalog(req):
    import tempfile
    from astropy.table import Table
    #from astrometry.util.fits import fits_table
    from django.http import HttpResponseRedirect
    from django.urls import reverse
    #from map.views import index
    from fastspecfit.webapp import settings

    if req.method != 'POST':
        return HttpResponse('POST only')
    print('Files:', req.FILES)
    cat = req.FILES['catalog']

    dirnm = settings.USER_QUERY_DIR
    if not os.path.exists(dirnm):
        try:
            os.makedirs(dirnm)
        except:
            pass
    f,tmpfn = tempfile.mkstemp(suffix='.fits', dir=dirnm)
    os.close(f)
    os.unlink(tmpfn)
    print('Saving to', tmpfn)
    with open(tmpfn, 'wb+') as destination:
        for chunk in cat.chunks():
            destination.write(chunk)
    print('Wrote', tmpfn)

    errtxt = ('<html><body>%s<p>Custom catalogs must be a binary FITS binary table with mandatory columns "SURVEY", "PROGRAM", "HEALPIX" and "TARGETID".</p></body></html>')

             #+'See <a href="https://www.legacysurvey.org/svtips/">Tips & Tricks</a> for some hints on how to produce such a catalog.</p></body></html>')

    T = None
    emsg = ''
    try:
        T = Table.read(tmpfn)
    except Exception as e:
        emsg = str(e)
    if T is None:
        try:
            # Try CSV...
            from astropy.table import Table
            T = Table.read(tmpfn)#, format='ascii')
        except Exception as e:
            emsg += '; ' + str(e)
    if T is None:
        return HttpResponse(errtxt % ('Error: '+emsg))

    # Rename and resave columns if necessary
    #if rename_cols(T):
    #    T.write_to(tmpfn)

    #cols = T.colnames
    #if not (('ra' in cols) and ('dec' in cols)):
    #    return HttpResponse(errtxt % '<p>Did not find column "RA" and "DEC" in table.</p>')

    catname = tmpfn.replace(dirnm, '').replace('.fits', '')
    if catname.startswith('/'):
        catname = catname[1:]

    #from map.views import my_reverse
    return HttpResponseRedirect(reverse(explore) + '?catalog={}'.format(catname))
    #return HttpResponseRedirect(req + '?catalog={}'.format(catname))

def main():
    from django.test import Client
    import time
    c = Client()
    t0 = time.process_time()
    wt0 = time.time()
    r = c.get('/explore.html')
    t1 = time.process_time()
    wt1 = time.time()
    print('Took', t1-t0, 'cpu', wt1-wt0, 'wall')

if __name__ == '__main__':
    # fix this
    import os, sys
    os.environ['DJANGO_SETTINGS_MODULE'] = 'fastspecfit.webapp.settings'
    import django
    django.setup()
    import logging
    lvl = logging.DEBUG
    logging.basicConfig(level=lvl, format='%(message)s', stream=sys.stdout)

    from django.test import Client
    from time import time
    c = Client()
    t0 = time()
    #r = c.get('/')
    r = c.get('/?survey__match=sv1&program__match=&tileid__match=&targetid__match=&healpix__match=&targetclass__match=#results')
    print('Took %.3f' % (time()-t0))
    f = open('debug.txt', 'wb')
    for x in r:
        f.write(x)
    f.close()

    from django.db import connection
    #print(connection.queries)
    times = []
    queries = []
    for q in connection.queries:
        queries.append(q['sql'])
        times.append(float(q['time']))
    import numpy as np
    times = np.array(times)
    I = np.argsort(times)
    for i in I:
        print('Time:', times[i], 'SQL:', queries[i])

    #main()
