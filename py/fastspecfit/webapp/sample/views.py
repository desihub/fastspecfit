#!/usr/bin/env python

"""Holds the functions that send http responses to the browser, including
rendering the html pages index.html, explore.html, and sample.html, or sending a
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
    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "SGA.webapp.settings")
    import django
    django.setup()

from django.shortcuts import render
from django.core.paginator import EmptyPage, PageNotAnInteger, Paginator
from django.http import HttpResponse

from SGA.webapp.sample.filters import SampleFilter
from SGA.webapp.sample.models import Sample

def explore(req):
    """Returns the explore.html file, or renders the explore.html page after it
    applies the filter, stores result to session, and sets up pagination.
    
    Args:
        req: the http request
        
    Returns: 
        File stream if user clicked download, otherwise render for explore.html

    """
    ##if download button was pressed return the selected subset of the FITS table.
    if req.method == 'POST':
        from SGA.webapp.load import DATADIR
        
        print('download: req.GET:', req.GET)
        query = pickle.loads(req.session['sample_query'])
        print('Query:', query)
        qs = Sample.objects.all()
        qs.query = query
        inds = qs.values_list('row_index')
        datafile = os.path.join(DATADIR, 'SGA-2020.fits')
        inds = np.array(inds)
        inds = inds[:,0]
        print('Query indices:', inds.shape)
        import fitsio
        fin = fitsio.FITS(datafile)
        hdu = fin['ELLIPSE']
        t = hdu.read(rows=inds)
        hdr = hdu.read_header()
        phdr = fin[0].read_header()
        hdu2 = fin['TRACTOR']
        t2 = hdu2.read(rows=inds)
        hdr2 = hdu2.read_header()
        #print('Read', t)
        fits = fitsio.FITS('mem://', 'rw')
        fits.write(None, header=phdr)
        fits.write(t, header=hdr, extname='ELLIPSE')
        fits.write(t2, header=hdr2, extname='TRACTOR')
        rawdata = fits.read_raw()
        fits.close()
        filename = 'sga-query.fits'
        res = HttpResponse(rawdata, content_type='application/fits')
        res['Content-Disposition'] = 'attachment; filename="%s"' % filename
        return res

    # Render the page based on new filter. Automatically sort by sga_id if no
    # other sort value given.
    sort = None
    if "sort" in req.GET:
        sort = req.GET.get('sort')

    queryset = None
    
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
            queryset = Sample.objects.all().annotate(
                r2=((F('ux')-x)*(F('ux')-x) +
                    (F('uy')-y)*(F('uy')-y) +
                    (F('uz')-z)*(F('uz')-z)))
            queryset = queryset.filter(r2__lt=r2)
            if sort is None:
                sort='r2'
            #queryset = sample_near_radec(cone_ra, cone_dec, cone_rad).order_by(sort)
        except ValueError:
            pass

    if queryset is None:
        queryset = Sample.objects.all()

    if sort is None:
        sort = 'sga_id'

    queryset = queryset.order_by(sort)

    #apply filter to Sample model, then store in queryset
    sample_filter = SampleFilter(req.GET, queryset)
    sample_filtered = sample_filter.qs

    #use pickle to serialize queryset (just the SQL query), and store in session
    req.session['sample_query'] = pickle.dumps(sample_filtered.query)
    
    #use django pagination functionality
    paginator = Paginator(sample_filtered, 50)
    page_num = req.GET.get('page')
    page = paginator.get_page(page_num)
    
    # Include pagination values we will use in html page in the return
    # statement.
    return render(req, 'explore.html', {'page': page, 'paginator': paginator,
                                        'cone_ra':cone_ra, 'cone_dec':cone_dec,
                                        'cone_rad':cone_rad_arcmin})
    
def group(req, group_name):
    # figure out the members of this group
    members = Sample.objects.all().filter(group_name=group_name)
    members.order_by('sga_id')
    nice_group_name = group_name.replace('_GROUP', ' Group')
    primary = [m for m in members if m.group_primary]
    primary = primary[0]

    result_index = req.GET.get('index', '-1')
    try:
        result_index = int(result_index, 10)
    except:
        result_index = -1

    has_next = has_prev = False
    if result_index > -1:
        i_next,_ = get_next_group(req, result_index)
        i_prev,_ = get_next_group(req, result_index, direction=-1)
        has_next = i_next is not None
        has_prev = i_prev is not None
    
    return render(req, 'group.html', {'group_name': group_name,
                                      'nice_group_name': nice_group_name,
                                      'primary': primary,
                                      'members': members,
                                      'result_index': result_index,
                                      'has_next': has_next,
                                      'has_prev': has_prev,})

def get_next_group(req, index, qs=None, direction=1):
    # "index" is actually 1-indexed...
    index -= 1
    if qs is None:
        query = pickle.loads(req.session['sample_query'])
        qs = Sample.objects.all()
        qs.query = query
    N = qs.count()
    if index >= N or index < 0:
        return None,None
    # Find the next group.
    gal = qs[index]
    grp = gal.group_name
    while True:
        index += direction
        if index >= N or index < 0:
            return None,None
        if qs[index].group_name != grp:
            return index+1, qs[index].group_name

def group_prev(req, index):
    from django.shortcuts import redirect
    from django.urls import reverse
    index = int(index,10)
    nextindex,nextgroup = get_next_group(req, index, direction=-1)
    if nextindex is None:
        return HttpResponse('bad index')
    return redirect(reverse(group, args=(nextgroup,)) + '?index=%i' % (nextindex))

def group_next(req, index):
    from django.shortcuts import redirect
    from django.urls import reverse
    index = int(index,10)
    nextindex,nextgroup = get_next_group(req, index)
    if nextindex is None:
        return HttpResponse('bad index')
    return redirect(reverse(group, args=(nextgroup,)) + '?index=%i' % (nextindex))

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

def sample_near_radec(ra, dec, rad, tablename='sample',
                      extra_where='', clazz=Sample):
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
    main()
