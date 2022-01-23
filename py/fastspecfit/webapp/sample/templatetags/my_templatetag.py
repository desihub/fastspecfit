"""
This document contains our custom template tags created for use in the legacyhalos html documents. 
Each function must be registered before use, then loaded using {% load my_templatetag %} within the html document.
These functions can then be called from within the html code.
"""
import os
import astropy.io.fits
import tempfile
from django import template
from django.http import HttpResponse, HttpResponseRedirect, HttpResponseBadRequest, HttpResponseForbidden, QueryDict, StreamingHttpResponse
from astropy.table import Table
import numpy as np

register = template.Library()

#try another decorator?
@register.simple_tag
def url_replace(req, field, value):
    """
    Replace the old GET value of desired field with a new value.
    
    Args:
        req: the http request
        field: the field to replace
        value: the new value
    
    Returns:
        The updated url with the new value
    """
    dict_ = req.GET.copy()
    dict_[field] = value
    return dict_.urlencode()

@register.simple_tag
def url_replace_sort(req, new_sort):
    """
    Replace the old GET value of sort with a new value, or negate it if they are equal to sort the opposite way.
    
    Args:
        req: the http request
        new_sort: the sort value a user clicked on
        
    Returns:
        The updated url with the new sort value
    """
    dict_ = req.GET.copy()
    if 'sort' in dict_ and dict_['sort'] is not "":
        current_sort = dict_['sort']
        if current_sort == new_sort:
            dict_['sort'] = '-' + new_sort
        else:
            dict_['sort'] = new_sort
    else:
        dict_['sort'] = new_sort
    return dict_.urlencode()
    

@register.simple_tag
def url_pull(req):
    """
    Return a string describing the search criteria used.
    
    Args:
        req: the http request
        
    Returns:
        Description of search criteria
    """
    dict_ = req.GET.copy()
    search = "Search Criteria:"
    entry = False
    if "mem_match_id__gte" in dict_:
        if dict_["mem_match_id__gte"] == "":
            search += " redMaPPer ID min: 45 |"
        else:
            search += " redMaPPer ID min: " + dict_["mem_match_id__gte"] + " |"
            entry = True
    if "mem_match_id__lte" in dict_:
        if dict_["mem_match_id__lte"] == "":
            search += "redMaPPer ID max: 695620 |"
        else:
            search += " redMaPPer ID max: " + dict_["mem_match_id__lte"] + " |"
            entry = True
    if "ra__gte" in dict_:
        if dict_["ra__gte"] ==  "":
            search += "\n RA min: 0 |"
        else:
            search += "\n RA min: " + dict_["ra__gte"] + " |"
            entry = True
    if "ra__lte" in dict_:
        if dict_["ra__lte"] == "":
            search += " RA high: 360 |"
        else:
            search += " RA high: " + dict_["ra__lte"] + " |"            
            entry = True
    if "dec__gte" in dict_:
        if dict_["dec__gte"] ==  "":
            search += " Dec min: -11 |"
        else:
            search += " Dec min: " + dict_["dec__gte"] + " |"
            entry = True
    if "dec__lte" in dict_:
        if dict_["dec__lte"] == "":
            search += " Dec max: 32 |"
        else:
            search += " Dec max: " + dict_["dec__lte"] + " |"
            entry = True
    if "z__gte" in dict_:
        if dict_["z__gte"] ==  "":
            search += "\n Redshift min: 0 |"
        else:
            search += "\n Redshift min: " + dict_["z__gte"] + " |"
            entry = True
    if "z__lte" in dict_:
        if dict_["z__lte"] == "":
            search += " Redshift max: 32 |"
        else:
            search += " Redshift max: " + dict_["z__lte"] + " |"
            entry = True
    if "la__gte" in dict_:
        if dict_["la__gte"] ==  "":
            search += " Richness min: -11 |"
        else:
            search += " Richness min: " + dict_["dec__gte"] + " |"
            entry = True
    if "la__lte" in dict_:
        if dict_["la__lte"] == "":
            search += " Redshift max: 32 |"
        else:
            search += " Redshift max: " + dict_["la__lte"] + " |"
            entry = True
    if not entry:
        search = "Showing all results"
    else:
        search = search[:-1]
    return search

@register.simple_tag
def photo_pull(req, id_num, img_name):
    """
    Creates path to image based on name and redMapper id number.
    
    Args:
        req: the http request
        id_num: the redmapperID of the image galaxy
        img_name: the name of the desired image
        
    Returns: 
        Path to desired image
    """
    path = "static/data/" + id_num + "/" + id_num + "-" + img_name 
    return path    

@register.simple_tag
def viewer_link(ra, dec):
    """
    Creates a string with the viewer link for desired galaxy.
    
    Args:
        ra: the ra value to use in link
        dec: the dec value to use in link
        
    Returns: 
        Viewer link based on ra and dec values
    """
    baseurl = 'http://legacysurvey.org/viewer/'
    viewer = '{}?ra={:.6f}&dec={:.6f}&zoom=15&layer=decals-dr5'.format(baseurl, ra, dec)
    return viewer

@register.simple_tag
def skyserver_link(sdss_objid):
    """
    Creates a string with skyserver link for desired galaxy.
    
    Args:
        sdss_objid -- the sdss_objid value to use in link
        
    Returns: 
        Viewer link based on sdss_objid value
    """
    return 'http://skyserver.sdss.org/dr14/en/tools/explore/summary.aspx?id=%d' % sdss_objid
