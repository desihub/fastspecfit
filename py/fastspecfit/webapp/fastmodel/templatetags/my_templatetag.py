"""This document contains our custom template tags created for use in the
generated HTML documents.  Each function must be registered before use, then
loaded using {% load my_templatetag %} within the HTML document.  These
functions can then be called from within the HTML code.

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
    """Replace the old GET value of desired field with a new value.

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
    """Replace the old GET value of sort with a new value, or negate it if they are
    equal to sort the opposite way.

    Args:
        req: the http request
        new_sort: the sort value a user clicked on

    Returns:
        The updated url with the new sort value

    """
    dict_ = req.GET.copy()
    if 'sort' in dict_ and dict_['sort'] != "":
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
    """Return a string describing the search criteria used.

    Args:
        req: the http request

    Returns:
        Description of search criteria

    """
    dict_ = req.GET.copy()
    search = "Search Criteria:"
    entry = False
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
    if not entry:
        search = "Showing all results"
    else:
        search = search[:-1]
    return search

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
    viewer = '{}?ra={:.6f}&dec={:.6f}&zoom=15&layer=ls-dr9'.format(baseurl, ra, dec)
    return viewer
