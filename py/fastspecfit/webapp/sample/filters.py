#!/usr/bin/env python

"""Custom filters for the Sample model, which work by selecting Sample objects
in the database based on meeting the desired criteria.

"""
import django_filters
from fastspecfit.webapp.sample.models import Sample

class SampleFilter(django_filters.FilterSet):
    """Custom filter for the Sample model.  Filter options include greater than or
    equal to, and less than or equal to on the following fields: ra, dec,
    targetid, ...

    The filter can be used in a form (see, e.g., list.html).

    """
    #field_name is the Sample object variable
    #lookup_expr is used to get ranges (currently using greater/less than or equal to  
    survey__match = django_filters.CharFilter(field_name='survey', lookup_expr='icontains')
    faprgrm__match = django_filters.CharFilter(field_name='faprgrm', lookup_expr='icontains')

    tileid__match = django_filters.CharFilter(field_name='tileid_list', lookup_expr='icontains')
    targetid__match = django_filters.CharFilter(field_name='targetid', lookup_expr='icontains')
    healpix__match = django_filters.CharFilter(field_name='healpix', lookup_expr='icontains')

    z__gte = django_filters.NumberFilter(field_name='z', lookup_expr='gte')
    z__lte = django_filters.NumberFilter(field_name='z', lookup_expr='lte')

    class Meta:
        model = Sample
        #add variable to fields[] if looking for exact match
        fields = []

        def id(self):
            return self.targetid
