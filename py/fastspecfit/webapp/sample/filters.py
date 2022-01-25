#!/usr/bin/env python

"""Custom filters for the Sample model, which work by selecting Sample objects
in the database based on meeting the desired criteria.

"""
import django_filters
from fastspecfit.webapp.sample.models import Sample

class SampleFilter(django_filters.FilterSet):
    """Custom filter for the Sample model.  Filter options include greater than or
    equal to, and less than or equal to on the following fields: ra, dec,
    sga_id, and diameter.

    The filter can be used in a form (see, e.g., list.html).

    """
    #field_name is the Sample object variable
    #lookup_expr is used to get ranges (currently using greater/less than or equal to  
    targetid__gte = django_filters.NumberFilter(field_name='targetid', lookup_expr='gte')
    targetid__lte = django_filters.NumberFilter(field_name='targetid', lookup_expr='lte')

    #galaxy__match = django_filters.CharFilter(field_name='galaxy', lookup_expr='icontains')
    #group__match = django_filters.CharFilter(field_name='group_name', lookup_expr='icontains')
    #
    #diam__gte = django_filters.NumberFilter(field_name='d26', lookup_expr='gte')
    #diam__lte = django_filters.NumberFilter(field_name='d26', lookup_expr='lte')
    #
    #groupdiam__gte = django_filters.NumberFilter(field_name='group_diam', lookup_expr='gte')
    #groupdiam__lte = django_filters.NumberFilter(field_name='group_diam', lookup_expr='lte')

    class Meta:
        model = Sample
        #add variable to fields[] if looking for exact match
        fields = []

        def id(self):
            return self.target_id
