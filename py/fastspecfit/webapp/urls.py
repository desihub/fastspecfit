"""URL Configuration

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/2.0/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  path('', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  path('', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.urls import include, path
    2. Add a URL to urlpatterns:  path('blog/', include('blog.urls'))

"""
from django.contrib import admin
from django.urls import re_path, path
from django.views.generic import TemplateView

import fastspecfit.webapp.fastmodel.views as fastmodel

urlpatterns = [
    re_path('^$', fastmodel.explore, name='index'),
    re_path(r'target/(?P<target_name>[^/]*)\Z', fastmodel.target, name='target'),
    re_path(r'target-prev/(\d+)$', fastmodel.target_prev, name='target-prev'),
    re_path(r'target-next/(\d+)$', fastmodel.target_next, name='target-next'),
    re_path(r'upload-catalog$', fastmodel.upload_catalog, name='upload-catalog'),
]
