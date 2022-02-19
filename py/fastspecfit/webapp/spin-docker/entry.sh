#!/bin/bash
export PYTHONPATH=$(pwd)/fastspecfit/py:$PYTHONPATH && \
cd /app/fastspecfit/py/fastspecfit/webapp
uwsgi --touch-reload wsgi.py --ini uwsgi.ini
