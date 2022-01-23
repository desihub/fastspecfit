#!/bin/bash
export PYTHONPATH=$(pwd)/SGA/py:$PYTHONPATH && \
cd /app/SGA/py/SGA/webapp
uwsgi --touch-reload wsgi.py --ini uwsgi.ini
