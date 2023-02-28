FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive
RUN apt -y update && apt install -y apt-utils && \
    apt install -y --no-install-recommends \
    python3 python3-pip less curl wget pkg-config git build-essential make gcc file

RUN apt install -y --no-install-recommends \
    emacs \
    patch \
    libcfitsio-dev \
    libcfitsio-bin \
    libbz2-dev \
    python3-dev \
    python3-pil \
    python3-tk

RUN for x in \
    setuptools \
    wheel \
    numpy \
    scipy \
    matplotlib \
    astropy \
    fitsio \
    Django \
    django-filter \
    ; do pip3 install $x; done
#&& rm -Rf /root/.cache/pip

RUN for x in \
    uWSGI \
    ; do pip3 install $x; done

RUN apt update && apt install -y --no-install-recommends \
    postgresql-client-12
RUN pip3 install --no-cache psycopg2-binary

RUN echo "desi:desi:44591:58102:desi,,,:/:/bin/bash" >> /etc/passwd
RUN echo "desi:x:58102:" >> /etc/group

RUN mkdir /app
WORKDIR /app
COPY entry.sh /app
ENTRYPOINT ["/app/entry.sh"]

#RUN git clone https://github.com/desihub/fastspecfit.git fastspecfit && \
#    git checkout webapp
