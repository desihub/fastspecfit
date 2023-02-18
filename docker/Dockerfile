FROM ubuntu:22.04

RUN mkdir -p /src
WORKDIR /src

RUN apt-get -y clean && apt -y update && apt install -y apt-utils && apt -y upgrade && echo 2

RUN DEBIAN_FRONTEND=noninteractive \
    apt install -y --no-install-recommends \
    build-essential \
    gfortran \
    wget \
    git \
    pkg-config \
    libbz2-dev \ 
    libgsl-dev \ 
    libssl-dev \
    libcfitsio-dev \
    libcfitsio-bin \
    # libhdf5 needed by h5py
    libhdf5-dev \
    # mpich \
    python3-dev \
    python3-pip \
    # python3-mpi4py \
    && apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# python = python3
RUN ln -s "$(which python3)" /usr/bin/python

ARG mpich=4.0.2
ARG mpich_prefix=mpich-$mpich

ENV FFLAGS="-fallow-argument-mismatch" \
    FCFLAGS="-fallow-argument-mismatch"

RUN wget --no-check-certificate -nv https://www.mpich.org/static/downloads/$mpich/$mpich_prefix.tar.gz \
    && tar xvzf $mpich_prefix.tar.gz \
    && cd $mpich_prefix \
    && ./configure \
    && make -j 4 \
    && make install \
    && make clean \
    && cd .. \
    && rm -rf $mpich_prefix

RUN /sbin/ldconfig

RUN for x in \
    wheel \
    setuptools \
    pytest \
    numpy \
    scipy \
    astropy \
    healpy \
    fitsio \
    numba \
    mpi4py \
    seaborn \
    matplotlib \
    ipython \
    ipykernel \
    h5py \
    ; do pip3 install $x; done \
    && rm -Rf /root/.cache/pip \
    && echo 2

ENV DESIUTIL_VERSION 3.2.6
ENV DESIMODEL_VERSION 0.18.0
ENV DESITARGET_VERSION 2.6.0
ENV DESISPEC_VERSION 0.57.0
ENV SPECLITE_VERSION v0.16
ENV FASTSPECFIT_VERSION 2.1.1

RUN pip3 install git+https://github.com/desihub/desiutil.git@${DESIUTIL_VERSION}#egg=desiutil
RUN pip3 install git+https://github.com/desihub/desimodel.git@${DESIMODEL_VERSION}#egg=desimodel
RUN pip3 install git+https://github.com/desihub/desitarget.git@${DESITARGET_VERSION}#egg=desitarget
RUN pip3 install git+https://github.com/desihub/desispec.git@${DESISPEC_VERSION}#egg=desispec
RUN pip3 install git+https://github.com/desihub/speclite.git@${SPECLITE_VERSION}#egg=speclite
RUN pip3 install git+https://github.com/desihub/fastspecfit.git@${FASTSPECFIT_VERSION}#egg=fastspecfit && echo 1.1

ENV PYTHONPATH /src/fastspecfit/py
ENV PATH /src/fastspecfit/bin:/opt/conda/bin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin

RUN mkdir /homedir && chmod 777 /homedir
ENV HOME /homedir

# set prompt and default shell
SHELL ["/bin/bash", "-c"]
ENTRYPOINT ["/bin/bash", "-c"]
CMD ["/bin/bash"]

RUN echo "export PS1='[container] \\u@\\h:\\w$ '" >> $HOME/.bashrc \
  # Create config files in $HOME
  && python -c "import astropy" \
  && python -c "import matplotlib.font_manager as fm; f = fm.FontManager()" \
  && ipython -c "print('hello')" \
  # Download astropy site locations and USNO sky model
  && python -c "from astropy.coordinates import EarthLocation; EarthLocation._get_site_registry(force_download=True)" \
  && python -c "from astropy.coordinates import EarthLocation, SkyCoord, AltAz; from astropy.time import Time; print(EarthLocation.of_site('ctio')); print(SkyCoord(180.,-45.,unit='deg').transform_to(AltAz(obstime=Time(56806.0, format='mjd'), location=EarthLocation.of_site('ctio'))))" \
  # Download astropy IERS leap-second list
  && python -c "from astropy.time import Time; Time.now()" \
  # Make astropy cache files readable!?!!
  && chmod -R a+rwX $HOME/.astropy \
  # Make ipython config files readable!?!!
  && chmod -R a+rwX $HOME/.ipython \
  && echo 2

LABEL Maintainer="John Moustakas jmoustakas@siena.edu"
