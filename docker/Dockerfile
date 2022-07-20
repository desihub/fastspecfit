FROM --platform=linux/amd64 ubuntu:22.04

RUN mkdir -p /src
WORKDIR /src

RUN apt-get -y clean && apt -y update && apt install -y apt-utils && apt -y upgrade && echo 1

RUN DEBIAN_FRONTEND=noninteractive \
    apt install -y --no-install-recommends \
    mpich \
    git \
    # wget \
    # emacs \
    less \
    python3-pip \
    # python3-numpy \
    # python3-scipy \
    # libcfitsio-dev \
    # libcfitsio-bin \
    && apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# RUN conda update --no-deps certifi
RUN for x in \
    wheel \
    pytest \
    numpy \
    scipy \
    matplotlib \
    astropy \
    numba \
    scikit-learn \
    healpy \
    seaborn \
    ipython \
    ipykernel \
    h5py \
    # setuptools \
    # pyyaml \
    # requests \
    # photutils \
    # configobj \
    # sqlalchemy \
    ; do pip3 install $x; done \
    && rm -Rf /root/.cache/pip \
    && echo 3

ENV DESIMODEL /src/desimodel
ENV DESIUTIL_VERSION 3.2.5
ENV FASTSPECFIT_VERSION main

# python = python3
RUN ln -s "$(which python3)" /usr/bin/python

# Install desiutil, fastspecfit, and then the other desihub dependencies.
RUN python -m pip install git+https://github.com/desihub/desiutil.git@${DESIUTIL_VERSION}#egg=desiutil
    
RUN git clone --branch dockerhell https://github.com/desihub/fastspecfit.git \
    && cd fastspecfit \
    && python -m pip install -r requirements.txt \
    && python setup.py install \
    && echo 1

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

# Curses upon you, astropy
RUN wget -O /opt/conda/lib/python3.9/site-packages/astropy/utils/iers/data/Leap_Second.dat \
    https://raw.githubusercontent.com/astropy/astropy/master/astropy/utils/iers/data/Leap_Second.dat \
    && rm $HOME/.wget-hsts

LABEL Maintainer="John Moustakas jmoustakas@siena.edu"
