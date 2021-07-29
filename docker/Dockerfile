FROM registry.services.nersc.gov/nersc/intel_cxx_fort_mpi_mkl_devel as builder

RUN apt-get -y clean && apt -y update && apt install -y apt-utils && echo 4

RUN DEBIAN_FRONTEND=noninteractive \
    apt install -y --no-install-recommends \
    make \
    git \
    file \
    patch \
    # pkg-config \
    wget \
    subversion \
    emacs \
    less \
    # gcc \
    zlib1g-dev \
    # gfortran \
    libbz2-dev \
    libgsl-dev \
    libcfitsio-dev \
    libcfitsio-bin \
    # python3 \
    # python3-dev \
    # python3-pip \
    # python3-pil \
    # python3-tk \
    # # Remove APT files
    && apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN echo 'PATH=/opt/conda/bin:$PATH' > /etc/profile.d/conda.sh

RUN wget -nv https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh && \
    bash miniconda.sh -b -p /opt/conda && \
    conda config --add channels intel && \
    conda config --set always_yes yes && \
    conda update conda && \
    conda install python=3.8 pip numpy scipy && \
    conda clean -a -y

ENV CC icc
ENV CXX icpc
ENV LDSHARED icc -shared
ENV CFLAGS -O3 -g -fPIC -std=gnu99 -pthread -x=haswell -ax=knl
ENV ARCH_FLAGS ""

ENV PYTHON python3
ENV PYTHON_CONFIG python3-config
ENV PYTHON_SCRIPT "/usr/bin/env python3"

RUN mkdir -p /src
WORKDIR /src

ENV BLAS -L${MKLROOT}/lib/intel64 -lmkl_rt -fopenmp -lpthread -limf -lsvml -ldl
ENV CPPFLAGS -I${MKLROOT}/include

RUN wget -nv https://www.mpich.org/static/downloads/3.3/mpich-3.3.tar.gz \
    && tar xvzf mpich-3.3.tar.gz \
    && cd mpich-3.3 \
    && ./configure \
    && make -j4 \
    && make install \
    && make clean \
    && cd .. && rm -Rf mpich-3.3
# defaults to installing in /usr/local

# mpi4py / mpicc can't handle -ax=knl
ENV CFLAGS -O3 -g -fPIC -std=gnu99 -pthread
RUN git clone https://github.com/mpi4py/mpi4py.git \
         && (cd mpi4py \
         && python setup.py build \
         && python setup.py install) \
         && rm -Rf mpi4py

# Pip installs -- built from source
RUN pip --no-cache-dir install fitsio

RUN /sbin/ldconfig

######## Stage 2 ########
FROM registry.services.nersc.gov/nersc/intel_cxx_fort_mpi_mkl_runtime

ENV PYTHON python3
ENV PYTHON_CONFIG python3-config
ENV PYTHON_SCRIPT "/usr/bin/env python3"
RUN echo 'PATH=/opt/conda/bin:$PATH' > /etc/profile.d/conda.sh

RUN rm /root/.profile

RUN apt -y update && apt install -y apt-utils \
 && DEBIAN_FRONTEND=noninteractive \
    apt install -y --no-install-recommends \
    make \
    git \
    file \
    pkg-config \
    wget \
    # gcc \
    zlib1g-dev \
    # gfortran \
    subversion \
    emacs \
    less \
    patch \
    libbz2-dev \
    libgsl-dev \
    libcfitsio-dev \
    libcfitsio-bin \
    imagemagick \
    # # Remove APT files
    && apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

COPY --from=builder /usr/local /usr/local
COPY --from=builder /opt/conda /opt/conda

RUN conda update --no-deps certifi
RUN for x in \
    setuptools \
    wheel \
    numpy \
    scipy \
    matplotlib \
    astropy \
    numba \
    pyyaml \
    requests \
    scikit-learn \
    healpy \
    seaborn \
    photutils \
    ipython \
    ipykernel \
    h5py \
    pytest \
    emcee \
    configobj \
    sqlalchemy \
    ; do pip3 install $x; done \
    && rm -Rf /root/.cache/pip

WORKDIR /src

# python = python3
RUN ln -s $(which python3) /usr/bin/python

RUN git clone https://github.com/dkirkby/speclite.git
RUN cd speclite \
    && python setup.py install

RUN for package in \
    desiutil \
    desimodel \
    desitarget \
    desispec \
    desisim \
    specsim; do \
      git clone https://github.com/desihub/$package \
      && cd $package \
      && python setup.py install; \
    done && echo 3

RUN git clone https://github.com/desihub/redrock
RUN cd redrock \
    && git clone https://github.com/desihub/redrock-templates py/redrock/templates \
    && python setup.py install \
    && echo 2

ENV DESIMODEL /src/desimodel

RUN mkdir $DESIMODEL \
    && cd $DESIMODEL \
    && svn export https://desi.lbl.gov/svn/code/desimodel/trunk/data \
    && echo 2

WORKDIR /src
RUN git clone https://github.com/desihub/fastspecfit.git fastspecfit \
    && cd fastspecfit && git describe --always

# update 
RUN cd /src/fastspecfit && git pull && git describe --always && echo 2

ENV PYTHONPATH /src/fastspecfit/py:$PYTHONPATH
ENV PATH /src/fastspecfit/bin:$PATH

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
  && echo 1

# Curses upon you, astropy
RUN wget -O /opt/conda/lib/python3.8/site-packages/astropy/utils/iers/data/Leap_Second.dat \
    https://raw.githubusercontent.com/astropy/astropy/master/astropy/utils/iers/data/Leap_Second.dat \
    && rm $HOME/.wget-hsts

LABEL Maintainer="John Moustakas jmoustakas@siena.edu"
