FROM ubuntu:20.04

RUN apt -y update && apt install -y apt-utils && echo 1

RUN DEBIAN_FRONTEND=noninteractive \
    apt install -y --no-install-recommends \
    make \
    gcc \
    gfortran \
    git \
    subversion \
    file \
    emacs \
    less \
    patch \
    pkg-config \
    wget \
    libcfitsio-dev \
    libcfitsio-bin \
    libbz2-dev \
    python3 \
    python3-dev \
    python3-pip \
    python3-pil \
    python3-tk \
    # # Remove APT files
    && apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Python related stuff
RUN echo "../site-packages" > /usr/local/lib/python3.8/dist-packages/site-packages.pth

# Pip installs
RUN for x in \
    setuptools \
    wheel \
    numpy \
    scipy \
    matplotlib \
    astropy \
    fitsio \
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
    corner \
    ppxf \
    fnnls \
    ; do pip3 install $x; done \
    && rm -Rf /root/.cache/pip

# python = python3
RUN ln -s /usr/bin/python3 /usr/bin/python

RUN mkdir -p /src
WORKDIR /src

# install FSPS
RUN git clone https://github.com/cconroy20/fsps
RUN cd fsps/src \
    && F90FLAGS="-fPIC" make

RUN git clone https://github.com/dkirkby/speclite.git
RUN cd speclite \
    && python setup.py install

ENV SPS_HOME /src/fsps

RUN git clone https://github.com/dfm/python-fsps.git \
    && cd python-fsps \
    && python setup.py install \
    && cd .. \
    && rm -rf python-fsps

RUN for package in \
    desiutil \
    desimodel \
    desitarget \
    desispec \
    desisim \
    prospect \
    specsim \
    redrock; do \
      git clone https://github.com/desihub/$package; \
      # cd $package; \
      # python setup.py develop;
    done

ENV PATH /src/desiutil/bin:/src/desimodel/bin:/src/desitarget/bin:/src/desispec/bin:/src/desisim/bin:/src/prospect/bin:/src/redrock/bin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
ENV PYTHONPATH /src/desiutil/py:/src/desimodel/py:/src/desitarget/py:/src/desispec/py:/src/desisim/py:/src/prospect/py:/src/redrock/py

RUN git clone https://github.com/desihub/redrock-templates

# RUN git clone https://github.com/desihub/fiberassign \
#     && cd fiberassign \
#     && python setup.py install

ENV DESIMODEL /src/desimodel

RUN cd $DESIMODEL \
    && svn export https://desi.lbl.gov/svn/code/desimodel/trunk/data

ENV DESI_ROOT /global/cfs/cdirs/desi
ENV RR_TEMPLATE_DIR /src/redrock-templates

RUN mkdir /homedir && chmod 777 /homedir
ENV HOME /homedir

# set prompt and default shell
SHELL ["/bin/bash", "-c"]

RUN echo "export PS1='[container] \\u@\\h:\\w$ '" >> $HOME/.bashrc \
  # Create config files in $HOME
  && python -c "import astropy" \
  && chmod -R a+rwX $HOME/.astropy \
  && python -c "import matplotlib.font_manager as fm; f = fm.FontManager()" \
  && ipython -c "print('hello')"
  
ENTRYPOINT ["/bin/bash", "-c"]
CMD ["/bin/bash"]
