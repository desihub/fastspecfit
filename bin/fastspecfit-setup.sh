#!/bin/bash
# Simple script to run fastspecfit in interactive / development mode. Run with
# % fastspecfit-setup.sh shifter
# % source fastspecfit-setup.sh env

if [[ $1 == "shifter" ]]; then
    # Load the desigal Docker container using shifter
    #SHIFTER=docker:desihub/fastspecfit:latest
    SHIFTER=docker:desihub/fastspecfit:2.1.1

    if [[ $NERSC_HOST == "perlmutter" ]]; then
        MPILIB=mpich
    else
        MPILIB=
    fi

    echo 'Updating and loading the shifter image '$SHIFTER
    echo 'Load the environment with: '
    echo '  source ./fastspecfit-setup.sh env'

    shifterimg pull $SHIFTER
    shifter --module=$MPILIB --image $SHIFTER bash

elif [[ $1 == "env" ]]; then
    #for package in fastspecfit; do
    #    echo Loading local check-out of $package
    #    export PATH=/global/homes/i/ioannis/code/desihub/$package/bin:$PATH
    #    export PYTHONPATH=/global/homes/i/ioannis/code/desihub/$package/py:$PYTHONPATH
    #done

    export DESI_SPECTRO_REDUX='/global/cfs/cdirs/desi/spectro/redux'
    export DUST_DIR='/global/cfs/cdirs/cosmo/data/dust/v0_1'
    export FPHOTO_DIR='/global/cfs/cdirs/desi/external/legacysurvey/dr9'
    export FTEMPLATES_DIR='/global/cfs/cdirs/desi/external/templates/fastspecfit'

    export OMP_NUM_THREADS=1
    export MKL_NUM_THREADS=1
    export KMP_AFFINITY=disabled
    export MPICH_GNI_FORK_MODE=FULLCOPY

    # Config directory nonsense
    export TMPCACHE=$(mktemp -d)
    mkdir $TMPCACHE/cache
    mkdir $TMPCACHE/config
    # astropy
    export XDG_CACHE_HOME=$TMPCACHE/cache
    export XDG_CONFIG_HOME=$TMPCACHE/config
    mkdir $XDG_CACHE_HOME/astropy
    cp -r $HOME/.astropy/cache $XDG_CACHE_HOME/astropy
    mkdir $XDG_CONFIG_HOME/astropy
    cp -r $HOME/.astropy/config $XDG_CONFIG_HOME/astropy
    # matplotlib
    export MPLCONFIGDIR=$TMPCACHE/matplotlib
    mkdir $MPLCONFIGDIR
    cp -r $HOME/.config/matplotlib $MPLCONFIGDIR
    # ipython
    export IPYTHONDIR=$TMPCACHE/ipython
    mkdir $IPYTHONDIR
    cp -r $HOME/.ipython $IPYTHONDIR  
else
    echo "Call this script with:"
    echo "  /path/to/setup/script/fastspecfit-setup.sh shifter"
    echo "  source /path/to/setup/script/fastspecfit-setup.sh env"
fi
