#!/bin/bash

# Simple script to run fastspecfit in interactive / development mode. Run with
# % fastspecfit-setup.sh shifter
# % source fastspecfit-setup.sh env

if [ $1 = "shifter" ]; then
    # Load the desigal Docker container using shifter
    SHIFTER=docker:desihub/fastspecfit:v0.0.2

    echo 'Updating and loading the shifter image '$SHIFTER
    echo 'Load the environment with: '
    echo '  source fastspecfit-setup.sh env'

    shifterimg pull $SHIFTER
    shifter --image $SHIFTER bash
elif [ $1 = "env" ]; then
    export PATH=/path/to/local/checkout/fastspecfit/bin:${PATH}
    export PYTHONPATH=/path/to/local/checkout/fastspecfit/py:${PYTHONPATH}

    export DESI_ROOT=/global/cfs/cdirs/desi
    export DESI_SPECTRO_REDUX=/global/cfs/cdirs/desi/spectro/redux
    export FASTSPECFIT_DATA=/path/to/output
    export FASTSPECFIT_HTML=/path/to/html/output
    export FASTSPECFIT_TEMPLATES=$DESI_ROOT/science/gqp/templates/SSP-CKC14z
    export DUST_DIR=/global/cfs/cdirs/cosmo/data/dust/v0_1

    export OMP_NUM_THREADS=1
    export MKL_NUM_THREADS=1
    export KMP_AFFINITY=disabled
    export MPICH_GNI_FORK_MODE=FULLCOPY
else
    echo "Please call with:"
    echo "  fastspecfit-setup.sh shifter"
    echo "  source fastspecfit-setup.sh env"
fi
