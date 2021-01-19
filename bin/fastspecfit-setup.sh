#!/bin/bash

# Simple script to run fastspecfit in interactive / development mode. Run with
# % fastspecfit-setup.sh shifter
# % source fastspecfit-setup.sh env

if [ $1 = "shifter" ]; then
    # Load the desigal Docker container using shifter
    SHIFTER=docker:desihub/fastspecfit:v0.0.1

    echo 'Updating and loading the shifter image '$SHIFTER
    echo 'Load the environment with: '
    echo '  source fastspecfit-setup.sh env'

    shifterimg pull $SHIFTER
    shifter --image $SHIFTER bash
elif [ $1 = "env" ]; then
    export PATH=/global/homes/i/ioannis/repos/desihub/fastspecfit/bin:${PATH}
    export PYTHONPATH=/global/homes/i/ioannis/repos/desihub/fastspecfit/py:${PYTHONPATH}

    export DESI_ROOT=/global/cfs/cdirs/desi
    export DUST_DIR=/global/cfs/cdirs/cosmo/data/dust/v0_1
    export FASTSPECFIT_DATA=/path/to/output
    export FASTSPECFIT_TEMPLATES=$DESI_ROOT/science/gqp/templates/SSP-CKC14z

    export OMP_NUM_THREADS=1
    export MKL_NUM_THREADS=1
    export KMP_AFFINITY=disabled
    export MPICH_GNI_FORK_MODE=FULLCOPY
else
    echo "Please call with:"
    echo "  fastspecfit-setup.sh shifter"
    echo "  source fastspecfit-setup.sh env"
fi
