#! /bin/bash

# Shell script for running mpi-fastspecfit with MPI+shifter at NERSC. Required arguments:
#   {1} stage [specfit, photfit, test]
#   {2} mp [should match the resources requested.]

# Example: build the coadds using 16 MPI tasks with 8 cores per node (and therefore 16*8/32=4 nodes)

#salloc -N 8 -C haswell -A desi -L cfs -t 02:00:00 --qos interactive --image=docker:desihub/fastspecfit:v0.0.1
#srun -n 8 -c 32 --kill-on-bad-exit=0 --no-kill shifter --module=mpich-cle6 /global/homes/i/ioannis/repos/desihub/fastspecfit/bin/mpi-fastspecfit.sh photfit 32 > photfit.log.1 2>&1 &
#srun -n 8 -c 32 --kill-on-bad-exit=0 --no-kill shifter --module=mpich-cle6 /global/homes/i/ioannis/repos/desihub/fastspecfit/bin/mpi-fastspecfit.sh specfit 32 > specfit.log.1 2>&1 &

# Grab the input arguments--
stage=$1
mp=$2

specprod=blanc

package=fastspecfit
export PATH=/global/homes/i/ioannis/repos/desihub/$package/bin:${PATH}
export PYTHONPATH=/global/homes/i/ioannis/repos/desihub/$package/py:${PYTHONPATH}
#export PYTHONPATH=/global/homes/i/ioannis/repos/git/fnnls/src/fnnls:${PYTHONPATH}
#export PYTHONPATH=/global/homes/i/ioannis/repos/git/fast-nnls:${PYTHONPATH}

export DESI_ROOT=/global/cfs/cdirs/desi
export DUST_DIR=/global/cfs/cdirs/cosmo/data/dust/v0_1
export FASTSPECFIT_DATA=/global/cfs/cdirs/desi/users/ioannis/fastspecfit
export FASTSPECFIT_TEMPLATES=$DESI_ROOT/science/gqp/templates/SSP-CKC14z

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export KMP_AFFINITY=disabled
export MPICH_GNI_FORK_MODE=FULLCOPY

if [ $stage = "test" ]; then
    time python /global/homes/i/ioannis/repos/desihub/fastspecfit/bin/mpi-fastspecfit --help
elif [ $stage = "specfit" ]; then
    time python /global/homes/i/ioannis/repos/desihub/fastspecfit/bin/mpi-fastspecfit --mp $mp --specprod $specprod --tile 80607 80608 80613
elif [ $stage = "photfit" ]; then
    time python /global/homes/i/ioannis/repos/desihub/fastspecfit/bin/mpi-fastspecfit --photfit --mp $mp --specprod $specprod --tile 80607 80608 80613
else
    echo "Unrecognized stage "$stage
fi
