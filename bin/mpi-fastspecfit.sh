#! /bin/bash

# Shell script for running mpi-fastspecfit with MPI+shifter at NERSC. Required arguments:
#   {1} stage [fastspec, fastphot, test]
#   {2} mp [should match the resources requested.]

# Example: build the coadds using 16 MPI tasks with 8 cores per node (and therefore 16*8/32=4 nodes)

#salloc -N 16 -C haswell -A desi -L cfs -t 02:00:00 --qos interactive --image=docker:desihub/fastspecfit:latest
#srun -n 16 -c 32 --kill-on-bad-exit=0 --no-kill shifter --module=mpich-cle6 /global/homes/i/ioannis/repos/desihub/fastspecfit/bin/mpi-fastspecfit.sh fastphot 32 > fastphot-everest.log.1 2>&1 &
#srun -n 16 -c 32 --kill-on-bad-exit=0 --no-kill shifter --module=mpich-cle6 /global/homes/i/ioannis/repos/desihub/fastspecfit/bin/mpi-fastspecfit.sh fastspec 32 > fastspec-everest.log.1 2>&1 &
#srun -n 16 -c 32 --kill-on-bad-exit=0 --no-kill shifter --module=mpich-cle6 /global/homes/i/ioannis/repos/desihub/fastspecfit/bin/mpi-fastspecfit.sh qafastspec 32 > qafastspec-everest.log.1 2>&1 &
#srun -n 16 -c 32 --kill-on-bad-exit=0 --no-kill shifter --module=mpich-cle6 /global/homes/i/ioannis/repos/desihub/fastspecfit/bin/mpi-fastspecfit.sh qafastphot 32 > qafastphot-everest.log.1 2>&1 &

# Grab the input arguments--
stage=$1
mp=$2

specprod=everest
coadd_type=cumulative

package=fastspecfit
#export PATH=/opt/conda/bin:$PATH # nersc hack!
export PATH=/global/homes/i/ioannis/repos/desihub/$package/bin:$PATH
export PYTHONPATH=/global/homes/i/ioannis/repos/desihub/$package/py:$PYTHONPATH

#for package in desispec desitarget desiutil desimodel speclite specsim; do
for package in desiutil desispec specter; do
    export PATH=/global/homes/i/ioannis/repos/desihub/$package/bin:$PATH
    export PYTHONPATH=/global/homes/i/ioannis/repos/desihub/$package/py:$PYTHONPATH
done

export DESI_ROOT='/global/cfs/cdirs/desi'
export DUST_DIR='/global/cfs/cdirs/cosmo/data/dust/v0_1'
export FASTSPECFIT_TEMPLATES='/global/cfs/cdirs/desi/science/gqp/templates/SSP-CKC14z'

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export KMP_AFFINITY=disabled
export MPICH_GNI_FORK_MODE=FULLCOPY

if [ $stage = "test" ]; then
    time python /global/homes/i/ioannis/repos/desihub/fastspecfit/bin/mpi-fastspecfit --help
elif [ $stage = "fastspec" ]; then
    time python /global/homes/i/ioannis/repos/desihub/fastspecfit/bin/mpi-fastspecfit --mp $mp --specprod $specprod --coadd-type $coadd_type --tile 80605 80606 80613 --base-datadir /global/cfs/cdirs/desi/spectro/fastspecfit
elif [ $stage = "fastphot" ]; then
    time python /global/homes/i/ioannis/repos/desihub/fastspecfit/bin/mpi-fastspecfit --fastphot --mp $mp --specprod $specprod --coadd-type $coadd_type --tile 80605 80606 80613 --base-datadir /global/cfs/cdirs/desi/spectro/fastspecfit
elif [ $stage = "qafastspec" ]; then
    time python /global/homes/i/ioannis/repos/desihub/fastspecfit/bin/mpi-fastspecfit --makeqa --mp $mp --specprod $specprod --coadd-type $coadd_type --tile 80605 80606 80613 --base-datadir /global/cfs/cdirs/desi/spectro/fastspecfit --base-htmldir /global/cfs/cdirs/desi/users/ioannis/fastspecfit
elif [ $stage = "qafastphot" ]; then
    time python /global/homes/i/ioannis/repos/desihub/fastspecfit/bin/mpi-fastspecfit --fastphot --makeqa --mp $mp --specprod $specprod --coadd-type $coadd_type --tile 80605 80606 80613 --base-datadir /global/cfs/cdirs/desi/spectro/fastspecfit --base-htmldir /global/cfs/cdirs/desi/users/ioannis/fastspecfit
else
    echo "Unrecognized stage "$stage
fi
