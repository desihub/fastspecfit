#! /bin/bash

# Shell script for running mpi-fastspecfit with MPI+shifter at NERSC. Required arguments:
#   {1} stage [fastspec, fastphot, test]
#   {2} mp [should match the resources requested.]

# Example: build the coadds using 16 MPI tasks with 8 cores per node (and therefore 16*8/32=4 nodes)

#salloc -N 16 -C haswell -A desi -L cfs -t 02:00:00 --qos interactive --image=docker:desihub/fastspecfit:v0.0.2
#srun -n 16 -c 32 --kill-on-bad-exit=0 --no-kill shifter --module=mpich-cle6 /global/homes/i/ioannis/repos/desihub/fastspecfit/bin/mpi-fastspecfit.sh fastphot 32 > fastphot-cascades.log.1 2>&1 &
#srun -n 16 -c 32 --kill-on-bad-exit=0 --no-kill shifter --module=mpich-cle6 /global/homes/i/ioannis/repos/desihub/fastspecfit/bin/mpi-fastspecfit.sh fastspec 32 > fastspec-cascades.log.1 2>&1 &

# Grab the input arguments--
stage=$1
mp=$2

specprod=cascades
coadd_type=deep

package=fastspecfit
export PATH=/opt/conda/bin:$PATH # nersc hack!
export PATH=/global/homes/i/ioannis/repos/desihub/$package/bin:$PATH
export PYTHONPATH=/global/homes/i/ioannis/repos/desihub/$package/py:$PYTHONPATH

for package in desispec desitarget desiutil desimodel speclite specsim; do
    export PATH=/global/homes/i/ioannis/repos/desihub/$package/bin:$PATH
    export PYTHONPATH=/global/homes/i/ioannis/repos/desihub/$package/py:$PYTHONPATH
done

export DESI_ROOT=/global/cfs/cdirs/desi
export DESI_SPECTRO_REDUX=/global/cfs/cdirs/desi/spectro/redux
export FASTSPECFIT_DATA=/global/cfs/cdirs/desi/spectro/fastspecfit
export FASTSPECFIT_HTML=/global/cfs/cdirs/desi/users/ioannis/fastspecfit
export FASTSPECFIT_TEMPLATES=$DESI_ROOT/science/gqp/templates/SSP-CKC14z
export DUST_DIR=/global/cfs/cdirs/cosmo/data/dust/v0_1

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export KMP_AFFINITY=disabled
export MPICH_GNI_FORK_MODE=FULLCOPY

if [ $stage = "test" ]; then
    time python /global/homes/i/ioannis/repos/desihub/fastspecfit/bin/mpi-fastspecfit --help
elif [ $stage = "fastspec" ]; then
    time python /global/homes/i/ioannis/repos/desihub/fastspecfit/bin/mpi-fastspecfit --mp $mp --specprod $specprod --coadd-type $coadd_type --tile 80619 80617 80643 80644 80660 80644 80645 80646 80653 80655
#--tile 80613 80608 80606 80610 80609 80605 80609 80607 80605
elif [ $stage = "fastphot" ]; then
    time python /global/homes/i/ioannis/repos/desihub/fastspecfit/bin/mpi-fastspecfit --fastphot --mp $mp --specprod $specprod --coadd-type $coadd_type --tile 80640 80641 80642 80643 80644 80645 80646 80647 80648 80649 80650 80651 80652 80653 80654 80655 80656 80657 80658 80659 80660 80661 80662 80663 80664 80665 80666 80611 80612 80613 80614 80742 80616 80617 80618 80619 80741 80740 80624 80626 80627 80628 80629 80632 80633 80635 80636 80637 80638 80639
#    80619 80617 80643 80644 80660 80644 80645 80646 80653 80655
#--tile 80613 80608 80606 80610 80609 80605 80609 80607 80605
else
    echo "Unrecognized stage "$stage
fi
