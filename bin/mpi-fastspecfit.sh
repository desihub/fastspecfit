#! /bin/bash

# Shell script for running mpi-fastspecfit with MPI+shifter at NERSC. Required arguments:
#   {1} stage [fastspec, fastphot, test]
#   {2} mp [should match the resources requested.]

# Example: build the coadds using 16 MPI tasks with 8 cores per node (and therefore 16*8/32=4 nodes)

# Perlmutter
#salloc -N 8 -C cpu -A desi -L cfs -t 04:00:00 --qos interactive --image=docker:desihub/fastspecfit:v1.0.0
#srun -n 16 -c 32 --kill-on-bad-exit=0 --no-kill shifter --module=mpich /global/homes/i/ioannis/code/desihub/fastspecfit/bin/mpi-fastspecfit.sh fastspec 32 sv2 dark > /global/cfs/cdirs/desi/spectro/fastspecfit/fuji/logs/fastspec-fuji-sv2-dark.log.1 2>&1 &
#srun -n 16 -c 32 --kill-on-bad-exit=0 --no-kill shifter --module=mpich /global/homes/i/ioannis/code/desihub/fastspecfit/bin/mpi-fastspecfit.sh fastphot 32 sv2 dark > /global/cfs/cdirs/desi/spectro/fastspecfit/fuji/logs/fastphot-fuji-sv2-dark.log.1 2>&1 &
#srun -n 16 -c 32 --kill-on-bad-exit=0 --no-kill shifter --module=mpich /global/homes/i/ioannis/code/desihub/fastspecfit/bin/mpi-fastspecfit.sh qafastspec 32 > qafastspec-fuji.log.1 2>&1 &
#srun -n 16 -c 32 --kill-on-bad-exit=0 --no-kill shifter --module=mpich /global/homes/i/ioannis/code/desihub/fastspecfit/bin/mpi-fastspecfit.sh qafastphot 32 > qafastphot-fuji.log.1 2>&1 &

# Cori
#salloc -N 8 -C haswell -A desi -L cfs -t 04:00:00 --qos interactive --image=docker:desihub/fastspecfit:v1.0.0
#srun -n 16 -c 32 --kill-on-bad-exit=0 --no-kill shifter --module=mpich-cle6 /global/homes/i/ioannis/code/desihub/fastspecfit/bin/mpi-fastspecfit.sh fastspec 32 sv2 dark > /global/cfs/cdirs/desi/spectro/fastspecfit/fuji/logs/fastspec-fuji-sv2-dark.log.1 2>&1 &
#srun -n 16 -c 32 --kill-on-bad-exit=0 --no-kill shifter --module=mpich-cle6 /global/homes/i/ioannis/code/desihub/fastspecfit/bin/mpi-fastspecfit.sh fastphot 32 sv2 dark > /global/cfs/cdirs/desi/spectro/fastspecfit/fuji/logs/fastphot-fuji-sv2-dark.log.1 2>&1 &
#srun -n 16 -c 32 --kill-on-bad-exit=0 --no-kill shifter --module=mpich-cle6 /global/homes/i/ioannis/code/desihub/fastspecfit/bin/mpi-fastspecfit.sh qafastspec 32 > qafastspec-fuji.log.1 2>&1 &
#srun -n 16 -c 32 --kill-on-bad-exit=0 --no-kill shifter --module=mpich-cle6 /global/homes/i/ioannis/code/desihub/fastspecfit/bin/mpi-fastspecfit.sh qafastphot 32 > qafastphot-fuji.log.1 2>&1 &

# Grab the input arguments--
stage=$1
mp=$2
survey=$3
program=$4
specprod=guadalupe

#coadd_type=cumulative

#for package in desispec desitarget desiutil desimodel speclite specsim; do
for package in fastspecfit; do
    echo Loading local check-out of $package
    export PATH=/global/homes/i/ioannis/code/desihub/$package/bin:$PATH
    export PYTHONPATH=/global/homes/i/ioannis/code/desihub/$package/py:$PYTHONPATH
done
#echo 'Loading local check-out of speclite'
#export PYTHONPATH=/global/homes/i/ioannis/code/desihub/speclite:$PYTHONPATH

#outdir_data=/global/cfs/cdirs/desi/users/ioannis/fastspecfit/vitiles
#outdir_html=/global/cfs/cdirs/desi/users/ioannis/fastspecfit/vitiles
outdir_data=/global/cfs/cdirs/desi/spectro/fastspecfit
outdir_html=/global/cfs/cdirs/desi/users/ioannis/fastspecfit

export DESI_ROOT='/global/cfs/cdirs/desi'
export DUST_DIR='/global/cfs/cdirs/cosmo/data/dust/v0_1'
export FASTSPECFIT_TEMPLATES='/global/cfs/cdirs/desi/science/gqp/templates/SSP-CKC14z'

# matplotlib
export TMPCACHE=$(mktemp -d)
export MPLCONFIGDIR=$TMPCACHE/matplotlib
mkdir $MPLCONFIGDIR
cp -r $HOME/.config/matplotlib $MPLCONFIGDIR

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export KMP_AFFINITY=disabled
export MPICH_GNI_FORK_MODE=FULLCOPY

if [ $stage = "test" ]; then
    time python /global/homes/i/ioannis/code/desihub/fastspecfit/bin/mpi-fastspecfit --help
elif [ $stage = "fastspec" ]; then
    #time python /global/homes/i/ioannis/code/desihub/fastspecfit/bin/mpi-fastspecfit --mp $mp --coadd-type $coadd_type --tile 80605 80606 80613 --outdir-data $outdir_data
    #time python /global/homes/i/ioannis/code/desihub/fastspecfit/bin/mpi-fastspecfit --mp $mp --coadd-type $coadd_type --tile 80648 80666 136 141 142 244 245 246 80605 80606 80613 --outdir-data $outdir_data
    #time python /global/homes/i/ioannis/code/desihub/fastspecfit/bin/mpi-fastspecfit --mp $mp --coadd-type $coadd_type --survey $survey --program $program --outdir-data $outdir_data
    time python /global/homes/i/ioannis/code/desihub/fastspecfit/bin/mpi-fastspecfit --mp $mp --survey $survey --program $program --outdir-data $outdir_data
elif [ $stage = "fastphot" ]; then
    #time python /global/homes/i/ioannis/code/desihub/fastspecfit/bin/mpi-fastspecfit --mp $mp --coadd-type $coadd_type --tile 136 141 142 244 245 246 80605 80606 80613 --outdir-data $outdir_data --fastphot --overwrite
    #time python /global/homes/i/ioannis/code/desihub/fastspecfit/bin/mpi-fastspecfit --mp $mp --coadd-type $coadd_type --tile 244 245 246 80605 80606 80613 --outdir-data $outdir_data --fastphot
    time python /global/homes/i/ioannis/code/desihub/fastspecfit/bin/mpi-fastspecfit --mp $mp --specprod $specprod --survey $survey --program $program --outdir-data $outdir_data --fastphot 
    #time python /global/homes/i/ioannis/code/desihub/fastspecfit/bin/mpi-fastspecfit --mp $mp --specprod $specprod --outdir-data $outdir_data --fastphot 
elif [ $stage = "qafastspec" ]; then
    time python /global/homes/i/ioannis/code/desihub/fastspecfit/bin/mpi-fastspecfit --makeqa --mp $mp --coadd-type $coadd_type --tile 80605 80606 80613 --outdir-data $outdir_data --outdir-html $outdir_html
elif [ $stage = "qafastphot" ]; then
    time python /global/homes/i/ioannis/code/desihub/fastspecfit/bin/mpi-fastspecfit --fastphot --makeqa --mp $mp --coadd-type $coadd_type --tile 80605 80606 80613 --outdir-data $outdir_data --outdir-html $outdir_html
else
    echo "Unrecognized stage "$stage
fi
