#! /bin/bash

# Shell script for running mpi-fastspecfit with MPI+shifter at NERSC. Required arguments:
#   {1} stage [fastspec, fastphot, test]
#   {2} mp [should match the resources requested.]

# Example: build the coadds using 16 MPI tasks with 8 cores per node (and therefore 16*8/32=4 nodes)

# Perlmutter
#salloc -N 4 -C cpu -A desi -L cfs -t 04:00:00 --qos interactive --image=docker:desihub/fastspecfit:2.0.0
#srun -n 4 -c 128 --kill-on-bad-exit=0 --no-kill shifter --module=mpich /global/homes/i/ioannis/code/desihub/fastspecfit/bin/mpi-fastspecfit.sh fastspec fuji 128 healpix sv1 > /global/cfs/cdirs/desi/spectro/fastspecfit/fuji/logs/fastspec-fuji-sv1.log.1 2>&1 &

#srun -n 16 -c 32 --kill-on-bad-exit=0 --no-kill shifter --module=mpich /global/homes/i/ioannis/code/desihub/fastspecfit/bin/mpi-fastspecfit.sh fastspec 32 sv2 dark > /global/cfs/cdirs/desi/spectro/fastspecfit/fuji/logs/fastspec-fuji-sv2-dark.log.1 2>&1 &
#srun -n 16 -c 32 --kill-on-bad-exit=0 --no-kill shifter --module=mpich /global/homes/i/ioannis/code/desihub/fastspecfit/bin/mpi-fastspecfit.sh fastphot 32 sv2 dark > /global/cfs/cdirs/desi/spectro/fastspecfit/fuji/logs/fastphot-fuji-sv2-dark.log.1 2>&1 &
#srun -n 16 -c 32 --kill-on-bad-exit=0 --no-kill shifter --module=mpich /global/homes/i/ioannis/code/desihub/fastspecfit/bin/mpi-fastspecfit.sh qafastspec 32 > qafastspec-fuji.log.1 2>&1 &
#srun -n 16 -c 32 --kill-on-bad-exit=0 --no-kill shifter --module=mpich /global/homes/i/ioannis/code/desihub/fastspecfit/bin/mpi-fastspecfit.sh qafastphot 32 > qafastphot-fuji.log.1 2>&1 &

# Cori
#salloc -N 32 -C haswell -A desi -L cfs -t 04:00:00 --qos interactive --image=docker:desihub/fastspecfit:2.0.0
#srun -n 32 -c 32 --kill-on-bad-exit=0 --no-kill shifter --module=mpich /global/homes/i/ioannis/code/desihub/fastspecfit/bin/mpi-fastspecfit.sh fastspec fuji 32 healpix sv1 > /global/cfs/cdirs/desi/spectro/fastspecfit/fuji/logs/fastspec-fuji-sv2sv3.log.1 2>&1 &

#srun -n 16 -c 32 --kill-on-bad-exit=0 --no-kill shifter --module=mpich-cle6 /global/homes/i/ioannis/code/desihub/fastspecfit/bin/mpi-fastspecfit.sh fastspec fuji 32 - sv1 > /global/cfs/cdirs/desi/spectro/fastspecfit/fuji/logs/fastspec-fuji-sv1.log.1 2>&1 &
#srun -n 16 -c 32 --kill-on-bad-exit=0 --no-kill shifter --module=mpich-cle6 /global/homes/i/ioannis/code/desihub/fastspecfit/bin/mpi-fastspecfit.sh fastphot 32 sv2 dark > /global/cfs/cdirs/desi/spectro/fastspecfit/fuji/logs/fastphot-fuji-sv2-dark.log.1 2>&1 &
#srun -n 16 -c 32 --kill-on-bad-exit=0 --no-kill shifter --module=mpich-cle6 /global/homes/i/ioannis/code/desihub/fastspecfit/bin/mpi-fastspecfit.sh qafastspec 32 > qafastspec-fuji.log.1 2>&1 &
#srun -n 16 -c 32 --kill-on-bad-exit=0 --no-kill shifter --module=mpich-cle6 /global/homes/i/ioannis/code/desihub/fastspecfit/bin/mpi-fastspecfit.sh qafastphot 32 > qafastphot-fuji.log.1 2>&1 &

codedir=/global/homes/i/ioannis/code/desihub
#codedir=/usr/local/bin
mpiscript=$codedir/fastspecfit/bin/mpi-fastspecfit

for package in fastspecfit desi; do
    echo Loading local check-out of $package
    export PATH=$codedir/$package/bin:$PATH
    export PYTHONPATH=$codedir/$package/py:$PYTHONPATH
done

outdir_data=/global/cfs/cdirs/desi/spectro/fastspecfit
outdir_html=/global/cfs/cdirs/desi/users/ioannis/fastspecfit

export DESI_ROOT='/global/cfs/cdirs/desi'
export DUST_DIR='/global/cfs/cdirs/cosmo/data/dust/v0_1'
export FASTSPECFIT_TEMPLATES='/global/cfs/cdirs/desi/science/gqp/templates/fastspecfit'

export TMPCACHE=$(mktemp -d)
export MPLCONFIGDIR=$TMPCACHE/matplotlib
mkdir $MPLCONFIGDIR
cp -r $HOME/.config/matplotlib $MPLCONFIGDIR

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export KMP_AFFINITY=disabled
export MPICH_GNI_FORK_MODE=FULLCOPY

stage=$1
specprod=$2
mp=$3
coadd_type=$4
survey=$5
program=$6

#echo stage=$stage
#echo specprod=$specprod
#echo mp=$mp
#echo coadd_type=$coadd_type
#echo survey=$survey
#echo program=$program

# petal 8 on tiles 80613, 80606, 80607
#args="--outdir-data $outdir_data --outdir-html $outdir_html --healpix 17684,17685,17686,17687,17692,7015,7020,7021,7022,7023,7026,7032,7015,7020,7021,7022,7023,7032"

# tiles 80613, 80606, 80607
args="--outdir-data $outdir_data --outdir-html $outdir_html --healpix 17680,17681,17682,17683,17689,17692,17682,17683,17688,17689,17692,17688,17689,17690,17691,17692,17689,17691,17692,17694,17692,17694,17716,17692,17693,17694,17695,17692,17693,17695,17686,17687,17692,17693,17684,17685,17686,17687,17692,17681,17683,17684,17686,17692,7017,7018,7019,7020,7022,7018,7019,7022,7104,7105,7019,7022,7104,7105,7108,7022,7105,7107,7108,7022,7023,7108,7109,7022,7023,7034,7109,7120,7022,7023,7032,7033,7034,7021,7022,7023,7032,7015,7020,7021,7022,7023,7026,7032,7017,7020,7022,7017,7018,7019,7020,7022,7018,7019,7022,7104,7105,7019,7022,7104,7105,7108,7022,7105,7107,7108,7022,7023,7108,7109,7022,7023,7034,7109,7120,7022,7023,7032,7033,7034,7021,7022,7023,7032,7015,7020,7021,7022,7023,7032,7017,7020,7022"

if [[ $stage == "fastphot" ]]; then
    args=$args" --fastphot"
fi
if [[ $stage == "makeqa" ]]; then
    args=$args" --makeqa"
fi
if [[ $specprod != " " ]] && [[ $specprod != "" ]] && [[ $specprod != "-" ]]; then
    args=$args" --specprod $specprod"
fi
if [[ $mp != " " ]] && [[ $mp != "" ]] && [[ $mp != "-" ]]; then
    args=$args" --mp $mp"
fi
if [[ $coadd_type != " " ]] && [[ $coadd_type != "" ]] && [[ $coadd_type != "-" ]]; then
    args=$args" --coadd-type $coadd_type"
else
    args=$args" --coadd-type healpix"
fi
if [[ $survey != " " ]] && [[ $survey != "" ]] && [[ $survey != "-" ]]; then
    args=$args" --survey $survey"
fi
if [[ ! -z $program ]] && [[ $program != "" ]] && [[ $program != "-" ]]; then
    args=$args" --program $program"
fi

echo python $mpiscript $args
time python $mpiscript $args
