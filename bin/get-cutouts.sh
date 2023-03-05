#! /bin/bash

#salloc -N 4 -C cpu -A desi -L cfs -t 04:00:00 --qos interactive --image=dstndstn/viewer-cutouts 
#srun -n 4 -c 128 shifter /global/homes/i/ioannis/code/desihub/fastspecfit/bin/get-cutouts.sh fuji 128 > /global/cfs/cdirs/desi/spectro/fastspecfit/fuji/v2.0/logs/cutouts-fuji-01.log 2>&1 &

codedir=/global/homes/i/ioannis/code/desihub
mpiscript=$codedir/fastspecfit/bin/get-cutouts

outdir_data=/global/cfs/cdirs/desi/spectro/fastspecfit
outdir_html=/global/cfs/cdirs/desi/users/ioannis/fastspecfit

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export KMP_AFFINITY=disabled
export MPICH_GNI_FORK_MODE=FULLCOPY

specprod=$1
mp=$2
coadd_type=$3
survey=$4
program=$5

args="--outdir-data $outdir_data --outdir-html $outdir_html"

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

echo $mpiscript $args
time $mpiscript $args
