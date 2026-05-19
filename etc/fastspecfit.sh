#!/bin/bash -l
# Run mpi-fastspecfit on a NERSC Perlmutter CPU node allocation.
# This script is intended to be called from a Slurm job file.
#
# Usage:
#   sh fastspecfit.sh N mp specprod stage [coadd_type [survey [program]]]
#
# Positional arguments:
#   N          number of nodes
#   mp         multiprocessing workers per MPI rank (1 = pure MPI)
#   specprod   spectroscopic production name (e.g. loa, iron)
#   stage      one of: fastspec, fastphot, makeqa
#   coadd_type healpix (default), cumulative, pernight, perexp
#   survey     comma-separated survey names, or '' for all (default: '')
#   program    comma-separated program names, or '' for all (default: '')
#
# Node/rank arithmetic (Perlmutter CPU, 128 physical cores per node):
#   ntasks           = 128 * N / mp    (MPI ranks)
#   cpus_per_task    = mp * 2          (physical cores * 2 for hyperthreading)
#   workers per node = (128 / mp) * mp = 128

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/fastspecfit-env.sh"

N=${1:?'N (nodes) is required'}
mp=${2:?'mp (workers per rank) is required'}
specprod=${3:?'specprod is required'}
stage=${4:-'fastspec'}
coadd_type=${5:-'healpix'}
survey=${6:-''}
program=${7:-''}

outdir_data=${OUTDIR_DATA:-${PSCRATCH}/fastspecfit/data}

args="--outdir-data=${outdir_data} --specprod=${specprod} --mp=${mp}"

case "${stage}" in
    fastphot) args+=" --fastphot" ;;
    makeqa)   args+=" --makeqa" ;;
esac

[[ -n "${coadd_type}" ]] && args+=" --coadd-type=${coadd_type}"
[[ -n "${survey}" ]]     && args+=" --survey=${survey}"
[[ -n "${program}" ]]    && args+=" --program=${program}"

ntasks=$(( 128 * N / mp ))

if [[ ${mp} -gt 1 ]]; then
    cpus_per_task=$(( mp * 2 ))
    cpu_bind="none"
else
    cpus_per_task=2
    cpu_bind="cores"
fi

mpiscript=$(type -p mpi-fastspecfit)

srun_args="--nodes=${N} --ntasks=${ntasks} --cpus-per-task=${cpus_per_task} --cpu-bind=${cpu_bind}"

cmd="time srun ${srun_args} ${mpiscript} ${args}"
echo "${cmd}"
${cmd}
