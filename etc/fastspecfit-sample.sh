#!/bin/bash -l
# Run mpi-fastspecfit with a --samplefile on a NERSC Perlmutter CPU node allocation.
# This script is intended to be called from a Slurm job file.
#
# Usage:
#   sh fastspecfit-sample.sh N mp specprod samplefile [stage [input_redshifts]]
#
# Positional arguments:
#   N               number of nodes
#   mp              multiprocessing workers per MPI rank (1 = pure MPI)
#   specprod        spectroscopic production name (e.g. loa, iron)
#   samplefile      full path to sample FITS file with {SURVEY,PROGRAM,HEALPIX,TARGETID}
#   stage           one of: fastspec (default), fastphot
#   input_redshifts pass --input-redshifts if 'true' (default: '')
#
# The samplefile is broadcast to all MPI ranks at startup; each rank processes
# its assigned healpix files, fitting only the targets listed in the samplefile.
# Per-file vs per-object parallelism is auto-detected from average targets/file.
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
samplefile=${4:?'samplefile is required'}
stage=${5:-'fastspec'}
input_redshifts=${6:-''}

outdir_data=${OUTDIR_DATA:-${PSCRATCH}/fastspecfit/data}

args="--outdir-data=${outdir_data} --specprod=${specprod} --mp=${mp}"
args+=" --samplefile=${samplefile}"

case "${stage}" in
    fastphot) args+=" --fastphot" ;;
esac

[[ "${input_redshifts}" == 'true' ]] && args+=" --input-redshifts"

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
