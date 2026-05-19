#!/bin/bash -l
#SBATCH --account=desi
#SBATCH --qos=regular
#SBATCH --constraint=cpu
#SBATCH --mail-user=jmoustakas@siena.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=16
#SBATCH --time=04:00:00
#SBATCH --job-name=fastspec-loa-sample
#SBATCH --output=fastspec-loa-sample-%j.log

# Edit the variables below before submitting.
N=16                    # number of nodes (match --nodes above)
mp=16                   # multiprocessing workers per MPI rank
specprod=loa            # spectroscopic production
samplefile=${PSCRATCH}/fastspecfit/samples/my-sample.fits
stage=fastspec          # fastspec or fastphot
input_redshifts=''      # set to 'true' to use redshifts from samplefile Z column

mpiscript=$(type -p mpi-fastspecfit)

# Warm up Numba cache on a single rank before the parallel run.
# Compiled artifacts go to $NUMBA_CACHE_DIR (shared Lustre), so all
# subsequent ranks load from cache rather than recompiling.
echo "Warming up Numba cache on one rank..."
srun --nodes=1 --ntasks=1 --cpus-per-task=2 --cpu-bind=cores \
    $mpiscript \
    --specprod=${specprod} --samplefile=${samplefile} \
    --mp=1 --ntargets=1 --nompi --nolog \
    --outdir-data=/tmp/fastspecfit-warmup --overwrite
echo "Warm-up complete."

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
sh "${SCRIPT_DIR}/fastspecfit-sample.sh" \
    ${N} ${mp} ${specprod} ${samplefile} ${stage} ${input_redshifts}
