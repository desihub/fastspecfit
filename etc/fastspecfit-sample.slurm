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

# Load the DESI software stack and set environment variables.
# SLURM_SUBMIT_DIR is set by Slurm to the directory from which sbatch was run.
# Submit this script from the directory that also contains fastspecfit-env.sh
# and fastspecfit-sample.sh, or set SCRIPT_DIR to an absolute path instead.
SCRIPT_DIR=${SLURM_SUBMIT_DIR}
source "${SCRIPT_DIR}/fastspecfit-env.sh"

# To use a local checkout instead of the module, uncomment and edit:
#export PYTHONPATH=/path/to/fastspecfit/py:$PYTHONPATH
#export PATH=/path/to/fastspecfit/bin:$PATH

# ---------------------------------------------------------------------------
# Edit the variables below before submitting.
# ---------------------------------------------------------------------------
N=16                    # number of nodes (must match --nodes above)
mp=16                   # multiprocessing workers per MPI rank
specprod=loa            # spectroscopic production
samplefile=${PSCRATCH}/fastspecfit/samples/my-sample.fits
stage=fastspec          # fastspec or fastphot
input_redshifts=''      # set to 'true' to use redshifts from samplefile Z column
# OUTDIR_DATA is set by fastspecfit-env.sh; override here if needed:
#export OUTDIR_DATA=/path/to/custom/output
# ---------------------------------------------------------------------------

mpiscript=$(type -p mpi-fastspecfit)

# Warm up Numba cache on a single rank before the parallel run.
# Compiled artifacts go to $NUMBA_CACHE_DIR (shared Lustre), so all
# subsequent ranks load from cache rather than recompiling.
# Use a single healpix pixel from the samplefile to guarantee exactly
# one file is processed (--ntargets limits per file, not globally).
echo "Warming up Numba cache on one rank..."
read wu_survey wu_program wu_healpix <<< $(python -c "
from astropy.table import Table; t = Table.read('${samplefile}')[0]
print(t['SURVEY'], t['PROGRAM'], t['HEALPIX'])")
srun --nodes=1 --ntasks=1 --cpus-per-task=2 --cpu-bind=cores \
    $mpiscript \
    --specprod=${specprod} --coadd-type=healpix \
    --survey=${wu_survey} --program=${wu_program} --healpix=${wu_healpix} \
    --mp=1 --ntargets=1 --nmonte=0 --nompi \
    --outdir-data=/tmp/fastspecfit-warmup --overwrite
echo "Warm-up complete."

sh "${SCRIPT_DIR}/fastspecfit-sample.sh" \
    ${N} ${mp} ${specprod} ${samplefile} ${stage} ${input_redshifts}
