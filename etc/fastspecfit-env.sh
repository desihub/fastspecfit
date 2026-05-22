#!/bin/bash
# Environment setup for FastSpecFit production runs at NERSC.
# Source this script; do not execute it directly.
#
# Usage:
#   source fastspecfit-env.sh [fastspecfit_version]
#
# The optional argument overrides the default FastSpecFit module version.
# All other module versions can be overridden by setting the corresponding
# environment variables before sourcing:
#   DESIUTIL_VERSION, DESISPEC_VERSION, DESITARGET_VERSION,
#   SPECLITE_VERSION

FASTSPECFIT_VERSION=${1:-"3.4.0"}

DESIUTIL_VERSION=${DESIUTIL_VERSION:-"3.6.1"}
DESISPEC_VERSION=${DESISPEC_VERSION:-"0.71.2"}
DESITARGET_VERSION=${DESITARGET_VERSION:-"4.7.2"}
SPECLITE_VERSION=${SPECLITE_VERSION:-"v1.0.0"}

echo "Loading DESI software stack"
source /dvs_or/common/software/desi/desi_environment.sh main
module swap desiutil/${DESIUTIL_VERSION}
module swap desispec/${DESISPEC_VERSION}
module swap desitarget/${DESITARGET_VERSION}
module swap speclite/${SPECLITE_VERSION}

echo "Loading FastSpecFit/${FASTSPECFIT_VERSION}"
module load fastspecfit/${FASTSPECFIT_VERSION}
# To use a development install instead:
#echo "Loading local check-out of FastSpecFit"
#export PYTHONPATH=/dvs_or/common/software/desi/users/ioannis/fastspecfit/py:$PYTHONPATH
#export PATH=/dvs_or/common/software/desi/users/ioannis/fastspecfit/bin:$PATH

export DESI_SPECTRO_REDUX=/dvs_ro/cfs/cdirs/desi/spectro/redux
export DUST_DIR=/dvs_ro/cfs/cdirs/cosmo/data/dust/v0_1
export FPHOTO_DIR=/dvs_ro/cfs/cdirs/desi/external/legacysurvey/dr9
export FTEMPLATES_DIR=/dvs_ro/cfs/cdirs/desi/public/external/templates/fastspecfit

# Output data directory. Override this variable before sourcing this script if
# the default location is not suitable (e.g. export OUTDIR_DATA=/path/to/output).
export OUTDIR_DATA=${OUTDIR_DATA:-${PSCRATCH}/fastspecfit/data}

# Shared Numba JIT cache: all MPI ranks and spawned workers read compiled
# artifacts from this directory.  Using $PSCRATCH ensures it is user-writable
# and visible from every compute node.  The .slurm job templates include a
# single-rank warm-up srun step (--ntasks=1 --ntargets=1 --nompi) before the
# full production run to pre-populate this cache and avoid a compilation
# stampede across all ranks at startup.
export NUMBA_CACHE_DIR=${NUMBA_CACHE_DIR:-${PSCRATCH}/fastspecfit/${FASTSPECFIT_VERSION}/numba-cache}
mkdir -p "${NUMBA_CACHE_DIR}"

echo "DESI_SPECTRO_REDUX=${DESI_SPECTRO_REDUX}"
echo "DUST_DIR=${DUST_DIR}"
echo "FPHOTO_DIR=${FPHOTO_DIR}"
echo "FTEMPLATES_DIR=${FTEMPLATES_DIR}"
echo "NUMBA_CACHE_DIR=${NUMBA_CACHE_DIR}"
echo "OUTDIR_DATA=${OUTDIR_DATA}"

export FASTSPECFIT_ENV_LOADED=1
