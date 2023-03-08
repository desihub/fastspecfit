#! /bin/bash
#SBATCH -A desi
#SBATCH -C cpu
#SBATCH -o /global/cfs/cdirs/desi/spectro/fastspecfit/fuji/v2.0/logs/fastspec-fuji-qa-03.log.%j
#SBATCH --mail-user=jmoustakas@siena.edu
#SBATCH --mail-type=ALL
#SBATCH -q regular
#SBATCH -N 32
#SBATCH -n 32
#SBATCH -t 02:00:00

# sbatch /global/homes/i/ioannis/code/desihub/fastspecfit/bin/mpi-fastspecfit.slurm

time srun -n 32 -c 128 /global/homes/i/ioannis/code/desihub/fastspecfit/bin/mpi-fastspecfit.sh makeqa fuji 128 healpix sv3
#time srun -n 64 -c 128 /global/homes/i/ioannis/code/desihub/fastspecfit/bin/mpi-fastspecfit.sh fastspec fuji 128 healpix
