#! /bin/bash
#SBATCH -A desi
#SBATCH -C cpu
#SBATCH -o /global/cfs/cdirs/desi/spectro/fastspecfit/iron/logs/fastspec-iron-main-dark-03.log.%j
#SBATCH --mail-user=jmoustakas@siena.edu
#SBATCH --mail-type=ALL
#SBATCH -q regular
#SBATCH -N 64
#SBATCH -n 64
#SBATCH -t 03:00:00

# sbatch /global/homes/i/ioannis/code/desihub/fastspecfit/bin/mpi-fastspecfit.slurm

#time srun -n 16 -c 128 /global/homes/i/ioannis/code/desihub/fastspecfit/bin/mpi-fastspecfit.sh makeqa fuji 128 healpix
time srun -n 64 -c 128 /global/homes/i/ioannis/code/desihub/fastspecfit/bin/mpi-fastspecfit.sh fastspec iron 128 healpix main dark
