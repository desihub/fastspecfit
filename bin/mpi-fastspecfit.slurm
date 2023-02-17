#! /bin/bash
#SBATCH -A desi
#SBATCH -C cpu
#SBATCH -o /global/cfs/cdirs/desi/spectro/fastspecfit/fuji/logs/fastspec-fuji-makeqa1.log.%j
#SBATCH --mail-user=jmoustakas@siena.edu
#SBATCH --mail-type=ALL
#SBATCH -q regular
#SBATCH -N 16
#SBATCH -n 16
#SBATCH -t 00:30:00

# sbatch /global/homes/i/ioannis/code/desihub/fastspecfit/bin/mpi-fastspecfit.slurm

#time srun -n 16 -c 128 /global/homes/i/ioannis/code/desihub/fastspecfit/bin/mpi-fastspecfit.sh makeqa fuji 128 healpix
time srun -n 16 -c 128 /global/homes/i/ioannis/code/desihub/fastspecfit/bin/mpi-fastspecfit.sh fastspec-test fuji 128 healpix
