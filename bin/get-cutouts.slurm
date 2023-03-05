#! /bin/bash
#SBATCH -A desi
#SBATCH -C cpu
#SBATCH -o /global/cfs/cdirs/desi/spectro/fastspecfit/fuji/logs/cutouts-fuji-01.log.%j
#SBATCH --mail-user=jmoustakas@siena.edu
#SBATCH --mail-type=ALL
#SBATCH -q regular
#SBATCH -N 8
#SBATCH -n 8
#SBATCH -t 01:00:00

# sbatch /global/homes/i/ioannis/code/desihub/fastspecfit/bin/get-cutouts.slurm

time srun -n 64 -c 128 /global/homes/i/ioannis/code/desihub/fastspecfit/bin/get-cutouts.sh fastspec fuji 128 healpix
