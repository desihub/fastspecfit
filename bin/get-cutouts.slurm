#! /bin/bash
#SBATCH -A desi
#SBATCH -C cpu
#SBATCH -o /global/cfs/cdirs/desi/spectro/fastspecfit/fuji/v2.0/logs/cutouts-fuji-07.log.%j
#SBATCH --mail-user=jmoustakas@siena.edu
#SBATCH --image=dstndstn/viewer-cutouts
#SBATCH --mail-type=ALL
#SBATCH -q regular
#SBATCH -N 2
#SBATCH -n 2
#SBATCH -t 04:00:00

# sbatch /global/homes/i/ioannis/code/desihub/fastspecfit/bin/get-cutouts.slurm

time srun -n 2 -c 128 /global/homes/i/ioannis/code/desihub/fastspecfit/bin/get-cutouts.sh fuji 128 healpix sv3 dark
