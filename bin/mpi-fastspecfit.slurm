#! /bin/bash
#SBATCH -A desi
#SBATCH --module=mpich-cle6
#SBATCH -C haswell
#SBATCH -L cfs,SCRATCH
#SBATCH -o fastspec-everest-main-dark.log.%j
#SBATCH --image=docker:desihub/fastspecfit:v0.2
#SBATCH --mail-user=jmoustakas@siena.edu
#SBATCH --mail-type=ALL
#SBATCH -p regular
#SBATCH -N 64
#SBATCH -n 64
#SBATCH -t 04:00:00

# sbatch mpi-fastspecfit.slurm

time srun -n 64 -c 32 --kill-on-bad-exit=0 --no-kill shifter /global/homes/i/ioannis/repos/desihub/fastspecfit/bin/mpi-fastspecfit.sh fastspec 32 main dark
