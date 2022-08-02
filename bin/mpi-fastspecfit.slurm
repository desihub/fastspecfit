#! /bin/bash
#SBATCH -A desi
#SBATCH --module=mpich
#SBATCH -C cpu
#SBATCH -L cfs
#SBATCH -o /global/cfs/cdirs/desi/spectro/fastspecfit/guadalupe/logs/fastspec-guadalupe-sv2.log.%j
#SBATCH --image=docker:desihub/fastspecfit:v1.0.0
#SBATCH --mail-user=jmoustakas@siena.edu
#SBATCH --mail-type=ALL
#SBATCH -q regular
#SBATCH -N 64
#SBATCH -n 64
#SBATCH -t 02:00:00

# sbatch mpi-fastspecfit.slurm

time srun -n 64 -c 128 --kill-on-bad-exit=0 --no-kill shifter /global/homes/i/ioannis/code/desihub/fastspecfit/bin/mpi-fastspecfit.sh fastspec guadalupe 128 healpix main
