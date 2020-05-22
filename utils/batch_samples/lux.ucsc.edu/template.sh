#!/bin/bash
#SBATCH --job-name=slugCode2          # Job name
#SBATCH --partition=lee               # queue for job submission
#SBATCH --account=lee                 # account for job submission
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=[[EMAIL]]         # Where to send mail
#SBATCH --ntasks=64                   # Number of MPI ranks
#SBATCH --nodes=2                     # Number of nodes
#SBATCH --ntasks-per-node=32          # How many tasks on each node
#SBATCH --time=12:00:00               # Time limit hrs:min:sec
#SBATCH --output=log/slug2_%j.log     # Standard output and error log

pwd; hostname; date

echo "Running program on $SLURM_JOB_NUM_NODES nodes with $SLURM_NTASKS total tasks, with each node getting $SLURM_NTASKS_PER_NODE running on cores."

# intel compiler stack
module load intel/ifort intel/impi intel/icc hdf5/1.10.6-parallel blas/gcc/64/3.8.0 lapack/gcc/64/3.8.0 slurm/18.08.4

mpirun -np 64 [[PATH/TO/SLUGCODE2]] [PATH/TO/INITFILE.init]

date
