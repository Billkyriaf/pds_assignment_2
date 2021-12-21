#!/bin/bash
#SBATCH --partition=testing
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=4
#SBATCH --cpus-per-task=4
#SBATCH --time=1:00

module load gcc/9.2.0 openmpi/3.1.3

mpirun --mca opal_warn_on_missing_libcuda 0 -hostfile ./hostfile build/mpi_main.out ./datasets/data.bin