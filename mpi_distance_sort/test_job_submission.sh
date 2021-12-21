#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --partition=testing
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=4

module load gcc openmpi


make build_mpi

mpirun -n 2 ./build/mpi_main.out data_file=datasets/bata.bin

