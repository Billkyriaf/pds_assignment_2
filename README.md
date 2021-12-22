<div id="top"></div>

<br />
<div align="center">
  <h1 align="center">Parallel and Distributed Systems Assignment 2</h1>
  <h3 align="center">Aristotle University of Thessaloniki</h3>
  <h4 align="center">School of Electrical & Computer Engineering</h4>
  <p align="center">
    Contributors: Kyriafinis Vasilis
    <br />
    Winter Semester 2021 - 2022
    <br />
    <br />
    <br />
    <br />
  </p>
</div>

- [1. About this project](#1-about-this-project)
- [2. Getting started](#2-getting-started)
- [3. Dependencies](#3-dependencies)
  - [3.1. Make](#31-make)
  - [3.2. OpenMPI](#32-openmpi)
- [4. Compile and run](#4-compile-and-run)

## 1. About this project

The objective of this project is to distribute n-dimensional points to multiple processes depending on the distance of each point from a randomly selected pivot point. The closest points to the pivot must end up to the first process and the furthest points from the pivot must end up to the last process. There is no need for the points of a single process to be sorted. 

The framework used for the implementation is openMPI. OpenMPI is a messaging system allowing processes running on different processors to communicate and exchange information.


## 2. Getting started

To setup this repository on your local machine run the following commands on the terminal:

```console
git clone https://github.com/Billkyriaf/pds_assignment_2.git
```

Or alternatively [*download*](https://github.com/Billkyriaf/pds_assignment_2/archive/refs/heads/main.zip) and extract the zip file of the repository
<br/>
<br/>

## 3. Dependencies
### 3.1. Make

This project uses make utilities to build and run the executables.

### 3.2. OpenMPI

The version of openMPI used in this project is `(Open MPI) 4.1.2`. You can download it and find more information [here](https://www.open-mpi.org/software/ompi/v4.1/).

## 4. Compile and run

To build the application `cd` to the [root](mpi_distance_sort) directory of the project. There you will find the [Makefile](mpi_distance_sort/Makefile). The Makefile provides all the required commands to build and run the application on your local machine. The targets provided from the Makefile are:

1. `make build_mpi`: Compiles all the files to a single executable. The produced binary is saved in the `build` directory
2. `make run_mpi data_file=path/to/binary/data`: Runs the application with `mpirun`. The `data_file=...` argument is required. 
3. `make debug_mpi`: Runs the application with `mpi_run` and `valgrind`.
4. `make clean`: Cleans the `build` directory.
5. `make create_datafile`: Creates a random binary data file with `dimension` `20` and `8388608` points. The file is saved in the `datasets` directory with name `data.bin`.


To adjust the number of hosts edit the slots in the [hostfile](mpi_distance_sort/hostfile).