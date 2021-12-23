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
  - [4.1. Run individual files](#41-run-individual-files)
  - [4.2. Run multiple files at once](#42-run-multiple-files-at-once)
  - [4.3. Create random data](#43-create-random-data)

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

### 4.1. Run individual files

To run individual data files a [Makefile](mpi_distance_sort/Makefile) is provided. To build and run the application `cd` to the [root](mpi_distance_sort) directory of the project and follow these steps:

1. `make clean`: Clean the `build` directory.
2. `make build_mpi`: Build the binary executable.
3. `make run_mpi args`: Run the executable with provided arguments.
      
      ```
      Mandatory arguments:

        argv[1]: The path to the data file
        argv[2]: One of the following:
          v: Verbose output for the executable
          nv: Non verbose output for the executable (Minimum information are printed to the console)  
          help: Displays help message and exits
      ```
To adjust the number of hosts edit the slots in the [hostfile](mpi_distance_sort/hostfile).

Examples:
      
`make run_mpi datasets/data.bin v  # This will run the program with verbose output`

`make run_mpi help  # This will display help message and exit`
      
 
### 4.2. Run multiple files at once

To run multiple files at once use the [`benchmarks.sh`](mpi_distance_sort/benchmarks.sh). Give `execute` rights to the script file with `chmod +x benchmarks.sh` and run the script with `./benchmarks.sh`.

The script accepts certain arguments described bellow:

```
Usage ./benchmarks.sh [OPTION]...

Runs a series of tests with the MPI executable and the data files found in the directory specified"
    
OPTIONS:
    -d {directory}              The directory containing the datasets. This is mandatory and the
                                data files must have .bin extension.

    -v true || false || one     Verbose output for the executable. Values:
                                    true: All runs are verbose 
                                    false: All runs are non verbose
                                    one: Only the first run is verbose (Default if not provided)

    -n {processes}              The number of processes. Must be a power of 2!! Default is 8.

    -h                          Help message for the executable

    -u                          Usage for this script

Example usage:
    ./benchmarks.sh -v -d ./datasets

```

### 4.3. Create random data

To to create random data files a [Makefile](mpi_distance_sort/Makefile) is provided. To build and run the application `cd` to the [root](mpi_distance_sort) directory of the project and follow these steps:

1. `make create_datafile`: Creates a random binary data file with `dimension` `20` and `8388608` points. The file is saved in the `datasets` directory with name `data.bin`. This file is quite big (1.3GB) so you may want to adjust this in the `Makefile`. Simply change the numbers in the `create_datafile` target (line 76). The first argument is the dimension of the space and the second is the number of points.
