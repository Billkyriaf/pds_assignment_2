#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#include "data/points.h"
#include "data/read_points.h"
#include "processes/masterProcess.h"
#include "processes/slaveProcess.h"
#include "main.h"


void sort_points(int master_rank, int min_rank, int max_rank, Info *info, MPI_Comm communicator){
    // Recursion finish statement
    if (min_rank - max_rank == 0){
        return;
    }

    if  (info->world_rank == master_rank){
        masterProcess(master_rank, min_rank, max_rank, info, communicator);
    } else {
        slaveProcess(master_rank, min_rank, max_rank, info, communicator);
    }



    // TODO call the function recursively with new masters for half of the processes
}

int main(int argc, char **argv) {

    // check for the command line arguments
    if (argc != 2){
        printf("Wrong number of arguments exiting...\n");
        return -1;
    }

    // Initialize the MPI communication
    MPI_Init(&argc, &argv);

    Info info;

    // Get the number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &info.world_size);

    // Get the rank of the process
    MPI_Comm_rank(MPI_COMM_WORLD, &info.world_rank);

    // Read data from the binary file
    info.points = readFromFile(
            argv[1],
            info.world_rank,
            info.world_size,
            &info.pointsDimension,
            &info.pointsPerProcess
    );

    // Create the pivot for every process. The pivot remains constant for the execution of the program, so it is only
    // selected once
    info.pivot.coordinates = (double *) malloc(info.pointsDimension * sizeof(double));


    // If this is the master process...
    if (info.world_rank == 0) {
        // ...pick a pivot point from the points owned.
        Point tmp = info.points[0];

        // Deep copy the tmp to pivot. This is required because the info.points array will change between function calls
        // but the pivot needs to remain constant
        for (int i = 0; i < info.pointsDimension; ++i) {
            info.pivot.coordinates[i] = tmp.coordinates[i];
        }

        // Update the distance to the pivot
        info.pivot.distance = 0.0;
    }

    // Broadcast the pivot to all the other processes
    // Send the coordinates of the pivot point to all the processes
    MPI_Bcast(
            info.pivot.coordinates,
            info.pointsDimension,
            MPI_DOUBLE,
            0,
            MPI_COMM_WORLD
    );

    // Call the function initially with master_rank = 0, rank = world_rank, min_rank = 0, max_rank = world_size - 1
    sort_points(
            0,
            0,
            info.world_size - 1,
            &info,
            MPI_COMM_WORLD
    );


    // Deallocate any allocated memory
    for (int i = 0; i < info.pointsPerProcess; ++i) {
        free(info.points[i].coordinates);
    }
    free(info.points);


    // Finalize the MPI environment.
    MPI_Finalize();

    return 0;
}