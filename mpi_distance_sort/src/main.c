#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#include "data/points.h"
#include "data/read_points.h"
#include "processes/masterProcess.h"
#include "processes/slaveProcess.h"
#include "main.h"


/**
 * Function that test the results of the program. The test is simple it calculates the distance and compares it to the
 * current median if the value is not correct it prints the error message.
 *
 * @param info The info struct of the process
 * @param min_rank  The minimum rank of the processes team
 * @param max_rank  The maximum rank of the processes team
 */
void testResult(Info *info, int min_rank, int max_rank){
    // Memory allocation for the distance vector
    double *distVector = (double*) malloc(info->pointsPerProcess * sizeof (double));


    // Start calculating the distances from the pivot for all the points
    for (int i = 0; i < info->pointsPerProcess; ++i) {
        distVector[i] = findDistance(&info->pivot, &info->points[i], info->pointsDimension);
    }

    if (info->world_rank >= (max_rank - min_rank + 1) / 2){
        for (int i = 0; i < info->pointsPerProcess; ++i) {
            if (distVector[i] < info->median){
                printf("Rank %d FAILED the test\n", info->initial_rank);
                return;
            }
        }

        printf("Rank %d PASS\n", info->initial_rank);

    } else {
        for (int i = 0; i < info->pointsPerProcess; ++i) {
            if (distVector[i] > info->median){
                printf("Rank %d FAILED the test\nMedian was: %.10f and distance was: %.10f\n", info->initial_rank, info->median, distVector[i]);
                return;
            }
        }

        printf("Rank %d PASS\n", info->initial_rank);
    }

    free(distVector);
}


/**
 * This is the main function that starts the shorting process. This function is call recursively until all the processes
 * have the correct points
 * @param master_rank The rank of the master process for this call
 * @param min_rank This smallest rank among the processes of the group. The initial group of processes contains every
 *                 process and it is gradually split to two, four, eight etc. parts with every call of the function
 * @param max_rank This is the biggest rank among the processes of the group
 * @param info The info struct of the current process
 * @param communicator The MPI communicator object. Initially it is MPI_COMM_WORLD
 */
void sort_points(int master_rank, int min_rank, int max_rank, Info *info, MPI_Comm communicator){
    // The master of every group calls the master function. Everyone else calls the slave function.
    if  (info->world_rank == master_rank){
        masterProcess(master_rank, min_rank, max_rank, info, communicator);
    } else {
        slaveProcess(master_rank, min_rank, max_rank, info, communicator);
    }

    // Function that tests the results for a given group of processes
    testResult(info, min_rank, max_rank);

    // Recursion finish statement.
    // The recursion stops if the function is called with only one process in the processes group
    if (max_rank - min_rank == 1){
        return;
    }


    // If the process belongs to the bigger half ...
    if (info->world_rank >= (max_rank - min_rank + 1) / 2){
        MPI_Comm upperComm;  // The half bigger ranks will go here

        // Call the split comm function with color
        MPI_Comm_split(communicator, 1, info->world_rank, &upperComm);

        // Get the number of processes
        MPI_Comm_size(upperComm, &info->world_size);

        // Get the rank of the process
        MPI_Comm_rank(upperComm, &info->world_rank);

//        printf("New Rank %d from %d upperComm. Old rank %d\n", info->world_rank, info->world_size, old_rank);

        sort_points(0, 0, info->world_size - 1, info, upperComm);

    } else {
        MPI_Comm lowerComm;  // The half bottom ranks will go here

        MPI_Comm_split(communicator, 2, info->world_rank, &lowerComm);

        // Get the number of processes
        MPI_Comm_size(lowerComm, &info->world_size);

        // Get the rank of the process
        MPI_Comm_rank(lowerComm, &info->world_rank);

//        printf("New Rank %d from %d lowerComm. Old rank %d\n", info->world_rank, info->world_size, old_rank);

        sort_points(0, 0, info->world_size - 1, info, lowerComm);
    }


    // FIXME fix the color property to work with more than 4 processes
    // FIXME terminate the communicators
}


/**
 * The main function of the program
 * @param argc The number of the cmd arguments
 * @param argv The cmd arguments. The function takes only one argument the path to the binary data file.
 * @return
 */
int main(int argc, char **argv) {

    // check for the command line arguments
    if (argc != 2){
        printf("Illegal number of arguments. The function takes one argument: The path to the binary data file. Exiting...\n");
        return -1;
    }

    // Initialize the MPI communication
    MPI_Init(&argc, &argv);

    Info info;  // The info struct holds useful information for every process that need to remain between function calls

    // Get the number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &info.world_size);

    // Get the rank of the process
    MPI_Comm_rank(MPI_COMM_WORLD, &info.world_rank);

    info.initial_rank = info.world_rank;  // Just for debugging

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