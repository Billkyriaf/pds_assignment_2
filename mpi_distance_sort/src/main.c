#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <limits.h>

#include "data/points.h"
#include "data/structs.h"
#include "data/read_points.h"
#include "processes/masterProcess.h"
#include "processes/slaveProcess.h"


/**
 * Function that test the intermediate results of the program. The test is simple it calculates the distance and
 * compares it to the current median if the value is not correct it prints the error message.
 *
 * @param info The info struct of the process
 * @param min_rank  The minimum rank of the processes team
 * @param max_rank  The maximum rank of the processes team
 */
void testResultIntermediate(Info *info, int min_rank, int max_rank){
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
 * This function tests the over all results. The method for testing is every process sends its biggest distance to the
 * next process. If the biggest distance of the ith process is smaller
 * that the smallest distance of the (i + 1)th process the problem is solved correctly.
 * @param info The info struct of every process
 */
void testOverAllResults(Info *info){
    // Distance does not need to be an array because we only care for the biggest and smallest distance
    double distance;

    // The biggest and smallest distance for every process
    double smallest = INT_MAX;
    double biggest = INT_MIN;

    // Start calculating the distances from the pivot for all the points
    for (int i = 0; i < info->pointsPerProcess; ++i) {
        distance = findDistance(&info->pivot, &info->points[i], info->pointsDimension);

        // At the same time figure out the biggest
        if (distance > biggest){
            biggest = distance;
        }

        // ...and smallest distance
        if (distance < smallest){
            smallest = distance;
        }
    }

    // If this is the first process it only needs to send the max distance to the next one
    if (info->initial_rank == 0){
        MPI_Send(&biggest, 1, MPI_DOUBLE, info->initial_rank + 1, 100, MPI_COMM_WORLD);

    } else if (info->initial_rank == info->initial_size - 1){  // If this is the last process it only needs to receive the max distance of the previous one
        double prev_max;
        MPI_Recv(&prev_max, 1, MPI_DOUBLE, info->initial_rank - 1, 100, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        if (prev_max < smallest){
            printf("        Rank %d PASS\n", info->initial_rank - 1);
            printf("        Rank %d PASS\n", info->initial_rank);

        } else {
            printf("        Rank %d FAIL\n", info->initial_rank);
        }

    } else {  // Every other process needs to send the max distance to the next process and receive the max distance from the previous
        MPI_Request request_1;
        MPI_Request request_2;
        double prev_max;

        MPI_Irecv(&prev_max, 1, MPI_DOUBLE, info->initial_rank - 1, 100, MPI_COMM_WORLD, &request_1);
        MPI_Isend(&biggest, 1, MPI_DOUBLE, info->initial_rank + 1, 100, MPI_COMM_WORLD, &request_2);

        MPI_Wait(&request_1, NULL);
        MPI_Wait(&request_2, NULL);

        if (prev_max < smallest){
            printf("        Rank %d PASS\n", info->initial_rank - 1);

        } else {
            printf("        Rank %d FAIL\n", info->initial_rank - 1);
        }
    }

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

    // Test the intermediate results.
    //testResultIntermediate(info, min_rank, max_rank);

    // Recursion finish statement.
    // The recursion stops if the function is called with two processes in the processes group
    if (max_rank - min_rank == 1){
        MPI_Comm_free(&communicator);
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

        //printf("New Rank %d from %d upperComm. Old rank %d\n", info->world_rank, info->world_size, old_rank);

        // Recursive call of the function
        sort_points(0, 0, info->world_size - 1, info, upperComm);

    } else {
        MPI_Comm lowerComm;  // The half bottom ranks will go here

        MPI_Comm_split(communicator, 2, info->world_rank, &lowerComm);

        // Get the number of processes
        MPI_Comm_size(lowerComm, &info->world_size);

        // Get the rank of the process
        MPI_Comm_rank(lowerComm, &info->world_rank);

        //printf("New Rank %d from %d lowerComm. Old rank %d\n", info->world_rank, info->world_size, old_rank);

        // Recursive call of the function
        sort_points(0, 0, info->world_size - 1, info, lowerComm);
    }

    // Free the communicators
    if (communicator != MPI_COMM_WORLD){
        MPI_Comm_free(&communicator);
    }

}

/**
 * Checks the Command line arguments and prints the welcome message
 * @param argc Number of arguments
 * @param argv Arguments
 * @param info The info struct of the program
 */
void checkArguments(int argc, char **argv, Info info){
    // Help message
    if (argc == 2 && strcmp(argv[1], "help") == 0){
        if (info.initial_rank == 0){
            printf(" __          __    _                                _           __  __  _____  _____   _    _        _        \n");
            printf(" \\ \\        / /   | |                              | |         |  \\/  ||  __ \\|_   _| | |  | |      | |       \n");
            printf("  \\ \\  /\\  / /___ | |  ___  ___   _ __ ___    ___  | |_  ___   | \\  / || |__) | | |   | |__| |  ___ | | _ __  \n");
            printf("   \\ \\/  \\/ // _ \\| | / __|/ _ \\ | '_ ` _ \\  / _ \\ | __|/ _ \\  | |\\/| ||  ___/  | |   |  __  | / _ \\| || '_ \\ \n");
            printf("    \\  /\\  /|  __/| || (__| (_) || | | | | ||  __/ | |_| (_) | | |  | || |     _| |_  | |  | ||  __/| || |_) |\n");
            printf("     \\/  \\/  \\___||_| \\___|\\___/ |_| |_| |_| \\___|  \\__|\\___/  |_|  |_||_|    |_____| |_|  |_| \\___||_|| .__/ \n");
            printf("                                                                                                       | |    \n");
            printf("                                                                                                       |_|    \n");

            printf("USAGE:\n\n");
            printf("    make run_mpi ./path/to/data [option]\n\n");
            printf("    Options:\n");
            printf("        v: Verbose. Print welcome message and other information\n");
            printf("        nv: Non verbose. Print only required information\n");
            printf("        help: Display help message\n\n\n");
            MPI_Finalize();
            exit(0);

        } else {
            MPI_Finalize();
            exit(0);
        }
    }

    if (argc != 3){
        // Error message wrong number of arguments
        if (info.initial_rank == 0){
            printf(" __          __    _                                _           __  __  _____  _____ \n");
            printf(" \\ \\        / /   | |                              | |         |  \\/  ||  __ \\|_   _|\n");
            printf("  \\ \\  /\\  / /___ | |  ___  ___   _ __ ___    ___  | |_  ___   | \\  / || |__) | | |  \n");
            printf("   \\ \\/  \\/ // _ \\| | / __|/ _ \\ | '_ ` _ \\  / _ \\ | __|/ _ \\  | |\\/| ||  ___/  | |  \n");
            printf("    \\  /\\  /|  __/| || (__| (_) || | | | | ||  __/ | |_| (_) | | |  | || |     _| |_ \n");
            printf("     \\/  \\/  \\___||_| \\___|\\___/ |_| |_| |_| \\___|  \\__|\\___/  |_|  |_||_|    |_____|\n\n\n");
            printf("Illegal number of arguments.\n\n");
            printf("USAGE:\n\n");
            printf("    make run_mpi ./path/to/data [option]\n\n");
            printf("    Options:\n");
            printf("        v: Verbose. Print welcome message and other information\n");
            printf("        nv: Non verbose. Print only required information\n");
            printf("        help: Display help message\n\n\n");
            MPI_Finalize();
            exit(0);

        } else {
            MPI_Finalize();
            exit(0);
        }

    } else {
        // Error message don't know arguments
        if (strcmp(argv[2], "v") != 0 && strcmp(argv[2], "nv") != 0 && strcmp(argv[2], "help") != 0){
            if (info.initial_rank == 0){
                printf(" __          __    _                                _           __  __  _____  _____ \n");
                printf(" \\ \\        / /   | |                              | |         |  \\/  ||  __ \\|_   _|\n");
                printf("  \\ \\  /\\  / /___ | |  ___  ___   _ __ ___    ___  | |_  ___   | \\  / || |__) | | |  \n");
                printf("   \\ \\/  \\/ // _ \\| | / __|/ _ \\ | '_ ` _ \\  / _ \\ | __|/ _ \\  | |\\/| ||  ___/  | |  \n");
                printf("    \\  /\\  /|  __/| || (__| (_) || | | | | ||  __/ | |_| (_) | | |  | || |     _| |_ \n");
                printf("     \\/  \\/  \\___||_| \\___|\\___/ |_| |_| |_| \\___|  \\__|\\___/  |_|  |_||_|    |_____|\n\n\n");
                printf("Can not understand %s argument\n\n", argv[2]);
                printf("USAGE:\n\n");
                printf("    make run_mpi ./path/to/data [option]\n\n");
                printf("    Options:\n");
                printf("        v: Verbose. Print welcome message and other information\n");
                printf("        nv: Non verbose. Print only required information\n");
                printf("        help: Display help message\n\n\n");
                MPI_Finalize();
                exit(0);

            } else {
                MPI_Finalize();
                exit(0);
            }

        } else if(strcmp(argv[2], "help") == 0){  // Welcome message only runs on verbose
            if (info.initial_rank == 0){
                printf(" __          __    _                                _           __  __  _____  _____   _    _        _        \n");
                printf(" \\ \\        / /   | |                              | |         |  \\/  ||  __ \\|_   _| | |  | |      | |       \n");
                printf("  \\ \\  /\\  / /___ | |  ___  ___   _ __ ___    ___  | |_  ___   | \\  / || |__) | | |   | |__| |  ___ | | _ __  \n");
                printf("   \\ \\/  \\/ // _ \\| | / __|/ _ \\ | '_ ` _ \\  / _ \\ | __|/ _ \\  | |\\/| ||  ___/  | |   |  __  | / _ \\| || '_ \\ \n");
                printf("    \\  /\\  /|  __/| || (__| (_) || | | | | ||  __/ | |_| (_) | | |  | || |     _| |_  | |  | ||  __/| || |_) |\n");
                printf("     \\/  \\/  \\___||_| \\___|\\___/ |_| |_| |_| \\___|  \\__|\\___/  |_|  |_||_|    |_____| |_|  |_| \\___||_|| .__/ \n");
                printf("                                                                                                       | |    \n");
                printf("                                                                                                       |_|    \n");

                printf("USAGE:\n\n");
                printf("    make run_mpi ./path/to/data [option]\n\n");
                printf("    Options:\n");
                printf("        v: Verbose. Print welcome message and other information\n");
                printf("        nv: Non verbose. Print only required information\n");
                printf("        help: Display help message\n\n\n");
                MPI_Finalize();
                exit(0);

            } else {
                MPI_Finalize();
                exit(0);
            }

        } else if(strcmp(argv[2], "v") == 0 && info.initial_rank == 0){  // Welcome message only runs on verbose
            printf(" __          __    _                                _           __  __  _____  _____ \n");
            printf(" \\ \\        / /   | |                              | |         |  \\/  ||  __ \\|_   _|\n");
            printf("  \\ \\  /\\  / /___ | |  ___  ___   _ __ ___    ___  | |_  ___   | \\  / || |__) | | |  \n");
            printf("   \\ \\/  \\/ // _ \\| | / __|/ _ \\ | '_ ` _ \\  / _ \\ | __|/ _ \\  | |\\/| ||  ___/  | |  \n");
            printf("    \\  /\\  /|  __/| || (__| (_) || | | | | ||  __/ | |_| (_) | | |  | || |     _| |_ \n");
            printf("     \\/  \\/  \\___||_| \\___|\\___/ |_| |_| |_| \\___|  \\__|\\___/  |_|  |_||_|    |_____|\n\n\n");


            printf("IMPORTANT NOTES:\n\n");
            printf("    1. The whole project is designed to work with datasets that the number of points they\n");
            printf("       have is a power of 2. You don't have to do anything yourself but keep in mind that\n");
            printf("       if you run the program it will read the power of 2 points that is the closest to the\n");
            printf("       number of points in the binary file.\n\n");
            printf("    2. The number of process must also be a power of 2. This you must ensure yourself from\n");
            printf("       the host file. If you don't change anything the default number of processes is 8.\n\n\n\n\n");
        }
    }
}

/**
 * The main function of the program
 * @param argc The number of the cmd arguments
 * @param argv The cmd arguments. The function takes only one argument the path to the binary data file.
 */
int main(int argc, char **argv) {

    // Initialize the MPI communication
    MPI_Init(&argc, &argv);

    double start, end;  // Time measuring

    Info info;  // The info struct holds useful information for every process that need to remain between function calls

    // Get the number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &info.world_size);

    // Get the rank of the process
    MPI_Comm_rank(MPI_COMM_WORLD, &info.world_rank);

    // The world size and rank for every process are depended on the communicator. As the communicator changes so do
    // those values too. The initial rank and size represent the values for each rank in the MPI_COM_WORLD communicator
    info.initial_rank = info.world_rank;  // Just for testing
    info.initial_size = info.world_size; // Just for testing

    // Check the cmd arguments
    checkArguments(argc, argv, info);

    // Read data from the binary file
    info.points = readFromFile(  // MEMORY
            argv[1],
            info.world_rank,
            info.world_size,
            &info.pointsDimension,
            &info.pointsPerProcess
    );

    // Time measurement starts after the points are read
    if (info.initial_rank == 0){
        start = MPI_Wtime();
    }


    // Create the pivot for every process. The pivot remains constant for the execution of the program, so it is only
    // selected once
    info.pivot.coordinates = (double *) malloc(info.pointsDimension * sizeof(double));  // MEMORY


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

    // Wait here for all the processes to finish before time measurement
    MPI_Barrier(MPI_COMM_WORLD);

    // Time measurement ends here
    if (info.initial_rank == 0){
        end = MPI_Wtime();

        printf("\nTime for execution: %.6f\n", end - start);

        // Start the timer for testing time
        start = MPI_Wtime();

        printf("\n    Validation results:\n\n");
    }

    // Wait here for all the processes to finish before time measurement
    MPI_Barrier(MPI_COMM_WORLD);

    // Test the final result
    testOverAllResults(&info);

    // Wait here for all the processes to finish before time measurement
    MPI_Barrier(MPI_COMM_WORLD);

    if (info.initial_rank == 0){
        end = MPI_Wtime();

        printf("\n    Time for validation: %.6f\n\n", end - start);
    }

    // Deallocate any allocated memory
    for (int i = 0; i < info.pointsPerProcess; ++i) {
        free(info.points[i].coordinates);
    }
    free(info.points);
    free(info.pivot.coordinates);


    // Finalize the MPI environment.
    MPI_Finalize();

    return 0;
}