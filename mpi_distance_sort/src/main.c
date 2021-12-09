#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <unistd.h>

#include "data/points.h"
#include "data/read_points.h"

typedef struct info{
    int world_size;
    int world_rank;

    Point *points;  // All the points for each process
    int pointsPerProcess;  // The total number of points each process has
    int pointsDimension;  // The dimension of the space for the points (number of coordinates)
} Info;

void sort_points(int rank, int master_rank, int min_rank, int max_rank, Info *info){
    // Recursion finish statement
    if (min_rank - max_rank == 0){
        return;
    }

    Point pivot;
    if  (rank == master_rank){
        // Allocate space for the distances from all the other processes
        double *distVector = (double*) malloc(info->pointsPerProcess * info->world_size * sizeof (double*));

        // Allocate space for all the requests
        MPI_Request *requests;
        requests = (MPI_Request *) malloc((max_rank - min_rank) * sizeof (MPI_Request));

        /*
        printf("Rank: %d\n", rank);
        for (int i = 0; i < info->pointsPerProcess; ++i) {
            printf("Point: %d\n", i);
            for (int j = 0; j < info->pointsDimension; ++j) {
                printf("  %f ", info->points[i].coordinates[j]);
            }
            printf("\n");
        }*/

        pivot = info->points[0];  // Pick a pivot point from the points owned if this is the master process

        // Send the coordinates of the pivot point to all the processes
        MPI_Bcast(
                pivot.coordinates,
                info->pointsDimension,
                MPI_DOUBLE,
                master_rank,
                MPI_COMM_WORLD
        );


        // Start receiving the distances from all the processes
        for (int i = 0; i < max_rank - min_rank; ++i) {
            MPI_Irecv(
                    distVector + (i + 1) * info->pointsPerProcess,
                    info->pointsPerProcess,
                    MPI_DOUBLE,
                    i + 1,
                    10,
                    MPI_COMM_WORLD,
                    requests + i
            );
        }

        // Start calculating the distances from the pivot for all the points
        for (int i = 0; i < info->pointsPerProcess; ++i) {
            distVector[i] = findDistance(&pivot, &info->points[i], info->pointsDimension);
            info->points->distance = distVector[i];
        }

        /*
        printf("Rank: %d\n", rank);
        for (int i = 0; i < info->pointsPerProcess; ++i) {
            printf("Distance: %f  ", distVectors[0][i]);
        }
        printf("\n");
        printf("\n");
        */

        // Wait for all the points to be received
        for (int i = 0; i < max_rank - min_rank; ++i) {
            MPI_Wait(requests + i, MPI_STATUS_IGNORE);
        }

        // Debug prints
        printf("Pivot: ");
        for (int i = 0; i < info->pointsDimension; ++i) {
            printf("%f ", pivot.coordinates[i]);
        }
        printf("\n");
        printf("\n");

        for (int i = 0; i <= max_rank; ++i) {
            printf("Distance table Rank %d: ", i);
            for (int j = 0; j < info->pointsPerProcess; ++j) {
                printf("%f ", distVector[i * info->pointsPerProcess + j]);
            }
            printf("\n");
        }

        // Calculate the median value and send it back to the other processes
        double median = findMedianValue(distVector, info->pointsPerProcess * info->world_size);
        printf("Median: %f\n", median);

        MPI_Bcast(
                &median,
                1,
                MPI_DOUBLE,
                master_rank,
                MPI_COMM_WORLD
        );

    } else {
        /*
        printf("Rank: %d\n", rank);
        for (int i = 0; i < info->pointsPerProcess; ++i) {
            printf("Point: %d\n", i);
            for (int j = 0; j < info->pointsDimension; ++j) {
                printf("  %f ", info->points[i].coordinates[j]);
            }
            printf("\n");
        }
        */

        // Memory allocation for the pivot point
        pivot.coordinates = (double *) malloc(info->pointsDimension * sizeof (double));

        // Memory allocation for the distance vector
        double *distVector = (double*) malloc(info->pointsPerProcess * sizeof (double));

        // Receive the coordinates of the pivot point from the master rank
        MPI_Bcast(
                pivot.coordinates,
                info->pointsDimension,
                MPI_DOUBLE,
                master_rank,
                MPI_COMM_WORLD
        );

        // Start calculating the distances from the pivot for all the points
        for (int i = 0; i < info->pointsPerProcess; ++i) {
            distVector[i] = findDistance(&pivot, &info->points[i], info->pointsDimension);
            info->points->distance = distVector[i];
        }

        /*
        printf("Rank: %d\n", rank);
        for (int i = 0; i < info->pointsPerProcess; ++i) {
            printf("Distance: %f  ", distVector[i]);
        }
        printf("\n");
        printf("\n");
        */

        // Send all the distances to the master process
        MPI_Send(
                distVector,
                info->pointsPerProcess,
                MPI_DOUBLE,
                master_rank,
                10,
                MPI_COMM_WORLD
        );

        // Receive the median value
        double median;
        MPI_Bcast(
                &median,
                1,
                MPI_DOUBLE,
                master_rank,
                MPI_COMM_WORLD
        );

        // TODO rearrange the elements of the the distVector with the median value received
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

    // Read from the file
    info.points = readFromFile(
            argv[1],
            info.world_rank,
            info.world_size,
            &info.pointsDimension,
            &info.pointsPerProcess
    );

    // Call the function initially with master_rank = 0, rank = world_rank, min_rank = 0, max_rank = world_size - 1
    sort_points(info.world_rank, 0, 0, info.world_size - 1, &info);

    // Finalize the MPI environment.
    MPI_Finalize();

    return 0;
}


//   {
//        volatile int i = 0;
//        char hostname[256];
//        gethostname(hostname, sizeof(hostname));
//        printf("PID %d on %s ready for attach\n", getpid(), hostname);
//        fflush(stdout);
//        while (0 == i)
//            sleep(5);
//    }



//    if (info.world_rank == 0) {
//        for (int i = 0; i < info.pointsPerProcess; ++i) {
//            printf("Rank: %d, Point %d, Coordinates: ", info.world_rank, i);
//            for (int j = 0; j < info.pointsDimension; ++j) {
//                printf("%f ", info.points[i].coordinates[j]);
//            }
//
//            printf("\n");
//        }
//    } else {
//        for (int i = 0; i < info.pointsPerProcess; ++i) {
//            printf("Rank: %d, Point %d, Coordinates: ", info.world_rank, i);
//            for (int j = 0; j < info.pointsDimension; ++j) {
//                printf("%f ", info.points[i].coordinates[j]);
//            }
//
//            printf("\n");
//        }
//    }