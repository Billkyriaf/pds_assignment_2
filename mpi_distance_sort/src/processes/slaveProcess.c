#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "points.h"
#include "masterProcess.h"

#include "slaveProcess.h"


void slaveProcess(int master_rank, int min_rank, int max_rank, Info *info, MPI_Comm communicator){
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

    // Memory allocation for the distance vector
    double *distVector = (double*) malloc(info->pointsPerProcess * sizeof (double));


    // Start calculating the distances from the pivot for all the points
    for (int i = 0; i < info->pointsPerProcess; ++i) {
        distVector[i] = findDistance(&info->pivot, &info->points[i], info->pointsDimension);
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
            communicator
    );

    // Receive the median value
    double median;
    MPI_Bcast(
            &median,
            1,
            MPI_DOUBLE,
            master_rank,
            communicator
    );

    // Rearrange the elements of the distVector with the median value received
    int i = 0;
    int j = 0;


    // TODO make this run parallel with prefix scan
    while (j < info->pointsPerProcess){
        if (distVector[j] < median) {
            // Swap
            double tmp = distVector[j];
            distVector[j] = distVector[i];
            distVector[i] = tmp;

            i++;
        }
        j++;
    }

    // We need to exchange the elements from i (inclusive) to the end (size - i elements) if this is a bottom process
    // Or if this is a top process the elements from start to i (exclusive). This is i elements
    // First step is to send to the main process the number of elements we want to exchange

    int pointsToSend;  // Number of points that need to be exchanged.

    // Bottom process
    if (info->world_rank < (max_rank - min_rank + 1) / 2) {
        pointsToSend = info->pointsPerProcess - i;

        // Send this information to the master process
        MPI_Send(&pointsToSend, 1, MPI_INT, master_rank, 20, communicator);

    } else { // Top process
        pointsToSend = i;

        // Send this information to the master process
        MPI_Send(&pointsToSend, 1, MPI_INT, master_rank, 20, communicator);
    }


}