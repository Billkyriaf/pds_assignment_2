#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "points.h"
#include "masterProcess.h"

void masterProcess(int master_rank, int min_rank, int max_rank, Info *info, MPI_Comm communicator){
    // Allocate space for the distances from all the other processes
    double *distVector = (double*) malloc(info->pointsPerProcess * (max_rank - min_rank + 1) * sizeof (double*));

    // Allocate space for all the requests
    MPI_Request *requests;
    requests = (MPI_Request *) malloc((max_rank - min_rank) * sizeof (MPI_Request));


    printf("Rank: %d\n", info->world_rank);
    for (int i = 0; i < info->pointsPerProcess; ++i) {
        printf("  Point: %d\n", i);

        for (int j = 0; j < info->pointsDimension; ++j) {
            printf("    %f ", info->points[i].coordinates[j]);
        }

        printf("\n");
    }

    // Start receiving the distances from all the processes. The receives start before the calculation of the
    // distances so if some processes finish before the master process the info send can begin
    for (int i = 0; i < max_rank - min_rank; ++i) {

        // Receive points from the other processes
        MPI_Irecv(
                distVector + (i + 1) * info->pointsPerProcess,
                info->pointsPerProcess,
                MPI_DOUBLE,
                i + 1,
                10,
                communicator,
                requests + i
        );
    }

    // Start calculating the distances from the pivot for all the points
    for (int i = 0; i < info->pointsPerProcess; ++i) {
        // The first points per process instances in the distVector are those of the master process
        distVector[i] = findDistance(&info->pivot, &info->points[i], info->pointsDimension);

        // Update the distance for every point
        info->points[i].distance = distVector[i];
    }


    printf("\n\nRank: %d distance: ", info->world_rank);
    for (int i = 0; i < info->pointsPerProcess; ++i) {
        printf("%f, ", distVector[i]);
    }


    // Wait for all the points to be received
    for (int i = 0; i < max_rank - min_rank; ++i) {
        MPI_Wait(requests + i, MPI_STATUS_IGNORE);
    }

    // Debug prints
    printf("\n\nPivot: ");
    for (int i = 0; i < info->pointsDimension; ++i) {
        printf("%f ", info->pivot.coordinates[i]);
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

    // Create a copy of the dist vector because find median value rearranges the vector
    double *tmpDistVector = (double*) malloc(info->pointsPerProcess * (max_rank - min_rank + 1) * sizeof (double*));

    for (int i = 0; i < info->pointsPerProcess * (max_rank - min_rank + 1); ++i) {
        tmpDistVector[i] = distVector[i];
    }

    // Calculate the median value and send it back to the other processes
    double median = findMedianValue(tmpDistVector, info->pointsPerProcess * info->world_size);
    printf("\n\nMedian: %f\n", median);

    free(tmpDistVector);

    MPI_Bcast(
            &median,
            1,
            MPI_DOUBLE,
            master_rank,
            communicator
    );

    // Clear all the requests
    free(requests);

    // Allocate new memory for the exchanges array ...
    // This array holds the number of points every process wants to send
    int *pointsToSend = (int *)malloc((max_rank - min_rank + 1) * sizeof(int));

    // ... and the requests array
    requests = (MPI_Request *)malloc((max_rank - min_rank) * sizeof (MPI_Request));

    // Start receiving the number of points that every process wants to send
    for (int i = 0; i < max_rank - min_rank; ++i) {
        MPI_Irecv(
                pointsToSend + i + 1,
                1,
                MPI_INT,
                i + 1,
                20,
                communicator,
                requests + i
        );
    }

    // Rearrange the elements of the info->points array with regard to the distVector with the median value received
    int indexI = 0;
    int indexJ = 0;

    // TODO make this run parallel with prefix scan
    while (indexJ < info->pointsPerProcess){
        if (distVector[indexJ] < median) {
            // Swap
            double tmp = distVector[indexJ];
            distVector[indexJ] = distVector[indexI];
            distVector[indexI] = tmp;

            //TODO rearrange the points too
            indexI++;
        }
        indexJ++;
    }

    printf("\n\nDistance after partition: ");
    for (int i = 0; i < info->pointsPerProcess * (max_rank - min_rank + 1); ++i) {
        printf("%f, ", distVector[i]);
    }
    printf("\n");
    printf("\n");

    // This is the number of points the master process wants to pointsToSend
    pointsToSend[0] = info->pointsPerProcess - indexI;

    // Wait for all the points to be received
    for (int k = 0; k < max_rank - min_rank; ++k) {
        MPI_Wait(requests + k, MPI_STATUS_IGNORE);
    }

    for (int i = 0; i < max_rank - min_rank + 1; ++i) {
        printf("Process %d wants to send %d points\n", i, pointsToSend[i]);
    }
    printf("\n");
    printf("\n");

    // Solve he exchanges

    /* This array holds the exchanges for every process the information is stored as follows:
     *
     *   exchanges[i]: This is an array of pointsToSend info for the ith process. There is one for every process.
     *
     *   exchanges[i][0]: This is the number of processes that the ith process needs to communicate with.
     *   exchanges[i][1]: This is the rank of the first process to communicate with.
     *   exchanges[i][2]: This is the number of the points to be exchanged between the ith and the exchanges[i][1] process
     *   .
     *   .
     *   .
     *   exchanges[i][j]: This is the rank of the jth process to communicate with.
     *   exchanges[i][j + 1]: This is the number of the points to be exchanged between the ith and the exchanges[i][j] process
     *
     *
     * Note that it is not mandatory for every process to pointsToSend points with all the other processes. This is actually
     * something to avoid. This is why the exchanges[i][0] shows the number of exchanges for every process.
     */
    int **exchanges = (int**)malloc((max_rank - min_rank + 1) * sizeof (int *));
    for (int k = 0; k <= max_rank - min_rank; ++k) {
        exchanges[k] = (int *)calloc((max_rank - min_rank) * 2 + 1, sizeof (int));
    }

    /* The pointsToSend vector is split in half and the bottom half wants to send data to the top half and vice versa
     * This means that in an MPI world with 8 processes if the process 0 wants to send to process 4 100 points the
     * process 4 will additionally send process 0 100 points.
     *
     * processes are matched in pairs initially
     *
     * e.g. For 8 processes
     *                              ┌───────────────┐
     *                      ┌───────┼───────┐       │
     *      processes:  0   1   2   3   4   5   6   7
     *                  └───────┼───────┘       │
     *                          └───────────────┘
     *
     *
     *
     *      If the process 4 can not take all the points from the process 1 process 1 will send also to the process 5
     *      etc. until all the processes are covered.
     */

    int indexTop = (max_rank - min_rank + 1) / 2;  // This is the rank of the first top process

    // For the bottom half of the processes...
    for (int k = 0; k < (max_rank - min_rank + 1) / 2; ++k) {
        while(pointsToSend[k] > 0){


            // There is a chance that the top process doesn't want to take any points, so we skip it
            while (pointsToSend[indexTop] == 0){
                indexTop++;
            }

            if(pointsToSend[k] == pointsToSend[indexTop]) {
                // Those are the indexes of the exchanges[k (+ indexTop)] that the new exchanges will be written
                int bottomInd = exchanges[k][0] * 2 + 1;
                int topInd = exchanges[indexTop][0] * 2 + 1;

                exchanges[k][0]++;  // add one more exchange to the list
                exchanges[indexTop][0]++;  // add one more exchange to the list


                // Fill the exchanges array
                exchanges[k][bottomInd] = indexTop;  // The rank of the process to send the points
                exchanges[k][bottomInd + 1] = pointsToSend[k];  // The number of points to send

                exchanges[indexTop][topInd] = k;  // The rank of the process to send the points
                exchanges[indexTop][topInd + 1] = pointsToSend[k];  // The number of points to send


                // Update the remaining points
                pointsToSend[k] = 0;  // This process has no more points to send
                pointsToSend[indexTop] = 0;  // This process has no more points to send
                indexTop++;  // incrementing the indexTop means that a process is covered

            } else if (pointsToSend[k] < pointsToSend[indexTop]){
                // Those are the indexes of the exchanges[k (+ indexTop)] that the new exchanges will be written
                int bottomInd = exchanges[k][0] * 2 + 1;
                int topInd = exchanges[indexTop][0] * 2 + 1;

                exchanges[k][0]++;  // add one more exchange to the list
                exchanges[indexTop][0]++;  // add one more exchange to the list

                // Fill the exchanges array
                exchanges[k][bottomInd] = indexTop;  // The rank of the process to send the points
                exchanges[k][bottomInd + 1] = pointsToSend[k];  // The number of points to send

                exchanges[indexTop][topInd] = k;  // The rank of the process to send the points
                exchanges[indexTop][topInd + 1] = pointsToSend[k];  // The number of points to send

                // Update the remaining points
                pointsToSend[indexTop] -= pointsToSend[k];  // Some points of this process are covered
                pointsToSend[k] = 0;  // This process has no more points to send

            } else if(pointsToSend[k] > pointsToSend[indexTop]){
                // Those are the indexes of the exchanges[k (+ indexTop)] that the new exchanges will be written
                int bottomInd = exchanges[k][0] * 2 + 1;
                int topInd = exchanges[indexTop][0] * 2 + 1;

                exchanges[k][0]++;  // add one more exchange to the list
                exchanges[indexTop][0]++;  // add one more exchange to the list

                // Fill the exchanges array
                exchanges[k][bottomInd] =indexTop;  // The rank of the process to send the points
                exchanges[k][bottomInd + 1] = pointsToSend[indexTop];  // The number of points to send

                exchanges[indexTop][topInd] = k;  // The rank of the process to send the points
                exchanges[indexTop][topInd + 1] = pointsToSend[indexTop];  // The number of points to send

                // Update the remaining points
                pointsToSend[k] -= pointsToSend[indexTop];  // Some points of this process are covered
                pointsToSend[indexTop] = 0;  // This process has no more points to send
                indexTop++;  // incrementing the indexTop means that a process is covered
            }
        }
    }

    for (int k = 0; k <= max_rank - min_rank; ++k) {
        printf("Rank %d wants to send ", k);
        for (int l = 1; l <= exchanges[k][0] * 2; l += 2) {
            printf("Rank %d, %d element(s), ", exchanges[k][l], exchanges[k][l + 1]);
        }
        printf("\n");
    }

}