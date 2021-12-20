#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <unistd.h>
#include "points.h"
#include "masterProcess.h"

void masterProcess(int master_rank, int min_rank, int max_rank, Info *info, MPI_Comm communicator){
    // Allocate space for the distances from all the other processes
    double *distVector = (double*) malloc(info->pointsPerProcess * (max_rank - min_rank + 1) * sizeof (double));  // MEMORY

    // Allocate space for all the requests objects
    MPI_Request *requests;
    requests = (MPI_Request *) malloc((max_rank - min_rank) * sizeof (MPI_Request));  // MEMORY



//    printf("\n\nRank: %d\n", info->world_rank);
//    for (int i = 0; i < info->pointsPerProcess; ++i) {
//        printf("  Point: %d\n", i);
//
//        for (int j = 0; j < info->pointsDimension; ++j) {
//            printf("    %f ", info->points[i].coordinates[j]);
//        }
//
//        printf("\n");
//    }



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



//    printf("\n\nRank: %d distance: ", info->world_rank);
//    for (int i = 0; i < info->pointsPerProcess; ++i) {
//        printf("%f, ", distVector[i]);
//    }
//    printf("\n");


    // Wait for all the distances to be received
    for (int i = 0; i < max_rank - min_rank; ++i) {
        MPI_Wait(requests + i, MPI_STATUS_IGNORE);
    }



    // Debug prints
//    printf("\n\nPivot: ");
//    for (int i = 0; i < info->pointsDimension; ++i) {
//        printf("%f ", info->pivot.coordinates[i]);
//    }
//    printf("\n");
//    printf("\n");
//
//    for (int i = 0; i <= max_rank; ++i) {
//        printf("Distance table Rank %d: ", i);
//        for (int j = 0; j < info->pointsPerProcess; ++j) {
//            printf("%f ", distVector[i * info->pointsPerProcess + j]);
//        }
//        printf("\n");
//    }



    // Calculate the median value and send it back to the other processes of the group
    double median = findMedianValue(distVector, info->pointsPerProcess * info->world_size);
    info->median = median;



//    printf("\n\nMedian: %f\n", median);



    // Broadcast the median to all the processes of the group
    MPI_Bcast(
            &median,
            1,
            MPI_DOUBLE,
            master_rank,
            communicator
    );

    // Clear all the requests
    free(requests);  // MEMORY free

    // Allocate new memory for the exchanges array ...
    // This array holds the number of points every process wants to send
    int *pointsToSend = (int *)malloc((max_rank - min_rank + 1) * sizeof(int));  // MEMORY

    // ... and the requests array
    requests = (MPI_Request *)malloc((max_rank - min_rank) * sizeof (MPI_Request));  // MEMORY

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

    // Rearrange the elements of the info->points array with regard to the distances with the median value received
    int indexI = 0;
    int indexJ = 0;

    // TODO make this run parallel with prefix scan
    while (indexJ < info->pointsPerProcess){
        if (info->points[indexJ].distance < median) {

            Point pnt = info->points[indexJ];
            info->points[indexJ] = info->points[indexI];
            info->points[indexI] = pnt;

            indexI++;
        }
        indexJ++;
    }

    // This is the number of points the master process wants to send
    pointsToSend[0] = info->pointsPerProcess - indexI;



//    printf("\n\nDistance after partition: ");
//    for (int i = 0; i < info->pointsPerProcess * (max_rank - min_rank + 1); ++i) {
//        printf("%f, ", distVector[i]);
//    }
//    printf("\n");
//    printf("\n");



    // Create a "send" and a "receive" buffer. Those buffers may seem like a waist of memory but in reality they are the
    // reason that the "send" and "receive" process can happen simultaneously

    // Stores temporarily the points that need to be sent
    double *sendVector = (double *) malloc(pointsToSend[0] * info->pointsDimension * sizeof(double));  // MEMORY

    // Stores temporarily the points that are received until they are permanently stored
    double *recVector = (double *) malloc(pointsToSend[0] * info->pointsDimension * sizeof(double));  // MEMORY

    // Because this is the master process and the master process will always be a bottom process the points that the
    // process needs to send will always be on the right after the rearrangement above
    for (int k = indexI; k < indexI + pointsToSend[0]; ++k) {

        // Copy all the coordinates to the sendVector
        for (int l = 0; l < info->pointsDimension; ++l) {
            sendVector[(k - indexI) * info->pointsDimension + l] = info->points[k].coordinates[l];
        }
    }

    // From this point and forward the number of points that every process wants to exchange is needed so we wait
    for (int k = 0; k < max_rank - min_rank; ++k) {
        MPI_Wait(requests + k, MPI_STATUS_IGNORE);
    }



//    for (int i = 0; i < max_rank - min_rank + 1; ++i) {
//        printf("Process %d wants to send %d points\n", i, pointsToSend[i]);
//    }
//    printf("\n");
//    printf("\n");



    // Solve he exchanges problem

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
    int **exchanges = (int**)calloc((max_rank - min_rank + 1), sizeof (int *));  // MEMORY
    for (int k = 0; k <= max_rank - min_rank; ++k) {
        exchanges[k] = (int *)calloc((max_rank - min_rank) * 2 + 1, sizeof (int));  // MEMORY
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



//    printf("\n");
//    printf("\n");
//    for (int k = 0; k <= max_rank - min_rank; ++k) {
//        printf("Rank %d wants to send ", k);
//        for (int l = 1; l <= exchanges[k][0] * 2; l += 2) {
//            printf("Rank %d, %d element(s), ", exchanges[k][l], exchanges[k][l + 1]);
//        }
//        printf("\n");
//    }



    // Clear all the requests
    free(requests);  // MEMORY free

    // ... and allocate the requests array again for the new requests
    requests = (MPI_Request *)malloc((max_rank - min_rank) * sizeof (MPI_Request));  // MEMORY

    // Send the solution to all the processes
    for (int i = 0; i < max_rank - min_rank; ++i) {
        MPI_Isend(
                exchanges[i + 1],
                exchanges[i + 1][0] * 2 + 1,
                MPI_INT,
                i + 1,
                30,
                communicator,
                &requests[i]
        );
    }

    // Wait for the requests to finish
    for (int k = 0; k < max_rank - min_rank; ++k) {
        MPI_Wait(&requests[k], MPI_STATUS_IGNORE);
    }

    // Clear all the requests
    free(requests); // MEMORY free


    // ... and allocate the requests array again for the new requests
    requests = (MPI_Request *) malloc(exchanges[0][0] * 2 * sizeof (MPI_Request)); // MEMORY

    // The offset is used to offset the recVector and sendVector so that they both can send and receive
    // from multiple processes at the same time
    int offset = 0;

    // Start the points exchange process
    for (int k = 0; k < exchanges[0][0]; ++k) {
        // Non-blocking receive of the points
        MPI_Irecv(
                recVector + offset,
                exchanges[0][k * 2 + 2] * info->pointsDimension,
                MPI_DOUBLE,
                exchanges[0][k * 2 + 1],
                40, communicator,
                requests + 2 * k
        );

        // Non-blocking send of the points
        MPI_Isend(
                sendVector + offset,
                exchanges[0][k * 2 + 2] * info->pointsDimension,
                MPI_DOUBLE,
                exchanges[0][k * 2 + 1],
                40,
                communicator,
                requests + 2 * k + 1
        );

        // Increment the offset by the number of coordinates written and read
        offset += exchanges[0][k * 2 + 2] * info->pointsDimension;
    }


    if (exchanges[0] > 0){
        // Wait for all the requests to complete
        for (int k = 0; k < exchanges[0][0] * 2; ++k) {
            MPI_Wait(requests + k, MPI_STATUS_IGNORE);
        }
    }


    // Permanently store the points received
    for (int k = indexI; k < info->pointsPerProcess; ++k) {
        for (int l = 0; l < info->pointsDimension; ++l) {
            info->points[k].coordinates[l] = recVector[(k - indexI) * info->pointsDimension + l];
        }
    }



//    printf("\n\nRank: %d after\n", info->initial_rank);
//    for (int i = 0; i < info->pointsPerProcess; ++i) {
//        printf("  Point: %d, Distance: %.10f\n", i, info->points[i].distance);
//        for (int j = 0; j < info->pointsDimension; ++j) {
//            printf("    %f ", info->points[i].coordinates[j]);
//        }
//        printf("\n");
//    }



    // free memory
    free(distVector);
    free(requests);
    free(pointsToSend);
    free(sendVector);
    free(recVector);
    for (int i = 0; i < max_rank - min_rank + 1; ++i) {
        free(exchanges[i]);
    }
    free(exchanges);
}