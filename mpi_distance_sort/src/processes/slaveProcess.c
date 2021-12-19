#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "points.h"
#include "masterProcess.h"

#include "slaveProcess.h"


void slaveProcess(int master_rank, int min_rank, int max_rank, Info *info, MPI_Comm communicator){



//    printf("\n\nRank: %d\n", info->world_rank);
//    for (int i = 0; i < info->pointsPerProcess; ++i) {
//        printf("  Point: %d\n", i);
//        for (int j = 0; j < info->pointsDimension; ++j) {
//            printf("    %f ", info->points[i].coordinates[j]);
//        }
//        printf("\n");
//    }



    // Memory allocation for the distance vector
    double *distVector = (double*) malloc(info->pointsPerProcess * sizeof (double));  // MEMORY


    // Start calculating the distances from the pivot for all the points
    for (int i = 0; i < info->pointsPerProcess; ++i) {
        distVector[i] = findDistance(&info->pivot, &info->points[i], info->pointsDimension);
        info->points[i].distance = distVector[i];
    }



//    printf("Rank: %d\n", rank);
//    for (int indexI = 0; indexI < info->pointsPerProcess; ++indexI) {
//        printf("Distance: %f  ", distVector[indexI]);
//    }
//    printf("\n");
//    printf("\n");



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

    // Update the median of the info struct (this is used for testing later)
    info->median = median;


    // Rearrange the elements of the info->points array with the median value received
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

    // We need to exchange the elements from indexI (inclusive) to the end (size - indexI elements) if this is a bottom process
    // Or if this is a top process the elements from start to indexI (exclusive). This is indexI elements
    // First step is to send to the main process the number of elements we want to exchange
    MPI_Request request;

    int pointsToSend;  // Number of points that need to be exchanged.
    double *sendVector;  // Stores the points that the process wants to exchange

    // If this is a bottom process it wants to exchange the larger points
    if (info->world_rank < (max_rank - min_rank + 1) / 2) {
        pointsToSend = info->pointsPerProcess - indexI;

        // Send this information to the master process
        MPI_Isend(
                &pointsToSend,
                1,
                MPI_INT,
                master_rank,
                20,
                communicator,
                &request
        );

        // Allocate memory for the sendVector
        sendVector = (double *) malloc(pointsToSend * info->pointsDimension * sizeof(double));  // MEMORY

        // Fill the sendVector with te coordinates of the points to send
        for (int k = indexI; k < indexI + pointsToSend; ++k) {
            for (int l = 0; l < info->pointsDimension; ++l) {
                sendVector[(k - indexI) * info->pointsDimension + l] = info->points[k].coordinates[l];
            }
        }


    } else { // If this is a top process it wants to exchange the smallest points
        pointsToSend = indexI;

        // Send this information to the master process
        MPI_Isend(
                &pointsToSend,
                1,
                MPI_INT,
                master_rank,
                20,
                communicator,
                &request
        );

        // Allocate memory for the sendVector
        sendVector = (double *) malloc(pointsToSend * info->pointsDimension * sizeof(double));  // MEMORY

        // Fill the sendVector with te coordinates of the points to send
        for (int k = 0; k < indexI; ++k) {
            for (int l = 0; l < info->pointsDimension; ++l) {
                sendVector[k * info->pointsDimension + l] = info->points[k].coordinates[l];
            }
        }
    }

    // There is no need to wait for the "send" above because we enter a "receive" that the master process will initiate
    // only if the above "send" is completed

    // Receive the exchanges vector
    MPI_Status status;

    // Probe the incoming message to get the size of the incoming vector
    MPI_Probe(master_rank, 30, communicator, &status);

    int nItems;  // The number of items in the exchange vector

    // Get the number of the incoming elements
    MPI_Get_count(&status, MPI_INT, &nItems);

    /* This vector holds the exchanges for this process the information is stored as follows:
     *
     *
     *   exchanges[0]: This is the number of processes that the ith process needs to communicate with.
     *   exchanges[1]: This is the rank of the first process to communicate with.
     *   exchanges[2]: This is the number of the points to be exchanged between the ith and the exchanges[1] process
     *   .
     *   .
     *   .
     *   exchanges[indexJ]: This is the rank of the jth process to communicate with.
     *   exchanges[indexJ + 1]: This is the number of the points to be exchanged between the ith and the exchanges[indexJ] process
     *
     */
    int *exchanges = (int *) malloc(nItems * sizeof(int));  // MEMORY

    // Receive the exchange solution from the master process
    MPI_Recv(
            exchanges,
            nItems,
            MPI_INT,
            master_rank,
            30,
            communicator,
            MPI_STATUS_IGNORE
    );


    // Allocate memory to keep track of all the requests
    MPI_Request *requests = (MPI_Request *) malloc(exchanges[0] * 2 * sizeof (MPI_Request));  // MEMORY


    // The recVector temporarily holds the received coordinates. This and the sendVector are necessary in order to
    // send and receive data at the same time and avoid race conditions.
    double *recVector = (double *) malloc(pointsToSend * info->pointsDimension * sizeof(double));  // MEMORY

    // The offset is used to offset the recVector and sendVector so that they both can send and receive
    // from multiple processes at the same time
    int offset = 0;

    // Start the points exchange process
    for (int k = 0; k < exchanges[0]; ++k) {
        // Non-blocking receive of the points
        MPI_Irecv(
                recVector + offset,
                exchanges[k * 2 + 2] * info->pointsDimension,
                MPI_DOUBLE,
                exchanges[k * 2 + 1],
                40,
                communicator,
                requests + 2 * k
        );

        // Non-blocking send of the points
        MPI_Isend(
                sendVector + offset,
                exchanges[k * 2 + 2] * info->pointsDimension,
                MPI_DOUBLE,
                exchanges[k * 2 + 1],
                40,
                communicator,
                requests + 2 * k + 1
        );

        // Increment the offset by the number of coordinates written and read
        offset += exchanges[k * 2 + 2] * info->pointsDimension;
    }

    // Wait for all the requests to complete
    for (int k = 0; k < exchanges[0] * 2; ++k) {
        MPI_Wait(requests + k, MPI_STATUS_IGNORE);
    }


    // Update the points array with the new points received

    // If this is a bottom process update the larger points
    if (info->world_rank < (max_rank - min_rank + 1) / 2) {
        for (int k = indexI; k < indexI + pointsToSend; ++k) {
            for (int l = 0; l < info->pointsDimension; ++l) {
                info->points[k].coordinates[l] = recVector[(k - indexI) * info->pointsDimension + l];
            }
        }

    } else { // If this is a top process update the lower points
        for (int k = 0; k < indexI; ++k) {
            for (int l = 0; l < info->pointsDimension; ++l) {
                info->points[k].coordinates[l] = recVector[k * info->pointsDimension + l];
            }
        }
    }



//    printf("\n\nRank: %d after\n", info->world_rank);
//    for (int i = 0; i < info->pointsPerProcess; ++i) {
//        printf("  Point: %d\n", i);
//        for (int j = 0; j < info->pointsDimension; ++j) {
//            printf("    %f ", info->points[i].coordinates[j]);
//        }
//        printf("\n");
//    }



    // Free memory
    free(distVector);
    free(sendVector);
    free(exchanges);
    free(requests);
    free(recVector);

}