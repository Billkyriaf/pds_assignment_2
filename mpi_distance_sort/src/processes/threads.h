#ifndef MPI_DISTANCE_SORT_THREADS_H
#define MPI_DISTANCE_SORT_THREADS_H

#include "../data/structs.h"

/**
 * Arguments for the runnable function that calculates the distances in parallel to save time
 */
typedef struct pthreadArgs {
    Info *info;  // Reference to the info struct of the process
    double *distPtr;  // Reference to the distance vector of the process

    int startIndex;  // The starting index of the distance vector for each thread
    int endIndex;  // The ending index of the distance vector for each thread
} PthreadArgs;

void *distanceRunnable(void *args);

#endif
