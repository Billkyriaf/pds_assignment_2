#include <pthread.h>

#include "threads.h"

/**
 * The runnable function for parallel computation of he distances
 * @param args The arguments struct for every thread
 */
void *distanceRunnable(void *args){

    // Type cast the arguments
    PthreadArgs *arguments = (PthreadArgs *) args;

    // This is the reference to the processes info struct. This variable is only here to make the code more readable.
    Info *info = (*arguments).info;

    // Start calculating the distances from the pivot for all the points
    for (int i = (*arguments).startIndex; i <= (*arguments).endIndex; ++i) {

        // Calculate the distance
        (*arguments).distPtr[i] = findDistance(&info->pivot, &info->points[i], info->pointsDimension);

        // Update the info struct
        info->points[i].distance = (*arguments).distPtr[i];
    }

    pthread_exit(NULL);
}