#include "points.h"
#include <math.h>
#include <stdlib.h>

/**
 * Gets two points of the struct Point and returns their Euclidean distance (almost read bellow). The return type is
 * float because there is no need for double precision.
 * Since we do not care for the exact distance but instead we want to compare the values for each point and since the
 * square root is a strictly increasing function we do not need to calculate the square root.
 *
 * @param pivot First point
 * @param target Second point
 * @param dimension The the dimension of the space
 * @return Euclidean distance (float)
 */
double findDistance(Point *pivot, Point *target, int dimension){
    double distance = 0;

    // add all the pairs that give the euclidean distance
    for (int coord = 0; coord < dimension; ++coord) {
        distance += pow((pivot->coordinates[coord] - target->coordinates[coord]), 2);
    }

    // set the target point distance
    target->distance = distance;

    return distance;
}

/**
 * Swaps two doubles in the memory
 * @param a First point to swap
 * @param b Second point to swap
 */
void swap(double* a, double* b)
{
    double temp = *a;
    *a = *b;
    *b = temp;
}

/**
 * Puts the pivot to the correct positions and returns the index. The pivot is previously swapped to the distVector[end]
 * @param distVector The vector of all the distances
 * @param start The starting index of the partition
 * @param end The ending index of the partition
 * @return The index of the pivot
 */
int Partition(double *distVector, int start, int end) {
    double pivotVal = distVector[end];
    int i = start, j = start;

    while (j < end) {
        if (distVector[j] < pivotVal) {
            swap(&distVector[i], &distVector[j]);
            i++;
        }
        j++;
    }
    swap(&distVector[i], &distVector[end]);
    return i;
}

/**
 * Picks a random point from the partition moves it to the and of the partition and calls Partition in order to find the
 * correct index for the selected value
 * @param distVector The vector of all the distances
 * @param start The starting index of the partition
 * @param end The ending index of the partition
 * @return The index of the pivot
 */
int randomPartition(double *distVector, int start, int end) {
    int n = end - start + 1;
    int pivot = rand() % n;
    swap(&distVector[start + pivot], &distVector[end]);
    return Partition(distVector, start, end);
}


/**
 * This function is called recursively until the median(s) is found.
 * @param distVector The vector of all the distances
 * @param start The starting index of the partition
 * @param end The ending index of the partition
 * @param k The middle index of the distVector
 * @param a The value of the median(s)
 * @param b The value of the median(s)
 */
void MedianUtil(double *distVector, int start, int end, int k, double *a, double *b) {

    // if start < end
    if (start <= end) {

        // Find the partition index
        int partitionIndex = randomPartition(distVector, start, end);

        // If partition index = k, then we found the median of odd number element in distVector
        if (partitionIndex == k) {
            *b = distVector[partitionIndex];
            if (*a != -1){
                return;
            }

        } else if (partitionIndex == k - 1) {
            // If index = k - 1, then we get a & b as middle element of distVector
            *a = distVector[partitionIndex];
            if (*b != -1) {
                return;
            }
        }

        // If partitionIndex >= k then find the index in first half of the distVector
        if (partitionIndex >= k) {
            return MedianUtil(distVector, start, partitionIndex - 1, k, a, b);

        } else {
            // If partitionIndex <= k then find the index in second half of the distVector
            return MedianUtil(distVector,partitionIndex + 1, end, k, a, b);
        }
    }
}

/**
 * Finds the median value of the vector using quick select algorith
 * @param distanceVector The vector of all the distances
 * @param numberOfPoints The size of the vector
 * @return The median value
 */
double findMedianValue(double *distanceVector, int numberOfPoints){
    double ans, a = -1.0, b = -1.0;

    // If numberOfPoints is odd
    if (numberOfPoints % 2 == 1) {
        MedianUtil(distanceVector, 0, numberOfPoints - 1,numberOfPoints / 2, &a, &b);
        ans = b;

    } else { // If numberOfPoints is even
        MedianUtil(distanceVector, 0, numberOfPoints - 1,numberOfPoints / 2, &a, &b);
        ans = (a + b) / 2;
    }

    return ans;
}