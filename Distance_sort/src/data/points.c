#include "points.h"
#include <math.h>

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
float findDistance(Point *pivot, Point *target, int dimension){
    float distance = 0;

    // add all the pairs that give the euclidean distance
    for (int coord = 0; coord < dimension; ++coord) {
        distance += powf((pivot->coordinates[coord] - target->coordinates[coord]), 2);
    }

    // set the target point distance
    target->distance = distance;

    // return a type cast to float because there is no need for double precision
    return distance;
}