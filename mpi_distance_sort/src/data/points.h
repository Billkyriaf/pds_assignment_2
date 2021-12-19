#ifndef DISTANCE_SORT_POINTS_H
#define DISTANCE_SORT_POINTS_H

/**
 * This struct represents a point.
 */
typedef struct point{
    double *coordinates;  // vector of the size of the space dimension for all the coordinates
    double distance;  // distance from the selected pivot
} Point;

double findDistance(Point *pivot, Point *target, int dimension);
double findMedianValue(double *distanceVector, int numberOfPoints);

#endif
