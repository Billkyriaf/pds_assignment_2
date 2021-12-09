#ifndef DISTANCE_SORT_POINTS_H
#define DISTANCE_SORT_POINTS_H

typedef struct point{
    double *coordinates;
    double distance;
} Point;

double findDistance(Point *pivot, Point *target, int dimension);
double findMedianValue(double *distanceVector, int numberOfPoints);

#endif
