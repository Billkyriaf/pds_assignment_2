#ifndef DISTANCE_SORT_POINTS_H
#define DISTANCE_SORT_POINTS_H

typedef struct point{
    float *coordinates;
    float distance;
} Point;

float findDistance(Point *pivot, Point *target, int dimension);

#endif
