#ifndef DISTANCE_SORT_READ_POINTS_H
#define DISTANCE_SORT_READ_POINTS_H

#include "points.h"

Point *readFromFile(char *fileName, int rank, int worldSize, int *pointsDimension, int *pointsPerProcess);


#endif
