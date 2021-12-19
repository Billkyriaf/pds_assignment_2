#ifndef MPI_DISTANCE_SORT_MAIN_H
#define MPI_DISTANCE_SORT_MAIN_H
/**
 * This structure passes constant information between the recursive call of the sort_points function, that the arguments
 * of the function are kept as few as possible for simplicity
 */
typedef struct info{
    int world_size;  // The size of the world (number of processes)
    int world_rank;  // The rank of the current process

    Point pivot;  // The initial pivot point selected by the master process
    Point *points;  // All the points for each process
    int pointsPerProcess;  // The total number of points each process has
    int pointsDimension;  // The dimension of the space for the points (number of coordinates)
} Info;


#endif
