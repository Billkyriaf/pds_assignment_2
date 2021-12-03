#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#include "points.h"
#include "read_points.h"


int main(int argc, char **argv) {
    // Initialize the MPI communication
    MPI_Init(&argc, &argv);

    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // Read the initial points
    Point *points = NULL;  // All the points for each process
    int pointsPerProcess;  // The total number of points each process has
    int pointsDimension;  // The dimension of the space for the points (number of coordinates)

    // Read from the file
    readFromFile(argv[1], points, world_rank, world_size, &pointsDimension, &pointsPerProcess);


    MPI_Request mpiRequest;
    MPI_Status mpiStatus;

    if (world_rank == 0) {
        for (int i = 0; i < pointsPerProcess; ++i) {
            printf("Rank: %d, Point %d coordinates: ",world_rank, i);
            for (int j = 0; j < pointsDimension; ++j) {
                printf("%f ", points[i].coordinates[j]);
            }

            printf("\n");
        }
    } else {
        for (int i = 0; i < pointsPerProcess; ++i) {
            printf("Rank: %d, Point %d coordinates: ",world_rank, i);
            for (int j = 0; j < pointsDimension; ++j) {
                printf("%f ", points[i].coordinates[j]);
            }

            printf("\n");
        }
    }

    // Finalize the MPI environment.
    MPI_Finalize();

    return 0;
}
