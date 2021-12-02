#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#include "points.h"
#include "read_points.h"


int main(int argc, char **argv) {
    printf("%s\n", argv[1]);
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
        printf("Master thread sending...\n");

        int *arr;
        arr = (int *)malloc(world_size * sizeof(int));

        for (int i = 0; i < world_size; i++) {
            arr[i] = i;
            MPI_Isend((arr + i), 1, MPI_INT, i, 10, MPI_COMM_WORLD, &mpiRequest);
        }

        MPI_Wait(&mpiRequest, &mpiStatus);

        free(arr);
        printf("Finished\n");

    } else {
        int a = 0;

        int receive[1];

        MPI_Recv(receive, 1, MPI_INT, 0, 11, MPI_COMM_WORLD, &mpiStatus);

        printf("Rank %d received: %d, a = %d\n", world_rank, receive[0], a);
    }

    // Finalize the MPI environment.
    MPI_Finalize();

    return 0;
}
