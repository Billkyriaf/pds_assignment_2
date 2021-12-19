#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#include "data/points.h"
#include "data/read_points.h"
#include "processes/masterProcess.h"
#include "processes/slaveProcess.h"
#include "main.h"


void sort_points(int master_rank, int min_rank, int max_rank, Info *info, MPI_Comm communicator){
    // Recursion finish statement
    if (min_rank - max_rank == 0){
        return;
    }

    if  (info->world_rank == master_rank){
        masterProcess(master_rank, min_rank, max_rank, info, communicator);
    } else {
        slaveProcess(master_rank, min_rank, max_rank, info, communicator);
    }



    // TODO call the function recursively with new masters for half of the processes
}

int main(int argc, char **argv) {

    // check for the command line arguments
    if (argc != 2){
        printf("Wrong number of arguments exiting...\n");
        return -1;
    }

    // Initialize the MPI communication
    MPI_Init(&argc, &argv);

    Info info;

    // Get the number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &info.world_size);

    // Get the rank of the process
    MPI_Comm_rank(MPI_COMM_WORLD, &info.world_rank);

    // Read from the file
    info.points = readFromFile(
            argv[1],
            info.world_rank,
            info.world_size,
            &info.pointsDimension,
            &info.pointsPerProcess
    );

    // Call the function initially with master_rank = 0, rank = world_rank, min_rank = 0, max_rank = world_size - 1
    sort_points(info.world_rank, 0, 0, info.world_size - 1, &info);

    // Finalize the MPI environment.
    MPI_Finalize();

    return 0;
}


//   {
//        volatile int i = 0;
//        char hostname[256];
//        gethostname(hostname, sizeof(hostname));
//        printf("PID %d on %s ready for attach\n", getpid(), hostname);
//        fflush(stdout);
//        while (0 == i)
//            sleep(5);
//    }



//    if (info.world_rank == 0) {
//        for (int i = 0; i < info.pointsPerProcess; ++i) {
//            printf("Rank: %d, Point %d, Coordinates: ", info.world_rank, i);
//            for (int j = 0; j < info.pointsDimension; ++j) {
//                printf("%f ", info.points[i].coordinates[j]);
//            }
//
//            printf("\n");
//        }
//    } else {
//        for (int i = 0; i < info.pointsPerProcess; ++i) {
//            printf("Rank: %d, Point %d, Coordinates: ", info.world_rank, i);
//            for (int j = 0; j < info.pointsDimension; ++j) {
//                printf("%f ", info.points[i].coordinates[j]);
//            }
//
//            printf("\n");
//        }
//    }