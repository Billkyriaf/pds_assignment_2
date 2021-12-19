#ifndef MPI_DISTANCE_SORT_MASTERPROCESS_H
#define MPI_DISTANCE_SORT_MASTERPROCESS_H

#include "../main.h"

void masterProcess(int master_rank, int min_rank, int max_rank, Info *info, MPI_Comm communicator);
#endif
