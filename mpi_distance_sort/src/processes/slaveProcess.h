#ifndef MPI_DISTANCE_SORT_SLAVEPROCESS_H
#define MPI_DISTANCE_SORT_SLAVEPROCESS_H

#include "../main.h"

void slaveProcess(int master_rank, int min_rank, int max_rank, Info *info, MPI_Comm communicator);
#endif
