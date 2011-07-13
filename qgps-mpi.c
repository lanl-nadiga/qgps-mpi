#include "qgps-mpi.h"

const int QGPS_NX = 64, QGPS_NY = 64;

int QGPS_COMM_WORLD;
int  qgps_current_task, qgps_master_task = 0, qgps_number_tasks;

qgps_block_t *qgps_blocks = 0;

int qgps_initialize(int argc, char **argv) {
        MPI_Init(&argc, &argv);
        MPI_Comm_dup(MPI_COMM_WORLD, &QGPS_COMM_WORLD);

        MPI_Comm_rank(QGPS_COMM_WORLD, &qgps_current_task);
        MPI_Comm_size(QGPS_COMM_WORLD, &qgps_number_tasks);

        qgps_blocks = calloc(qgps_number_tasks, sizeof(qgps_block_t));
        if(!qgps_blocks)
                return 1;

        return 0;
}

int qgps_cleanup() {
        MPI_Finalize();

        return 0;
}
