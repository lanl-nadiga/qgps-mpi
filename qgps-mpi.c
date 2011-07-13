#include "qgps-mpi.h"
#include <stdlib.h>

const int QGPS_NX = 64, QGPS_NY = 64;

const MPI_Comm QGPS_COMM_WORLD;
int  qgps_current_task, qgps_master_task = 0, qgps_number_tasks;

qgps_block_t *qgps_blocks = NULL;
qgps_block_t *qgps_transpose_blocks = NULL;

int qgps_initialize_mpi(int argc, char **argv);
int qgps_initialize_blocks();
int qgps_cleanup_mpi();

int qgps_initialize(int argc, char **argv) {
        if (qgps_initialize_mpi(argc, argv))
                return 1;
        if (qgps_initialize_blocks())
                return 1;

        return 0;
}

int qgps_cleanup() {
        if (qgps_cleanup_mpi())
                return 1;

        return 0;
}

int qgps_initialize_mpi(int argc, char **argv) {
        MPI_Init(&argc, &argv);
        MPI_Comm_dup(MPI_COMM_WORLD, (MPI_Comm*)&QGPS_COMM_WORLD);

        MPI_Comm_rank(QGPS_COMM_WORLD, &qgps_current_task);
        MPI_Comm_size(QGPS_COMM_WORLD, &qgps_number_tasks);

        return 0;
}
int qgps_initialize_blocks(int argc, char **argv) {
        qgps_blocks = calloc(qgps_number_tasks, sizeof(qgps_block_t));
        if(!qgps_blocks)
                return 1;

        qgps_transpose_blocks = calloc(qgps_number_tasks, sizeof(qgps_block_t));
        if(!qgps_transpose_blocks)
                return 1;

        return 0;
}
int qgps_cleanup_mpi() {
        MPI_Finalize();

        return 0;
}
