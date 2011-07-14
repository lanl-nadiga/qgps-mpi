#include "qgps-mpi.h"
#include "step.h"
#include <stdlib.h>

const int QGPS_NX = 64, QGPS_NY = 64;

MPI_Comm QGPS_COMM_WORLD;
int  qgps_current_task, qgps_master_task = 0, qgps_number_tasks;

qgps_block_t *qgps_blocks = NULL;
qgps_block_t * qgps_current_block = NULL;
qgps_block_t * qgps_current_transpose_block = NULL;
qgps_block_t *qgps_transpose_blocks = NULL;
ptrdiff_t qgps_local_size;

fftw_plan qgps_plan, qgps_inverse_plan;

int qgps_initialize_mpi(int argc, char **argv);
int qgps_initialize_blocks();
int qgps_cleanup_mpi();
int qgps_broadcast_block(qgps_block_t *block, int src_task);
int qgps_initialize_fftw();

int qgps_initialize(int argc, char **argv) {
        if (qgps_initialize_mpi(argc, argv))
                return 1;
        if (qgps_initialize_blocks())
                return 1;
        if (qgps_initialize_fftw())
                return 1;
        if (qgps_step_init())
                return 1;

        return 0;
}

int qgps_cleanup() {
        if (qgps_step_free())
                return 1;
        if (qgps_cleanup_mpi())
                return 1;

        return 0;
}

void qgps_exit() {
        qgps_cleanup();
        exit(EXIT_FAILURE);
}

int qgps_initialize_mpi(int argc, char **argv) {
        MPI_Init(&argc, &argv);
        MPI_Comm_dup(MPI_COMM_WORLD, &QGPS_COMM_WORLD);

        MPI_Comm_rank(QGPS_COMM_WORLD, &qgps_current_task);
        MPI_Comm_size(QGPS_COMM_WORLD, &qgps_number_tasks);

        return 0;
}

int qgps_initialize_blocks(int argc, char **argv) {
        qgps_blocks = calloc(qgps_number_tasks, sizeof(qgps_block_t));
        if (!qgps_blocks)
                return 1;

        qgps_transpose_blocks = calloc(qgps_number_tasks, sizeof(qgps_block_t));
        if (!qgps_transpose_blocks)
                return 1;

        return 0;
}

int qgps_cleanup_mpi() {
        MPI_Finalize();

        return 0;
}

int qgps_broadcast_block(qgps_block_t *block, int src_task) {
        MPI_Bcast(block, sizeof(qgps_block_t), MPI_BYTE, src_task,
                                                QGPS_COMM_WORLD);
        MPI_Barrier(QGPS_COMM_WORLD);

        return 0;
}

int qgps_initialize_fftw() {
        ptrdiff_t alloc_local, local_n0, local_0_start;
        qgps_block_t* block;

        fftw_mpi_init();

        qgps_local_size = fftw_mpi_local_size_2d(QGPS_NX, QGPS_NY/2 + 1,
                                                        QGPS_COMM_WORLD,
                                                        &local_n0,
                                                        &local_0_start);

        qgps_plan         = fftw_mpi_plan_dft_r2c_2d(QGPS_NX, QGPS_NY, 
                                                                NULL, NULL,
                                                                QGPS_COMM_WORLD,
                                                                FFTW_ESTIMATE);

        qgps_inverse_plan = fftw_mpi_plan_dft_c2r_2d(QGPS_NX, QGPS_NY,
                                                                NULL, NULL,
                                                                QGPS_COMM_WORLD,
                                                                FFTW_ESTIMATE);


        qgps_current_block = &qgps_blocks[qgps_current_task];
        qgps_current_transpose_block = &qgps_transpose_blocks[qgps_current_task];

        // FIXME: x_length and y_length are not well-understood

        // init block information for this task
        block = qgps_current_block;
        block->id       = qgps_current_task;
        block->x_begin  = local_0_start;
        block->x_end    = block->x_begin + local_n0;
        block->x_length = qgps_local_size / QGPS_NY;
        block->y_begin  = 0;
        block->y_end    = QGPS_NY;
        block->y_length = QGPS_NY;
        block->size     = qgps_local_size;

        block = qgps_current_transpose_block;
        block->id       = qgps_current_task;
        block->x_begin  = 0;
        block->x_end    = QGPS_NX;
        block->x_length = QGPS_NX;
        block->y_begin  = local_0_start;
        block->y_end    = block->y_begin + local_n0 / 2 + 1;
        block->y_length = qgps_local_size / QGPS_NX;
        block->size     = qgps_local_size;

        MPI_Barrier(QGPS_COMM_WORLD);

        // broadcast block information
        for (int task = 0; task < qgps_number_tasks; task++) {
                qgps_broadcast_block(&(qgps_blocks[task]), task);
        }

        return 0;
}

