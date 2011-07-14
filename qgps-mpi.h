/* 
 * Header file for exported local variables
 */

#ifndef _QGPS_MPI_H
#define _QGPS_MPI_H
#include <mpi.h>
#include <complex.h>
#include <fftw3-mpi.h>

/* domain size */
extern const int QGPS_NX, QGPS_NY;

/* local domain block */
typedef struct {
        int id;

        int x_begin, x_end, x_length;
        int y_begin, y_end, y_length;

        int size;
} qgps_block_t;

extern qgps_block_t *qgps_blocks;
extern qgps_block_t *qgps_transpose_blocks;

/* MPI task variables */
extern MPI_Comm QGPS_COMM_WORLD;
extern int qgps_current_task, qgps_master_task, qgps_number_tasks;

extern qgps_block_t * qgps_current_block;
extern qgps_block_t * qgps_current_transpose_block;
extern ptrdiff_t qgps_local_size;

/* FFTW task variables */
extern fftw_plan qgps_plan, qgps_inverse_plan;

int qgps_initialize(int argc, char **argv);
int qgps_cleanup();

#endif
