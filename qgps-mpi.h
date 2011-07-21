/* 
 * Header file for exported local variables
 */

#ifndef _QGPS_MPI_H
#define _QGPS_MPI_H
#include <mpi.h>
#include <complex.h>
#include <fftw3-mpi.h>
#include "step.h"
#include "qgps-output.h"
#include "qgps-input.h"

/* local domain block */
typedef struct {
        int id;

        int x_begin, x_end, x_length;
        int y_begin, y_end, y_length;

        int size;
} qgps_block_t;

extern qgps_block_t * qgps_real_blocks;
extern qgps_block_t * qgps_complex_blocks;
extern qgps_block_t * qgps_current_real_block;
extern qgps_block_t * qgps_current_complex_block;

int qgps_local_nx();
int qgps_local_ny();
int qgps_x(int i);
int qgps_y(int i);
int qgps_x_t(int i);
int qgps_y_t(int i);
int qgps_index(int x, int y);
int qgps_index_t(int y, int x);
complex qgps_z(int i);

extern ptrdiff_t qgps_local_size;

/* MPI task variables */
extern MPI_Comm QGPS_COMM_WORLD;
extern int qgps_current_task, qgps_master_task, qgps_number_tasks;
extern int qgps_is_master_task;


/* FFTW execute functions */
int qgps_dft_c2r(const complex *in, double *out);
int qgps_dft_r2c(const double *in, complex *out);
int qgps_transpose_r(double *data);

int qgps_initialize(int argc, char **argv);
int qgps_cleanup();

void qgps_exit(int status);

#endif
