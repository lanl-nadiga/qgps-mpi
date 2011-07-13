/* 
 * Header file for exported local variables
 */

#ifndef _QGPS_MPI_H
#define _QGPS_MPI_H

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

/* MPI task variables */
extern const int QGPS_COMM_WORLD;
extern int qgps_current_task, qgps_master_task, qgps_number_tasks;

int qgps_initialize(int argc, char **argv);
int qgps_cleanup();

#endif
