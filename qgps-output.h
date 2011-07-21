#ifndef _QGPS_OUTPUT_H
#define _QGPS_OUTPUT_H
#include <mpi.h>

int qgps_output();

MPI_File qgps_output_open_file(char *filename);
int qgps_output_close_file(MPI_File *file);
int qgps_output_write_complex(MPI_File file, complex *data);
int qgps_output_write_real(MPI_File file, double *data);
int qgps_output_complex(char *filename, complex *data);
int qgps_output_real(char *filename, double *data);

char *qgps_output_filename();

#endif
