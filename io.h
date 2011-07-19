#ifndef _QGPS_IO_H
#define _QGPS_IO_H

int qgps_output_open();
int qgps_output_write();
int qgps_output_close();
int qgps_output();

char *qgps_output_filename();

extern MPI_File qgps_output_file;

#endif
