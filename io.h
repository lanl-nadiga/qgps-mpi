#ifndef _QGPS_IO_H
#define _QGPS_IO_H

int qgps_open();
int qgps_write();
int qgps_close();

char *qgps_output_filename();

extern MPI_File qgps_output_file;

#endif
