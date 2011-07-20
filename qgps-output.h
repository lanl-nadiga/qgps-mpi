#ifndef _QGPS_OUTPUT_H
#define _QGPS_OUTPUT_H

int qgps_output_open();
int qgps_output_write();
int qgps_output_close();
int qgps_output();

char *qgps_output_filename();

#endif
