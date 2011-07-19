#ifndef _QGPS_INPUT_H
#define _QGPS_INPUT_H

#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <iniparser.h>
#include "step.h"
#include "qgps-mpi.h"

/* initialization type */
extern qgps_init_type_t qgps_init_type;
/* initialization data */
extern complex *qgps_init_data;
/* output directory */
extern char *qgps_output_directory;
/* input configuration file */
extern char *qgps_configuration_file;
/* domain size */
extern int qgps_nx, qgps_ny;

extern double qgps_time_step;

int qgps_configure(int argc, char **argv);

#endif

