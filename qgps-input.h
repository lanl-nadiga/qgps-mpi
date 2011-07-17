#ifndef _QGPS_INPUT_H
#define _QGPS_INPUT_H

#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <iniparser.h>
#include "step.h"
#include "qgps-mpi.h"

extern complex *qgps_init_data;
extern char *qgps_output_directory;
extern char *qgps_configuration_file;
extern int qgps_nx, qgps_ny;

int qgps_configure(int argc, char **argv);

#endif

