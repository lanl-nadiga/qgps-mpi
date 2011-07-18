#ifndef _STEP_H
#define _STEP_H

#include "qgps-mpi.h"

// omega initialization types
typedef enum {
        QGPS_INIT_RESTART,
        QGPS_INIT_DELTA_K,
} qgps_init_type_t;

qgps_init_type_t qgps_init_type_parse(const char *name);
extern complex *qgps_init_data;

extern qgps_init_type_t qgps_init_type;

// time step
extern double qgps_time;
extern const double qgps_time_start;
extern const double qgps_time_end;
extern double qgps_time_step;

// vorticity in Fourier space
extern complex  *omega;

int qgps_step_init();
int qgps_step_free();
int qgps_step();

double l2_norm_squared(complex *f);

#endif

