#ifndef _STEP_H
#define _STEP_H

#include "qgps-mpi.h"

// time step
extern double qgps_time;
extern const double qgps_time_start;
extern const double qgps_time_end;
extern const double qgps_time_step;

// gradient of stream function in spectral space
extern complex  *psi_x,
                *psi_y;
// vorticity in spectral space
extern complex  *omega;

int qgps_step_init();
int qgps_step_free();
int qgps_step();

#endif

