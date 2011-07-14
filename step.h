#ifndef _STEP_H
#define _STEP_H

#include "qgps-mpi.h"

// time step
extern double time;
const double time_start = 0.0;
const double time_end   = 1.0;
const double time_step  = 0.1;

// gradient of stream function in spectral space
extern complex  *psi_x,
                *psi_y;
// vorticity in spectral space
extern complex  *omega;

int qgps_step_init();
int qgps_step_free();
int qgps_step();

#endif

