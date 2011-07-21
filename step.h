#include <complex.h>
#ifndef _STEP_H
#define _STEP_H

// time step
extern double qgps_time;
extern const double qgps_time_start;
extern const double qgps_time_end;

// vorticity in Fourier space
extern complex *omega;

extern int *qgps_kx, *qgps_ky, *qgps_k_sq;
extern complex *qgps_k;

int qgps_step_init();
int qgps_step_free();
int qgps_step();

int gradient(complex *f, complex *dfdx, complex *dfdy);

double l2_norm_squared(complex *f);
double complex_global_max_squared(complex *f);

#endif

