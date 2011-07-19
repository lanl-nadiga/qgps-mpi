#include <stdlib.h>
#include <string.h>
#include "step.h"

const double qgps_time_start = 0.0;
const double qgps_time_end   = 1.0;
double qgps_time_step  = 0.001;
double qgps_time = 0;

complex *omega;

int gradient(complex *f, complex *dfdx, complex *dfdy);
int calc_vel(complex *vorticity, complex *uvel, complex *vvel);
int advection(complex *advt, complex *tracer, complex *uvel, complex *vvel);

int init_omega(qgps_init_type_t init_type);

int qgps_step_init() {
        omega = fftw_alloc_complex(qgps_local_size);
        if (!omega)
                return 1;

        qgps_time = qgps_time_start;

        init_omega(QGPS_INIT_DELTA_K);

        return 0;
}

int qgps_step_free() {
        fftw_free(omega);
        omega = NULL;

        return 0;
}


int qgps_step() {
        //Compute RK4 time step

        static complex *omega_t = NULL,
                       *work    = NULL,
                       *uvel    = NULL,
                       *vvel    = NULL;

        MPI_Barrier(QGPS_COMM_WORLD);
        if (!omega_t)
                omega_t = fftw_alloc_complex(qgps_local_size);
        if (!work)
                work = fftw_alloc_complex(qgps_local_size);
        if (!uvel)
                uvel = fftw_alloc_complex(qgps_local_size);
        if (!vvel)
                vvel = fftw_alloc_complex(qgps_local_size);

        calc_vel(omega, uvel, vvel);
        advection(work, omega, uvel, vvel);

        for (int idx = 0; idx < qgps_local_size; idx++) {
                omega_t[idx] = -work[idx] / 6.0;
                work[idx] = omega[idx] - work[idx] * qgps_time_step / 2.0;
        }
        calc_vel(work, uvel, vvel);
        advection(work, work, uvel, vvel);

        for (int idx = 0; idx < qgps_local_size; idx++) {
                omega_t[idx] -= work[idx] / 3.0;
                work[idx] = omega[idx] - work[idx] * qgps_time_step / 2.0;
        }
        calc_vel(work, uvel, vvel);
        advection(work, work, uvel, vvel);

        for (int idx = 0; idx < qgps_local_size; idx++) {
                omega_t[idx] -= work[idx] / 3.0;
                work[idx] = omega[idx] - work[idx] * qgps_time_step;
        }
        calc_vel(work, uvel, vvel);
        advection(work, work, uvel, vvel);

        for (int idx = 0; idx < qgps_local_size; idx++) {
                omega_t[idx] -= work[idx] / 6.0;
                omega[idx] += omega_t[idx] * qgps_time_step;
        }

        qgps_time += qgps_time_step;

        return 0;
}

int gradient(complex *f, complex *dfdx, complex *dfdy) {
        int idx;
        int k1, k2;

        int ib = qgps_current_complex_block->x_begin;
        int jb = qgps_current_complex_block->y_begin;
        int nx = qgps_current_complex_block->x_length;
        int ny = qgps_current_complex_block->y_length;

        for(int i = 0; i < nx; i++) {
        for(int j = 0; j < ny; j++) {
                idx = j*nx + i;

                // Calculate the local wavenumbers
                if(i <= nx/2) {
                  k1 = i+ib;
                }
                else {
                  k1 = i+ib - QGPS_NX;
                }
                k2 = j+jb;

                dfdx[idx] = I * k1*f[idx];
                dfdy[idx] = I * k2*f[idx];
        }}

        return 0;
}

double cabs_sqr(complex c) {
        return creal(c)*creal(c) + cimag(c)*cimag(c);
}

double l2_norm_squared(complex *f) {
        double tmp, norm = 0.0;
        int idx, k1, k2;

        int ib = qgps_current_complex_block->x_begin;
        int jb = qgps_current_complex_block->y_begin;
        int nx = qgps_current_complex_block->x_length;
        int ny = qgps_current_complex_block->y_length;

        for(int i = 0; i < nx; i++) {
        for(int j = 0; j < ny; j++) {
                idx = j*nx + i;
                k1 = i+ib;
                k2 = j+jb;

                if(k2 == 0) {
                  norm += cabs_sqr(f[idx]);
                }
                else if(QGPS_NX%2 == 1 && k1 == QGPS_NX/2 + 1) {
                  norm += cabs_sqr(f[idx]);
                }
                else if(QGPS_NY%2 == 1 && k2 == QGPS_NY/2 + 1) {
                  norm += cabs_sqr(f[idx]);
                }
                else {
                  norm += 2.0*cabs_sqr(f[idx]);
                }
        }}

        tmp = norm;
        MPI_Reduce(&tmp, &norm, 1, MPI_DOUBLE, MPI_SUM,
                                        qgps_master_task,
                                        QGPS_COMM_WORLD);
        MPI_Bcast(&norm, 1, MPI_DOUBLE, qgps_master_task,
                                        QGPS_COMM_WORLD);

        return norm;
}


int calc_vel(complex *vorticity, complex *uvel, complex *vvel) {
        int idx;
        int k1, k2, k_sq;

        int ib = qgps_current_complex_block->x_begin;
        int jb = qgps_current_complex_block->y_begin;
        int nx = qgps_current_complex_block->x_length;
        int ny = qgps_current_complex_block->y_length;

        // \omega = -\Delta\psi

        for(int i = 0; i < nx; i++) {
        for(int j = 0; j < ny; j++) {
                idx = j*nx + i;

                // Calculate the local wavenumbers
                if(i <= nx/2) {
                  k1 = i+ib;
                }
                else {
                  k1 = i+ib - QGPS_NX;
                }
                k2 = j+jb;
                k_sq = (k1 * k1) + (k2 * k2);

                if(k_sq > 0) {
                        vvel[idx] = -I * k1*omega[idx]/(double)k_sq;
                        uvel[idx] =  I * k2*omega[idx]/(double)k_sq;
                }
                else {
                        uvel[idx] = 0.0;
                        vvel[idx] = 0.0;
                }
        }}

        return 0;
}

int advection(complex *advt, complex *tracer, complex *uvel, complex *vvel) {
        /* Takes in the variables
         *  tracer      (a tracer quantity, such as vorticity)
         *  uvel        (x directional flow velocity)
         *  vvel        (y directional flow velocity)
         *  and produces the variable
         *  advt        (tracer advection)
         */

         static double   *uvel_real = NULL,      // real x directional velocity
                        *vvel_real = NULL,      // real y directional velocity
                        *dtdx_real = NULL,      // real x gradient of tracer
                        *dtdy_real = NULL,      // real y gradient of tracer
                        *advt_real = NULL;      // real tracer advection

         static complex  *dtdx = NULL,   // complex x gradient of tracer
                        *dtdy = NULL;   // complex y gradient of tracer

        // initialize the arrays used for advection calculation
        if(!dtdx)
                dtdx = fftw_alloc_complex(qgps_local_size);
        if(!dtdy)
                dtdy = fftw_alloc_complex(qgps_local_size);
        if(!dtdx_real)
                dtdx_real = fftw_alloc_real(qgps_local_size*2);
        if(!dtdy_real)
                dtdy_real = fftw_alloc_real(qgps_local_size*2);
        if(!uvel_real)
                uvel_real = fftw_alloc_real(qgps_local_size*2);
        if(!vvel_real)
                vvel_real = fftw_alloc_real(qgps_local_size*2);
        if(!advt_real)
                advt_real = fftw_alloc_real(qgps_local_size*2);

        // calculate the velocity in real space
        qgps_dft_c2r(uvel, uvel_real);
        qgps_dft_c2r(vvel, vvel_real);

        // compute the gradient of tracer
        gradient(tracer, dtdx, dtdy);

        // compute the tracer gradient in real space
        qgps_dft_c2r(dtdx, dtdx_real);
        qgps_dft_c2r(dtdy, dtdy_real);

        // calculate the advection in real space
        for (int idx = 0; idx < 2*qgps_local_size; idx++) {
                advt_real[idx]  = uvel_real[idx] * dtdx_real[idx]
                                + vvel_real[idx] * dtdy_real[idx];
        }

        // convert the advection to Fourier space
        qgps_dft_r2c(advt_real, advt);

        return 0;
}

void qgps_init_delta_k() {
        for(int idx = 0; idx < qgps_local_size; idx++) {
                omega[idx] = 0.0;
        }

        if(qgps_current_real_block->x_begin == 0) {
                omega[QGPS_NX+1] = 0.5 + 0.5*I;
                omega[2*QGPS_NX+1] = -0.5 + 0.5*I;
        }
}


int init_omega(qgps_init_type_t init_type) {
        switch (init_type) {
                case QGPS_INIT_RESTART:
                        fprintf(stderr,"init option not supported yet.\n");
                        qgps_exit(EXIT_FAILURE);
                        break;
                case QGPS_INIT_DELTA_K:
                        qgps_init_delta_k();
                        break;
                default:
                        fprintf(stderr,"unknown init option.\n");
                        qgps_exit(EXIT_FAILURE);
                        break;
        }
}

qgps_init_type_t qgps_init_type_parse(const char *string) {
        if(!strcmp("delta", string))
                return QGPS_INIT_DELTA_K;
        else {
                fprintf(stderr, "Unknown init type %s\n", string);
                qgps_exit(EXIT_FAILURE);
                return -1;
        }
}
