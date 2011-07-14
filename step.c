#include <stdlib.h>
#include "step.h"

const double qgps_time_start = 0.0;
const double qgps_time_end   = 1.0;
const double qgps_time_step  = 0.1;
double qgps_time = 0;

complex *psi_x;
complex *psi_y;
complex *omega;

int advection(complex *tracer_advt, complex *tracer);

int update_psi();

int qgps_step_init() {
        psi_x = fftw_alloc_complex(qgps_local_size);
        if (!psi_x)
                return 1;

        psi_y = fftw_alloc_complex(qgps_local_size);
        if (!psi_y)
                return 1;

        omega = fftw_alloc_complex(qgps_local_size);
        if (!omega)
                return 1;

        qgps_time = qgps_time_start;

        return 0;
}

int qgps_step_free() {
        fftw_free(psi_x);
        psi_x = NULL;

        fftw_free(psi_y);
        psi_y = NULL;

        fftw_free(omega);
        omega = NULL;

        return 0;
}


int qgps_step() {
        //Compute RK4 time step

        static complex *omega_t = NULL,
                       *work    = NULL;

        if (!omega_t)
                omega_t = fftw_alloc_complex(qgps_local_size);

        if (!work)
                work = fftw_alloc_complex(qgps_local_size);

        advection(work,omega);

        for (int idx = 0; idx < qgps_local_size; idx++) {
                omega_t[idx] = -work[idx] / 6.0;
                work[idx] = omega[idx] - work[idx] * qgps_time_step / 2.0;
        }
        advection(work,work);

        for (int idx = 0; idx < qgps_local_size; idx++) {
                omega_t[idx] -= work[idx] / 3.0;
                work[idx] = omega[idx] - work[idx] * qgps_time_step / 2.0;
        }
        advection(work,work);

        for (int idx = 0; idx < qgps_local_size; idx++) {
                omega_t[idx] -= work[idx] / 3.0;
                work[idx] = omega[idx] - work[idx] * qgps_time_step;
        }
        advection(work,work);

        for (int idx = 0; idx < qgps_local_size; idx++) {
                omega_t[idx] -= work[idx] / 6.0;
                omega[idx] += omega_t[idx] * qgps_time_step;
        }

        qgps_time += qgps_time_step;

        return 0;
}

int update_psi() {
        int     nx = qgps_current_transpose_block->x_length,
                ny = qgps_current_transpose_block->y_length,
                xb = qgps_current_transpose_block->x_begin,
                yb = qgps_current_transpose_block->y_begin;

        int k1, k2, k_sq, idx;

        for (int x = 0; x < nx; x++)
        for (int y = 0; y < ny; y++) {
                idx = y*nx + x;
                k1 = x + xb;
                k2 = y + yb;
                k_sq = (k1 * k1) + (k2 * k2);
                if (k_sq > 0) {
                        psi_x[idx] = I * (k1 / (double)k_sq) * omega[idx];
                        psi_y[idx] = I * (k2 / (double)k_sq) * omega[idx];
                }
                else {
                        // What do we do here?
                        psi_x[idx] = 0.0;
                        psi_y[idx] = 0.0;
                }
        }

        return 0;
}

int advection(complex *tracer_advt, complex *tracer) {
        int     nx = qgps_current_transpose_block->x_length,
                ny = qgps_current_transpose_block->y_length,
                xb = qgps_current_transpose_block->x_begin,
                yb = qgps_current_transpose_block->y_begin;

        static complex *tracer_kx = NULL,
                       *tracer_ky = NULL;

        if (!tracer_kx)
                tracer_kx = fftw_alloc_complex(qgps_local_size);
        if (!tracer_ky)
                tracer_ky = fftw_alloc_complex(qgps_local_size);

        static double *tracer_x = NULL,
                      *tracer_y = NULL;

        if (!tracer_x)
                tracer_x = fftw_alloc_real(qgps_local_size * 2);
        if (!tracer_y)
                tracer_y = fftw_alloc_real(qgps_local_size * 2);

        /* v_vel has reversed sign */
        static double *u_vel = NULL,
                      *v_vel = NULL;

        if (!u_vel)
                u_vel = fftw_alloc_real(qgps_local_size * 2);
        if (!v_vel)
                v_vel = fftw_alloc_real(qgps_local_size * 2);


        static double *advt_real = NULL;

        if (!advt_real)
                advt_real = fftw_alloc_real(qgps_local_size * 2);

        // Calculate the velocity in real space
        fftw_mpi_execute_dft_c2r(qgps_inverse_plan, psi_x, v_vel);
        fftw_mpi_execute_dft_c2r(qgps_inverse_plan, psi_y, u_vel);

        // Compute the gradient of tracer
        for (int x = 0; x < nx; x++)
        for (int y = 0; y < ny; y++) {
                int idx = y*nx + x;
                int k1 = x + xb;
                int k2 = y + yb;
                if (k1 > 0 || k2 > 0) {
                        tracer_kx[idx] = -I * k1 * tracer[idx];
                        tracer_ky[idx] = -I * k2 * tracer[idx];
                }
                else {
                        // What do we do here?
                        tracer_kx[idx] = 0.0;
                        tracer_ky[idx] = 0.0;
                }
        }

        // Compute inverse fft of tracer gradient
        fftw_mpi_execute_dft_c2r(qgps_inverse_plan, tracer_kx, tracer_x);
        fftw_mpi_execute_dft_c2r(qgps_inverse_plan, tracer_ky, tracer_y);

        // Calculate the advection in real space
        for (int idx = 0; idx < qgps_local_size; idx++) {
                advt_real[idx] = u_vel[idx] * tracer_x[idx] - v_vel[idx] * tracer_y[idx];
        }

        // Compute the fft of advection
        fftw_mpi_execute_dft_r2c(qgps_plan, advt_real, tracer_advt);

        return 0;
}

