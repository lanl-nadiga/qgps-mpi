#include <stdlib.h>
#include "step.h"

const double time_start = 0.0;
const double time_end   = 1.0;
const double time_step  = 0.1;
double time = 0;

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

        time = time_start;

        return 0;
}

int qgps_step_free() {
        free(psi_x);
        psi_x = NULL;

        free(psi_y);
        psi_y = NULL;

        free(omega);
        omega = NULL;

        return 0;
}


int qgps_step() {
        //Compute RK4 time step

        int idx;

        complex *omega_t;

        complex *work;

        omega_t = fftw_alloc_complex(qgps_local_size);
        work = fftw_alloc_complex(qgps_local_size);

        advection(work,omega);

        for(idx = 0; idx < qgps_local_size; idx++) {
                omega_t[idx] = -work[idx]/6.0;
                work[idx] = omega[idx] - work[idx]*time_step/2.0;
        }
        advection(work,work);

        for(idx = 0; idx < qgps_local_size; idx++) {
                omega_t[idx] -= work[idx]/3.0;
                work[idx] = omega[idx] - work[idx]*time_step/2.0;
        }
        advection(work,work);

        for(idx = 0; idx < qgps_local_size; idx++) {
                omega_t[idx] -= work[idx]/3.0;
                work[idx] = omega[idx] - work[idx]*time_step;
        }
        advection(work,work);

        for(idx = 0; idx < qgps_local_size; idx++) {
                omega_t[idx] -= work[idx]/6.0;
                omega[idx] += omega_t[idx]*time_step;
        }

        time += time_step;

        free(omega_t);
        free(work);

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

        complex *tracer_kx = fftw_alloc_complex(qgps_local_size),
                *tracer_ky = fftw_alloc_complex(qgps_local_size);

        double *tracer_x = fftw_alloc_real(qgps_local_size * 2),
               *tracer_y = fftw_alloc_real(qgps_local_size * 2);

        /* v_vel has reversed sign */
        double *u_vel = fftw_alloc_real(qgps_local_size * 2),
               *v_vel = fftw_alloc_real(qgps_local_size * 2);

        double *advt_real = fftw_alloc_real(qgps_local_size * 2);

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

        free(tracer_kx);
        free(tracer_ky);
        free(tracer_x);
        free(tracer_y);
        free(u_vel);
        free(v_vel);

        return 0;
}

