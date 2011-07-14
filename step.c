#include <stdlib.h>
#include "step.h"

double dt = 0.1;

complex *psi_x;
complex *psi_y;
complex *omega;

int advection(complex *tracer_advt, complex *tracer);

int update_psi();

int qgps_step_init() {
        int local_size = qgps_blocks[qgps_current_task].size;

        psi_x = (complex *) fftw_alloc_complex(local_size);
        psi_y = (complex *) fftw_alloc_complex(local_size);
        omega = (complex *) fftw_alloc_complex(local_size);
}

int qgps_step_free() {
        free(psi_x);
        free(psi_y);
        free(omega);
}


int qgps_step() {
        return 0;
}

int update_psi() {
        int     nx = qgps_transpose_blocks[qgps_current_task].x_length,
                ny = qgps_transpose_blocks[qgps_current_task].y_length,
                xb = qgps_transpose_blocks[qgps_current_task].x_begin,
                yb = qgps_transpose_blocks[qgps_current_task].y_begin;

        int k1, k2, k_sq, idx;

        for(int x = 0; x < nx; x++)
        for(int y = 0; y < ny; y++) {
                idx = y*nx + x;
                k1 = x + xb;
                k2 = y + yb;
                k_sq = (k1*k1 + k2*k2);
                if(k_sq > 0) {
                        psi_x[idx] = I*(k1 / (double)k_sq) * omega[idx];
                        psi_y[idx] = I*(k2 / (double)k_sq) * omega[idx];
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
        int local_size = qgps_blocks[qgps_current_task].size;

        int     nx = qgps_transpose_blocks[qgps_current_task].x_length,
                ny = qgps_transpose_blocks[qgps_current_task].y_length,
                xb = qgps_transpose_blocks[qgps_current_task].x_begin,
                yb = qgps_transpose_blocks[qgps_current_task].y_begin;

        complex *tracer_kx, *tracer_ky;

        double *tracer_x, *tracer_y;

        double *u_vel, *v_vel; // v_vel has sign reversed

        double *advt_real;

        tracer_kx = fftw_alloc_complex(local_size);
        tracer_ky = fftw_alloc_complex(local_size);
        tracer_x = fftw_alloc_real(local_size*2);
        tracer_y = fftw_alloc_real(local_size*2);
        u_vel = fftw_alloc_real(local_size*2);
        v_vel = fftw_alloc_real(local_size*2);
        advt_real = fftw_alloc_real(local_size*2);

        // Calculate the velocity in real space
        fftw_mpi_execute_dft_c2r(inverse_plan, psi_x, v_vel);
        fftw_mpi_execute_dft_c2r(inverse_plan, psi_y, u_vel);

        // Compute the gradient of tracer
        for(int x = 0; x < nx; x++)
        for(int y = 0; y < ny; y++) {
                int idx = y*nx + x;
                int k1 = x + xb;
                int k2 = y + yb;
                if(k1 > 0 || k2 > 0) {
                        tracer_kx[idx] = -I*k1 * tracer[idx];
                        tracer_ky[idx] = -I*k2 * tracer[idx];
                }
                else {
                        // What do we do here?
                        tracer_kx[idx] = 0.0;
                        tracer_ky[idx] = 0.0;
                }
        }

        // Compute inverse fft of tracer gradient
        fftw_mpi_execute_dft_c2r(inverse_plan, tracer_kx, tracer_x);
        fftw_mpi_execute_dft_c2r(inverse_plan, tracer_ky, tracer_y);

        // Calculate the advection in real space
        for(int idx = 0; idx < local_size; idx++) {
                advt_real[idx] = u_vel[idx]*tracer_x[idx] - v_vel[idx]*tracer_y[idx];
        }

        // Compute the fft of advection
        fftw_mpi_execute_dft_r2c(plan, advt_real, tracer_advt);

        free(tracer_kx);
        free(tracer_ky);
        free(tracer_x);
        free(tracer_y);
        free(u_vel);
        free(v_vel);

        return 0;
}

