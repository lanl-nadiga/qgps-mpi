#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "qgps-mpi.h"
#include "qgps-input.h"
#include "diagnostics.h"
#include "step.h"

#define M_PI 3.1415926535897932384626433832795028841971693993751058209749445923

const double qgps_time_start = 0.0;
const double qgps_time_end   = 100.0;
const double qgps_error_max  = 5e-3;
const double qgps_time_step_min = 1e-15;
const double qgps_time_step_max = 1e-1;
double qgps_time = 0;

// FIXME: our scaling parameter for viscous forcing
double k_dispersion = 0;
double dispersion_exponent = 8;
int dispersion_cutoff = 1;

complex *omega;

int *qgps_kx, *qgps_ky, *qgps_k_sq;    // Array of wave numbers
complex *qgps_k = NULL;
double *cwgt;           // weighting for global sums

int calc_vel(complex *vorticity, complex *uvel, complex *vvel);
int advection(complex *advt, complex *tracer, complex *uvel, complex *vvel);
int qgps_rk54(complex *omega_t, double *time_step);
int viscous_forcing(complex *tracer);
int cutoff_high_frequencies(complex *tracer);

int init_k_dispersion();
int init_omega(qgps_init_type_t init_type);

int qgps_step_init() {
        omega = fftw_alloc_complex(qgps_local_size);
        if (!omega)
                return 1;

        qgps_time = qgps_time_start;

        // initialize the wave number arrays
        // !! must be done before call to init_omega
        qgps_k       = calloc(qgps_local_size, sizeof(complex));
        qgps_kx      = calloc(qgps_local_size,sizeof(int));
        qgps_ky      = calloc(qgps_local_size,sizeof(int));
        qgps_k_sq    = calloc(qgps_local_size,sizeof(int));
        cwgt    = calloc(qgps_local_size,sizeof(double));

        int ib = qgps_current_complex_block->x_begin;
        int jb = qgps_current_complex_block->y_begin;
        int nx = qgps_current_complex_block->x_length;
        int ny = qgps_current_complex_block->y_length;

        for(int i = 0; i < nx; i++)
        for(int j = 0; j < ny; j++) {
                int idx = j*nx + i;

                if(i <= nx/2) {
                  qgps_kx[idx] = i+ib;
                }
                else {
                  qgps_kx[idx] = i+ib - qgps_nx;
                }

                qgps_ky[idx] = j+jb;

                qgps_k[idx] = qgps_kx[idx] + I * qgps_ky[idx];
                qgps_k_sq[idx] = qgps_kx[idx]*qgps_kx[idx] + qgps_ky[idx]*qgps_ky[idx];

                if(qgps_ky[idx] == 0) {
                        cwgt[idx] = 1.0;
                }
                else if(qgps_nx%2 == 1&& qgps_kx[idx] == qgps_nx/2) {
                        cwgt[idx] = 1.0;
                }
                else if(qgps_ny%2 == 1&& qgps_ky[idx] == qgps_ny/2) {
                        cwgt[idx] = 1.0;
                }
                else {
                        cwgt[idx] = 2.0;
                }
        }

        init_k_dispersion();

        // Initialize the vorticity
        init_omega(QGPS_INIT_PATCHES);

        return 0;
}

int init_k_dispersion() {

        if(dispersion_cutoff) {
                dispersion_cutoff = ((qgps_ny/2)/3)*2 + 1;
                k_dispersion = dispersion_cutoff;
                k_dispersion *= 0.90;
        }
        else {
                k_dispersion = qgps_ny/2;
                k_dispersion *= 0.65;
        }

        return 0;
}

int qgps_step_free() {
        fftw_free(omega);
        omega = NULL;

        free(qgps_kx);
        free(qgps_ky);
        free(qgps_k_sq);
        free(cwgt);

        return 0;
}


int qgps_step() {
        /*
         * update the vorticity and time to qgps_time + qgps_time_step using
         * by breaking up qgps_time_step into smaller time intervals based on
         * the value of qgps_error_max
         */

        static double *cached_dt_guess = NULL;
        if (!cached_dt_guess) {
                cached_dt_guess = malloc(sizeof(double));
                *cached_dt_guess = 0.1*qgps_time_step;
        }

        double dt_guess = 0.1*qgps_time_step;
        double dt_total = 0.0,
               dt_max   = 0.0;

        static complex *omega_t = NULL;

        if(!omega_t)
                omega_t = fftw_alloc_complex(qgps_local_size);

        while(dt_total < qgps_time_step - 1e-15) {
                // calculate the forcing based on dt_guess
                qgps_rk54(omega_t, &dt_guess);

                if(dt_guess < qgps_time_step_min) {
                        dt_guess = qgps_time_step_min;
                }
                else if(dt_guess > qgps_time_step_max) {
                        dt_guess = qgps_time_step_max;
                }
                if(dt_guess + dt_total > qgps_time_step)  {
                        dt_guess = qgps_time_step - dt_total;
                }

                // update omega based on forcing and dt_guess
                for(int idx = 0; idx < qgps_local_size; idx++) {
                        omega[idx] += dt_guess*omega_t[idx];
                }
                dt_total += dt_guess;
                viscous_forcing(omega);
                qgps_diagnostics_update();
        }
        qgps_time += qgps_time_step;
        *cached_dt_guess = dt_guess;

        cutoff_high_frequencies(omega);
}

int qgps_rk54(complex *omega_t, double *dt) {
        /*
         * uses the Runge–Kutta–Fehlberg (rk54) algorithm to compute the
         * forcing and optimal time step on the fly based on a user-defined
         * error tolerance level.
         */
        double time_step = *dt;

        const double    rk4b1 = 25.0/216.0,
                        rk4b3 = 1408.0/2565.0,
                        rk4b4 = 2197.0/4104.0,
                        rk4b5 = -1.0/5.0;

        const double    rk5b1 = 16.0/135.0,
                        rk5b3 = 6656.0/12825.0,
                        rk5b4 = 28561.0/56430.0,
                        rk5b5 = -9.0/50.0,
                        rk5b6 = 2.0/55.0;

        const double    rk5a21= 1.0/4.0,
                        rk5a31= 3.0/32.0,
                        rk5a32= 9.0/32.0,
                        rk5a41= 1932.0/2197.0,
                        rk5a42= -7200.0/2197.0,
                        rk5a43= 7296.0/2197.0,
                        rk5a51= 439.0/216.0,
                        rk5a52= -8.0,
                        rk5a53= 3680.0/513.0,
                        rk5a54= -845.0/4104.0,
                        rk5a61= -8.0/27.0,
                        rk5a62= 2.0,
                        rk5a63= -3544.0/2565.0,
                        rk5a64= 1859.0/4104.0,
                        rk5a65= -11.0/40.0;

        static complex *work1   = NULL,
                       *work2   = NULL,
                       *work3   = NULL,
                       *work4   = NULL,
                       *work5   = NULL,
                       *work6   = NULL,
                       *uvel    = NULL,
                       *vvel    = NULL;

        MPI_Barrier(QGPS_COMM_WORLD);
        if (!work1)
                work1 = fftw_alloc_complex(qgps_local_size);
        if (!work2)
                work2 = fftw_alloc_complex(qgps_local_size);
        if (!work3)
                work3 = fftw_alloc_complex(qgps_local_size);
        if (!work4)
                work4 = fftw_alloc_complex(qgps_local_size);
        if (!work5)
                work5 = fftw_alloc_complex(qgps_local_size);
        if (!work6)
                work6 = fftw_alloc_complex(qgps_local_size);
        if (!uvel)
                uvel = fftw_alloc_complex(qgps_local_size);
        if (!vvel)
                vvel = fftw_alloc_complex(qgps_local_size);

        for (int idx = 0; idx < qgps_local_size; idx++) {
                work1[idx]      = omega[idx];
        }
        calc_vel(work1, uvel, vvel);
        advection(work1, work1, uvel, vvel);

        for (int idx = 0; idx < qgps_local_size; idx++) {
                work2[idx]      = work1[idx] * rk5a21;
                work2[idx]      = omega[idx] - time_step*work2[idx];
        }
        calc_vel(work2, uvel, vvel);
        advection(work2, work2, uvel, vvel);

        for (int idx = 0; idx < qgps_local_size; idx++) {
                work3[idx]      = work1[idx] * rk5a31
                                + work2[idx] * rk5a32;
                work3[idx]      = omega[idx] - time_step*work3[idx];
        }
        calc_vel(work3, uvel, vvel);
        advection(work3, work3, uvel, vvel);

        for (int idx = 0; idx < qgps_local_size; idx++) {
                work4[idx]      = work1[idx] * rk5a41
                                + work2[idx] * rk5a42
                                + work3[idx] * rk5a43;
                work4[idx]      = omega[idx] - time_step*work4[idx];
        }
        calc_vel(work4, uvel, vvel);
        advection(work4, work4, uvel, vvel);

        for (int idx = 0; idx < qgps_local_size; idx++) {
                work5[idx]      = work1[idx] * rk5a51
                                + work2[idx] * rk5a52
                                + work3[idx] * rk5a53
                                + work4[idx] * rk5a54;
                work5[idx]      = omega[idx] - time_step*work5[idx];
        }
        calc_vel(work5, uvel, vvel);
        advection(work5, work5, uvel, vvel);

        for (int idx = 0; idx < qgps_local_size; idx++) {
                work6[idx]      = work1[idx] * rk5a61
                                + work2[idx] * rk5a62
                                + work3[idx] * rk5a63
                                + work4[idx] * rk5a64
                                + work5[idx] * rk5a65;
                work6[idx]      = omega[idx] - time_step*work6[idx];
        }
        calc_vel(work6, uvel, vvel);
        advection(work6, work6, uvel, vvel);

        for (int idx = 0; idx < qgps_local_size; idx++) {
                // Calculate the RK5 estimate for the forcing
                work6[idx]      = work1[idx] * rk5b1
                                + work3[idx] * rk5b3
                                + work4[idx] * rk5b4
                                + work5[idx] * rk5b5
                                + work6[idx] * rk5b6;
                work6[idx] *= -1.0;

                // Calculate the RK4 estimate of the forcing
                omega_t[idx]    = work1[idx] * rk4b1
                                + work3[idx] * rk4b3
                                + work4[idx] * rk4b4
                                + work5[idx] * rk4b5;
                omega_t[idx] *= -1.0;

                // Calculate the difference between RK4 and RK5
                work6[idx] -= omega_t[idx];
        }

        // Calculate the optimal time step
        double err_max_sq = complex_global_max_squared(work6);

        double s = pow(0.25*qgps_error_max*qgps_error_max/err_max_sq,0.125);
        if ( s < 0.1 ) {
                s = 0.1;
        }
        else if ( s > 2.0 ) {
                s = 2.0;
        }

        time_step *= s;

        fprintf(stderr,"OPTIMAL T: %1.16lf %1.16lf %1.16lf %1.16f\n", *dt, s, time_step, err_max_sq);

        *dt = time_step;
        return 0;
}

int viscous_forcing(complex *tracer) {
        const qgps_block_t const *b = qgps_current_complex_block;
        const double k_dsp_sq = k_dispersion * k_dispersion;

        for (int i = 0; i < b->x_length; i++)
        for (int j = 0; j < b->y_length; j++) {
                int idx = j * b->x_length + i;

                tracer[idx] /= 1.0 + pow(qgps_k_sq[idx] / k_dsp_sq,
                                        dispersion_exponent);

        }
}

int cutoff_high_frequencies(complex *tracer) {
        if(! dispersion_cutoff) return 0;

        for(int idx = 0; idx < qgps_local_size; idx++) {
                if(qgps_kx[idx] > dispersion_cutoff || qgps_ky[idx] > dispersion_cutoff)
                        tracer[idx] = 0.0;
        }

        return 0;
}


int gradient(complex *f, complex *dfdx, complex *dfdy) {
        int idx;

        int nx = qgps_current_complex_block->x_length;
        int ny = qgps_current_complex_block->y_length;

        for(int i = 0; i < nx; i++)
        for(int j = 0; j < ny; j++) {
                idx = j*nx + i;

                dfdx[idx] = I * qgps_kx[idx]*f[idx];
                dfdy[idx] = I * qgps_ky[idx]*f[idx];
        }

        return 0;
}

int laplacian(complex *f, complex *delf) {
        int idx;

        int nx = qgps_current_complex_block->x_length;
        int ny = qgps_current_complex_block->y_length;

        for(int i = 0; i < nx; i++)
        for(int j = 0; j < ny; j++) {
                idx = j*nx + i;

                delf[idx] = - qgps_k_sq[idx]*f[idx];
        }

        return 0;
}

double cabs_sqr(complex c) {
        return creal(c)*creal(c) + cimag(c)*cimag(c);
}

double l2_norm_squared(complex *f) {
        double tmp, norm = 0.0;
        int idx;

        int nx = qgps_current_complex_block->x_length;
        int ny = qgps_current_complex_block->y_length;

        // integrate on the local task
        for(int i = 0; i < nx; i++) {
        for(int j = 0; j < ny; j++) {
                idx = j*nx + i;

                norm += cabs_sqr(f[idx])*cwgt[idx];
        }}

        // sum across tasks
        tmp = norm;
        MPI_Reduce(&tmp, &norm, 1, MPI_DOUBLE, MPI_SUM,
                                        qgps_master_task,
                                        QGPS_COMM_WORLD);
        MPI_Bcast(&norm, 1, MPI_DOUBLE, qgps_master_task,
                                        QGPS_COMM_WORLD);
        // integration pwnd
        return norm;
}

complex complex_integral(complex *f) {
        double tmp, total = 0.0;
        int idx;

        int nx = qgps_current_complex_block->x_length;
        int ny = qgps_current_complex_block->y_length;

        // integrate on the local task
        for(int i = 0; i < nx; i++) {
        for(int j = 0; j < ny; j++) {
                idx = j*nx + i;

                total += f[idx]*cwgt[idx];
        }}

        // sum across tasks
        tmp = total;
        MPI_Reduce(&tmp, &total, 1, MPI_COMPLEX, MPI_SUM,
                                        qgps_master_task,
                                        QGPS_COMM_WORLD);
        MPI_Bcast(&total, 1, MPI_DOUBLE, qgps_master_task,
                                        QGPS_COMM_WORLD);
        // integration pwnd
        return total;
}

double complex_global_max_squared(complex *f) {
        double tmp, max = 0.0;
        int idx;

        int nx = qgps_current_complex_block->x_length;
        int ny = qgps_current_complex_block->y_length;

        // integrate on the local task
        for(int i = 0; i < nx; i++) {
        for(int j = 0; j < ny; j++) {
                idx = j*nx + i;

                tmp = cabs_sqr(f[idx]);

                if(tmp > max) max = tmp;
        }}

        // find max of all tasks
        tmp = max;
        MPI_Reduce(&tmp, &max, 1, MPI_DOUBLE, MPI_MAX,
                                        qgps_master_task,
                                        QGPS_COMM_WORLD);
        MPI_Bcast(&max, 1, MPI_DOUBLE, qgps_master_task,
                                        QGPS_COMM_WORLD);
        // global max pwnd
        return max;
}

int calc_vel(complex *vorticity, complex *uvel, complex *vvel) {
        int idx;

        int nx = qgps_current_complex_block->x_length;
        int ny = qgps_current_complex_block->y_length;

        // \omega = -\Delta\psi

        for(int i = 0; i < nx; i++)
        for(int j = 0; j < ny; j++) {
                idx = j*nx + i;

                if(qgps_k_sq[idx] > 0) {
                        uvel[idx] = -I * qgps_ky[idx]*omega[idx]/(double)qgps_k_sq[idx];
                        vvel[idx] =  I * qgps_kx[idx]*omega[idx]/(double)qgps_k_sq[idx];
                }
                else {
                        uvel[idx] = 0.0;
                        vvel[idx] = 0.0;
                }
        }

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

         static double  *uvel_real = NULL,      // real x directional velocity
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

        // compute the gradient of tracer
        gradient(tracer, dtdx, dtdy);

        // smooth fields
        cutoff_high_frequencies(uvel);
        cutoff_high_frequencies(vvel);
        cutoff_high_frequencies(dtdx);
        cutoff_high_frequencies(dtdy);

        // calculate the velocity in real space
        qgps_dft_c2r(uvel, uvel_real);
        qgps_dft_c2r(vvel, vvel_real);

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
        /*
         * Initialize a delta funciton at a specific wave number
         */
        int idx, pad = 2 - qgps_ny%2;

        double a, b;

        int qgps_kx = 4, qgps_ky = 4;;

        complex k_amp = 0.5 + 0.5*I;

        int ib = qgps_current_real_block->x_begin;
        int jb = qgps_current_real_block->y_begin;
        int ie = qgps_current_real_block->x_end;
        int je = qgps_current_real_block->y_end;

        double *omega_real = fftw_alloc_real(qgps_local_size*2);

        for(int i = ib; i < ie; i++)
        for(int j = jb; j < je; j++) {
                a = 2*M_PI*i*qgps_kx/(double)qgps_nx;
                b = 2*M_PI*j*qgps_ky/(double)qgps_ny;
                idx = (i-ib)*(qgps_ny + pad) + (j-jb);
                omega_real[idx] = 2.0*creal(k_amp)*cos(a)*cos(b)
                                - 2.0*creal(k_amp)*sin(a)*sin(b)
                                - 2.0*cimag(k_amp)*sin(a)*cos(b)
                                - 2.0*cimag(k_amp)*cos(a)*sin(b);
        }

        qgps_dft_r2c(omega_real,omega);

        fftw_free(omega_real);
}

void qgps_init_patches() {
        /*
         * Initialize a delta funciton at a specific wave number
         */
        int idx, pad = 2 - qgps_ny%2;

        double a, b;

        int ib = qgps_current_real_block->x_begin;
        int jb = qgps_current_real_block->y_begin;
        int ie = qgps_current_real_block->x_end;
        int je = qgps_current_real_block->y_end;

        double *omega_real = fftw_alloc_real(qgps_local_size*2);

        for(int i = ib; i < ie; i++)
        for(int j = jb; j < je; j++) {
                a = 2*M_PI*i*3/(double)qgps_nx;
                b = 2*M_PI*j*3/(double)qgps_ny;
                idx = (i-ib)*(qgps_ny + pad) + (j-jb);
                omega_real[idx] = (cos(a) + sin(a))
                                * (cos(b) + sin(b));
        }

        qgps_dft_r2c(omega_real,omega);

        complex *specific_vorticity = fftw_alloc_complex(qgps_local_size);
        specific_vorticity[0] = 0;
        for (int i = 1; i < qgps_local_size; i++)
                specific_vorticity[i] = omega[i] / qgps_k[i];

        double total_energy = l2_norm_squared(specific_vorticity) / 2;

        omega[0] /= sqrt(total_energy);
        for (int i = 1; i < qgps_local_size; i++)
                omega[i] /= sqrt(total_energy);

        fftw_free(omega_real);
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
                case QGPS_INIT_PATCHES:
                        qgps_init_patches();
                        break;
                default:
                        fprintf(stderr,"unknown init option.\n");
                        qgps_exit(EXIT_FAILURE);
                        break;
        }
}

