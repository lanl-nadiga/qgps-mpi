#include <stdlib.h>
#include "qgps-mpi.h"

MPI_File qgps_output_file;

int qgps_output() {
        return qgps_output_open() || qgps_output_write() || qgps_output_close();
}
int qgps_transpose_r(double *data) {
        static fftw_plan plan = NULL;
        if (!plan)
                plan = fftw_mpi_plan_transpose(qgps_ny, qgps_nx,
                                data, data, QGPS_COMM_WORLD,
                                FFTW_ESTIMATE | FFTW_MPI_TRANSPOSED_IN);

        fftw_mpi_execute_r2r(plan, data, data);

        return 0;
}

char *qgps_output_filename() {
        static char *s = NULL;
        if (!s) {
                s = malloc(strlen(qgps_output_directory) + strlen("/qgps-mpi.0000.bin") + 1);
                strcpy(s, qgps_output_directory);
        }

        sprintf(s + strlen(qgps_output_directory), "/qgps-mpi.%.4i.bin", (int)(qgps_time / qgps_time_step));
        return s;
}

int qgps_output_open() {
        MPI_File_open(QGPS_COMM_WORLD, qgps_output_filename(),
                        MPI_MODE_CREATE | MPI_MODE_WRONLY,
                        MPI_INFO_NULL, &qgps_output_file);

        MPI_File_set_view(qgps_output_file, qgps_current_real_block->y_begin * qgps_nx,
                        MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL);
}

int qgps_output_close() {
        return MPI_File_close(&qgps_output_file);
}

int qgps_output_write() {
        static double *omega_real = NULL;

        if (!omega_real)
                omega_real = fftw_alloc_real(qgps_local_size * 2);

        qgps_dft_c2r(omega, omega_real);
        qgps_transpose_r(omega_real);

        qgps_block_t *b = qgps_current_real_block;

        int pad = 2 - qgps_ny % 2;
        for (int j = b->y_begin; j < b->y_end; j++) {
                int idx = (j - b->y_begin) * (qgps_nx + pad);

                MPI_File_write(qgps_output_file, &omega_real[idx], qgps_nx,
                                MPI_DOUBLE, MPI_STATUS_IGNORE);
        }
        return 0;
}
