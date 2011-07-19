#include <stdlib.h>
#include "qgps-mpi.h"

MPI_File qgps_output_file;

char *qgps_output_filename() {
        static char *s = NULL;
        if (!s)
                s = malloc(256);

        sprintf(s, "./qgps-mpi.%.4i.bin", (int)(qgps_time / qgps_time_step));
        return s;
}

int qgps_open() {
        MPI_File_open(QGPS_COMM_WORLD, qgps_output_filename(),
                        MPI_MODE_CREATE | MPI_MODE_WRONLY,
                        MPI_INFO_NULL, &qgps_output_file);

        MPI_File_set_view(qgps_output_file, qgps_current_real_block->y_begin * QGPS_NX,
                        MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL);
}

int qgps_close() {
        return MPI_File_close(&qgps_output_file);
}

int qgps_write() {
        static double *omega_real = NULL;

        if (!omega_real)
                omega_real = fftw_alloc_real(qgps_local_size * 2);

        qgps_dft_c2r(omega, omega_real);

        MPI_File_write(qgps_output_file, omega_real, qgps_local_size * 2,
                        MPI_DOUBLE, MPI_STATUS_IGNORE);
        return 0;
}
