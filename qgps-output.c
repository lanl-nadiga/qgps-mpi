#include <stdlib.h>
#include <mpi.h>
#include "qgps-mpi.h"

int qgps_output() {
        return qgps_output_complex(qgps_output_filename(), omega);
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

MPI_File qgps_output_open_file(char *filename) {
        MPI_File file = NULL;

        if (MPI_File_open(QGPS_COMM_WORLD, filename,
                        MPI_MODE_CREATE | MPI_MODE_WRONLY,
                        MPI_INFO_NULL, &file))
                return NULL;

        if (MPI_File_set_view(file, qgps_current_real_block->x_begin * qgps_nx * sizeof(double),
                        MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL))
                return NULL;

        return file;
}

int qgps_output_close_file(MPI_File *file) {
        return MPI_File_close(file);
}

int qgps_output_real(char *filename, double *data) {
        MPI_File file = qgps_output_open_file(filename);

        if (!file)
                return 1;

        if (qgps_output_write_real(file, data))
                return 1;

        return qgps_output_close_file(&file);
}
int qgps_output_complex(char *filename, complex *data) {
        MPI_File file = qgps_output_open_file(filename);

        if (!file)
                return 1;

        if (qgps_output_write_complex(file, data))
                return 1;

        return qgps_output_close_file(&file);
}
int qgps_output_write_complex(MPI_File file, complex *data) {
        static double *data_real = NULL;

        if (!data_real)
                data_real = fftw_alloc_real(qgps_local_size * 2);

        if (!data_real)
                return 1;

        if (qgps_dft_c2r(data, data_real))
                return 1;

        return qgps_output_write_real(file, data_real);
}
int qgps_output_write_real(MPI_File file, double *data) {
        if (qgps_transpose_r(data))
                return 1;

        qgps_block_t *b = qgps_current_real_block;

        if (!b)
                return 1;

        int pad = 2 - qgps_ny % 2;
        for (int j = b->x_begin; j < b->x_end; j++) {
                int idx = (j - b->x_begin) * (qgps_nx + pad);

                if (MPI_File_write(file, &data[idx], qgps_nx,
                                MPI_DOUBLE, MPI_STATUS_IGNORE))
                        return 1;
        }
        return 0;
}
