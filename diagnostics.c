#include <stdio.h>
#include "qgps-mpi.h"
#include "diagnostics.h"

double qgps_total_energy = 0,
       qgps_enstrophy = 0;

complex *temp = NULL;

int qgps_update_enstrophy();
int qgps_update_total_energy();

int qgps_diagnostics_update() {
        if (!temp)
                temp = fftw_alloc_complex(qgps_local_size);

        qgps_update_total_energy();
        qgps_update_enstrophy();

        return 0;
}

int qgps_update_total_energy() {
        complex *const specific_vorticity = temp;

        for (int i = 1; i < qgps_local_size; i++)
                specific_vorticity[i] = omega[i] / qgps_k[i]; 

        qgps_total_energy = l2_norm_squared(specific_vorticity) / 2;

        return 0;
}
int qgps_update_enstrophy() {
        qgps_enstrophy = l2_norm_squared(omega);

        return 0;
}
