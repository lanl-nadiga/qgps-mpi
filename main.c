#include "qgps-mpi.h"
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv) {
        if (qgps_initialize(argc, argv)) {
                fprintf(stderr, "failed to initialize qgps\n");
                return EXIT_FAILURE;
        }

        qgps_step();

        if (qgps_cleanup()) {
                fprintf(stderr, "failed to cleanup qgps\n");
                return EXIT_FAILURE;
        }
        return 0;
}

