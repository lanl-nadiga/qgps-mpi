#include "qgps-mpi.h"
#include "step.h"
#include <string.h>
#include <stdlib.h>

const int QGPS_NX = 64, QGPS_NY = 64;

MPI_Comm QGPS_COMM_WORLD;
int  qgps_current_task, qgps_master_task = 0, qgps_number_tasks;

qgps_block_t * qgps_real_blocks                 = NULL;
qgps_block_t * qgps_complex_blocks              = NULL;
qgps_block_t * qgps_current_real_block          = NULL;
qgps_block_t * qgps_current_complex_block       = NULL;

ptrdiff_t qgps_local_size;

int qgps_initialize_mpi(int argc, char **argv);
int qgps_initialize_blocks();
int qgps_cleanup_mpi();
int qgps_broadcast_block(qgps_block_t *block, int src_task);
int qgps_initialize_fftw();

int qgps_dft_c2r(const complex *in, double *out) {
        static complex *temporary = NULL;
        if (!temporary)
                temporary = fftw_alloc_complex(qgps_local_size);

        static fftw_plan plan = NULL;
        if (!plan)
                plan = fftw_mpi_plan_dft_c2r_2d(QGPS_NX, QGPS_NY,
                                                temporary, out,
                                                QGPS_COMM_WORLD,
                                                FFTW_MEASURE);

        memcpy(temporary, in, sizeof(complex) * qgps_local_size);

        fftw_mpi_execute_dft_c2r(plan, temporary, out);
        return 0;
}

int qgps_dft_r2c(const double *in, complex *out) {
        static double *temporary = NULL;
        if (!temporary)
                temporary = fftw_alloc_real(qgps_local_size*2);

        static fftw_plan plan = NULL;
        if (!plan)
                plan = fftw_mpi_plan_dft_r2c_2d(QGPS_NX, QGPS_NY,
                                                temporary, out,
                                                QGPS_COMM_WORLD,
                                                FFTW_MEASURE);

        memcpy(temporary, in, sizeof(double) * qgps_local_size * 2);

        fftw_mpi_execute_dft_r2c(plan, temporary, out);
        return 0;
}

int qgps_initialize(int argc, char **argv) {
//        if (qgps_configure(argc, argv))
//                return 1;
        if (qgps_initialize_mpi(argc, argv))
                return 1;
        if (qgps_initialize_blocks())
                return 1;
        if (qgps_initialize_fftw())
                return 1;
        if (qgps_step_init())
                return 1;

        return 0;
}

int qgps_cleanup() {
        if (qgps_step_free())
                return 1;
        if (qgps_cleanup_mpi())
                return 1;

        return 0;
}

void qgps_exit() {
        qgps_cleanup();
        exit(EXIT_FAILURE);
}

int qgps_initialize_mpi(int argc, char **argv) {
        MPI_Init(&argc, &argv);
        MPI_Comm_dup(MPI_COMM_WORLD, &QGPS_COMM_WORLD);

        MPI_Comm_rank(QGPS_COMM_WORLD, &qgps_current_task);
        MPI_Comm_size(QGPS_COMM_WORLD, &qgps_number_tasks);

        return 0;
}

int qgps_initialize_blocks(int argc, char **argv) {
        qgps_real_blocks = calloc(qgps_number_tasks, sizeof(qgps_block_t));
        if (!qgps_real_blocks)
                return 1;

        qgps_complex_blocks = calloc(qgps_number_tasks, sizeof(qgps_block_t));
        if (!qgps_complex_blocks)
                return 1;

        return 0;
}

int qgps_cleanup_mpi() {
        MPI_Finalize();

        return 0;
}

int qgps_broadcast_block(qgps_block_t *block, int src_task) {
        MPI_Bcast(block, sizeof(qgps_block_t), MPI_BYTE, src_task,
                                                QGPS_COMM_WORLD);
        MPI_Barrier(QGPS_COMM_WORLD);

        return 0;
}

int qgps_initialize_fftw() {
        ptrdiff_t alloc_local, local_n0, local_0_start;
        qgps_block_t* block;

        fftw_mpi_init();

        qgps_local_size = fftw_mpi_local_size_2d(QGPS_NX, QGPS_NY/2 + 1,
                                                        QGPS_COMM_WORLD,
                                                        &local_n0,
                                                        &local_0_start);

        qgps_current_real_block    = &qgps_real_blocks[qgps_current_task];
        qgps_current_complex_block = &qgps_complex_blocks[qgps_current_task];

        // init block information for this task
 
        /*
         * REAL BLOCKS:
         *
         * Each real block contains a sequence of consecutive *columns* of the
         * real space, where each individual column is terminated by cells of
         * padding (junk data).  The padding takes up 1 cell if NY is odd and
         * two cells if NY is even.
         *
         * For example, if r(i,j) represents the i'th column, j'th row value of
         * the real global array of size NX by NY (even), and the specific block
         * starts at the k'th column, then the data layout of the local array is
         *
         * local(0)     = r(k  ,0)
         * local(1)     = r(k  ,1)
         * local(2)     = r(k  ,2)
         * ...
         * local(NY-1)  = r(k  ,NY-1)
         * local(NY)    = pad
         * local(NY+1)  = pad
         * local(NY+2)  = r(k+1,0)
         * local(NY+3)  = r(k+1,1)
         * local(NY+4)  = r(k+1,2)
         * ...
         * etc.
         *
         * More specifically, for 0 <= j < NY and i less than the number of
         * rows contained within the current local block
         *
         * local((NY+2)*i+j) = r(k+i,j)
         *
         * (important: the columns are iterated first in this scheme!)
         *
         */
        
        block = qgps_current_real_block;
        block->id       = qgps_current_task;
        block->x_begin  = local_0_start;
        block->x_end    = block->x_begin + local_n0;
        block->x_length = qgps_local_size / QGPS_NY;
        block->y_begin  = 0;
        block->y_end    = QGPS_NY;
        block->y_length = QGPS_NY;
        block->size     = 2*qgps_local_size;

        /*
         * COMPLEX BLOCKS:
         *
         * Each complex block contains a sequence of consecutive *rows* of the
         * complex space.  Unlike with the real blocks, there is no padding.
         * The complex space itself contains only NX by NY/2 + 1 complex
         * numbers, exploiting the fact that Fourier transform F of a real
         * function f satisfies the relation
         * 
         * F(i,j) = conjugate ( F(NX-i,NY-j) )
         *
         * For example, if c(i,j) represents the i'th column, j'th row value of
         * the complex global array of size NX by NY/2 + 1, and the specific
         * block starts at the k'th row, then the data layout of the local
         * array is
         *
         * local(0)     = c(0  ,k)
         * local(1)     = c(1  ,k)
         * local(2)     = c(2  ,k)
         * ...
         * local(NX-1)  = c(NX-1 ,k)
         * local(NX)    = c(0    ,k+1)
         * local(NX+1)  = c(1    ,k+1)
         * local(NX+2)  = c(2    ,k+1)
         * ...
         * etc.
         *
         * More specifically, for 0 <= j < NY and i less than the number of
         * rows contained within the current local block
         *
         * local(NX*j+i) = c(i,k+j)
         *
         * (important: the number of complex values per column is one more than
         *             half the number of reals per column! )
         *
         */
        block = qgps_current_complex_block;
        block->id       = qgps_current_task;
        block->x_begin  = 0;
        block->x_end    = QGPS_NX;
        block->x_length = QGPS_NX;
        block->y_begin  = local_0_start;
        block->y_end    = block->y_begin + local_n0 / 2 + 1;
        block->y_length = qgps_local_size / QGPS_NX;
        block->size     = qgps_local_size;

        MPI_Barrier(QGPS_COMM_WORLD);

        // broadcast block information
        for (int task = 0; task < qgps_number_tasks; task++) {
                qgps_broadcast_block(&(qgps_real_blocks[task]), task);
        }

        for (int task = 0; task < qgps_number_tasks; task++) {
                qgps_broadcast_block(&(qgps_complex_blocks[task]), task);
        }

        return 0;
}

