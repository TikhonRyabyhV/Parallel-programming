#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#define ISIZE 5000
#define JSIZE 5000

int main(int argc, char **argv)
{
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int i_start = rank * (ISIZE / size);
    int i_end = (rank + 1) * (ISIZE / size);
    if (rank == size - 1) i_end = ISIZE;
    int i_count = i_end - i_start;

    double *a_local = (double *)malloc(i_count * JSIZE * sizeof(double));

    int i, j;
    for (i = 0; i < i_count; i++) {
        int global_i = i_start + i;
        for (j = 0; j < JSIZE; j++) {
            a_local[i * JSIZE + j] = 10 * global_i + j;
        }
    }

    double t_start = MPI_Wtime();

    for (i = 0; i < i_count; i++) {
        for (j = 0; j < JSIZE; j++) {
            a_local[i * JSIZE + j] = sin(2 * a_local[i * JSIZE + j]);
        }
    }

    double t_end = MPI_Wtime();

    if (rank == 0) {
        printf("Reference MPI Time: %f\n", t_end - t_start);
        
        FILE *ff = fopen("result_ref_mpi.txt", "w");
        
        for (i = 0; i < i_count; i++) {
            for (j = 0; j < JSIZE; j++) {
                fprintf(ff, "%f ", a_local[i * JSIZE + j]);
            }
            fprintf(ff, "\n");
        }

        double *recv_buf = (double *)malloc(((ISIZE / size) + 1) * JSIZE * sizeof(double));
        MPI_Status status;

        for (int r = 1; r < size; r++) {
            int r_i_start = r * (ISIZE / size);
            int r_i_end = (r + 1) * (ISIZE / size);
            if (r == size - 1) r_i_end = ISIZE;
            int r_count = r_i_end - r_i_start;

            MPI_Recv(recv_buf, r_count * JSIZE, MPI_DOUBLE, r, 0, MPI_COMM_WORLD, &status);

            for (i = 0; i < r_count; i++) {
                for (j = 0; j < JSIZE; j++) {
                    fprintf(ff, "%f ", recv_buf[i * JSIZE + j]);
                }
                fprintf(ff, "\n");
            }
        }
        free(recv_buf);
        fclose(ff);
    } else {
        MPI_Send(a_local, i_count * JSIZE, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }

    free(a_local);
    MPI_Finalize();
    return 0;
}

