#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#define ISIZE 5000
#define JSIZE 5000

int main(int argc, char **argv) {
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (8 % size != 0) {
        if (rank == 0) {
            printf("Error: Number of processes must be a divisor of 8 (1, 2, 4, 8) for this optimization.\n");
        }
        MPI_Finalize();
        return 1;
    }

    int local_cols_count = JSIZE / size;

    double *local_a = (double *)malloc(local_cols_count * ISIZE * sizeof(double));

    for (int j_loc = 0; j_loc < local_cols_count; j_loc++) {
        int global_j = j_loc * size + rank;
        for (int i = 0; i < ISIZE; i++) {
            local_a[j_loc * ISIZE + i] = 10 * i + global_j;
        }
    }

    // start timer
    double tstart = MPI_Wtime();

    // external cycle (j, columns)
    for (int j_loc = 0; j_loc < local_cols_count; j_loc++) {
        int global_j = j_loc * size + rank;

        if (global_j < 8) continue;

        int stride = 8 / size;
        int prev_j_loc = j_loc - stride;

        // internal cycle (i, lines)
        for (int i = 1; i < ISIZE; i++) {
            double prev_val = local_a[prev_j_loc * ISIZE + (i - 1)];
            
            local_a[j_loc * ISIZE + i] = sin(5 * prev_val);
        }
    }

    // end timer
    double tend = MPI_Wtime();

    // collect data
    double *full_data = NULL;
    if (rank == 0) {
        full_data = (double *)malloc(ISIZE * JSIZE * sizeof(double));
    }

    MPI_Gather(local_a, local_cols_count * ISIZE, MPI_DOUBLE,
               full_data, local_cols_count * ISIZE, MPI_DOUBLE,
               0, MPI_COMM_WORLD);

    if (rank == 0) {
        // exec time info
        printf("Task 1e finished. Time: %f\n", tend - tstart);

        FILE *ff = fopen("result_mpi_opt.txt", "w");
        for (int i = 0; i < ISIZE; i++) {
            for (int j = 0; j < JSIZE; j++) {
                int owner_rank = j % size;
                int owner_col_idx = j / size;
                
                long rank_offset = (long)owner_rank * (local_cols_count * ISIZE);
                long local_offset = (long)owner_col_idx * ISIZE + i;
                
                fprintf(ff, "%f ", full_data[rank_offset + local_offset]);
            }
            fprintf(ff, "\n");
        }
        fclose(ff);
        free(full_data);
    }

    free(local_a);
    MPI_Finalize();
    return 0;
}
