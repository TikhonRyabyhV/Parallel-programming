//
// (i, j) depends on (i - 1, j - 8) => D = (1, 8), True Dependence,
// so we use column decomposition
//
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#define ISIZE 5000
#define JSIZE 5000
#define GHOST 8

int main(int argc, char **argv)
{
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int j_start = rank * (JSIZE / size);
    int j_end = (rank + 1) * (JSIZE / size);
    if (rank == size - 1) j_end = JSIZE;
    int j_count = j_end - j_start;

    double *a_local = (double *)malloc(ISIZE * j_count * sizeof(double));
    
    int i, j;
    for (i = 0; i < ISIZE; i++) {
        for (j = j_start; j < j_end; j++) {
            a_local[i * j_count + (j - j_start)] = 10 * i + j;
        }
    }

    double *ghost_recv = (double*)malloc(GHOST * sizeof(double));
    double *ghost_send = (double*)malloc(GHOST * sizeof(double));

    double t_start = MPI_Wtime();

    for (i = 1; i < ISIZE; i++) {
        if (rank < size - 1) {
             for (int k = 0; k < GHOST; k++) {
                 ghost_send[k] = a_local[(i-1) * j_count + (j_count - GHOST + k)];
             }
             MPI_Send(ghost_send, GHOST, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
        }
        
        if (rank > 0) {
             MPI_Recv(ghost_recv, GHOST, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        for (j = j_start; j < j_end; j++) {
            int local_j = j - j_start;
            
            if (j >= 8) {
                double val_prev;
                if (local_j >= 8) {
                    val_prev = a_local[(i-1) * j_count + (local_j - 8)];
                } else {
                    val_prev = ghost_recv[local_j]; 
                }
                a_local[i * j_count + local_j] = sin(5 * val_prev);
            }
        }
    }

    double t_end = MPI_Wtime();
    
    if (rank == 0) {
        FILE *ff = fopen("result.txt", "w");
        double *row_buffer = (double*)malloc(JSIZE * sizeof(double));
        
        for (i = 0; i < ISIZE; i++) {
            for (int k = 0; k < j_count; k++) {
                row_buffer[k] = a_local[i * j_count + k];
            }
            
            int current_pos = j_count;
            for (int r = 1; r < size; r++) {
                int r_j_start = r * (JSIZE / size);
                int r_j_end = (r + 1) * (JSIZE / size);
                if (r == size - 1) r_j_end = JSIZE;
                int r_count = r_j_end - r_j_start;
                
                MPI_Recv(row_buffer + current_pos, r_count, MPI_DOUBLE, r, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                current_pos += r_count;
            }
            
            for (j = 0; j < JSIZE; j++) {
                fprintf(ff, "%f ", row_buffer[j]);
            }
            fprintf(ff, "\n");
        }
        fclose(ff);
        free(row_buffer);
        printf("Task 1e finished. Time: %f\n", t_end - t_start);
        
    } else {
        for (i = 0; i < ISIZE; i++) {
             MPI_Send(a_local + i * j_count, j_count, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
        }
    }

    free(a_local);
    free(ghost_recv);
    free(ghost_send);
    MPI_Finalize();
    return 0;
}

