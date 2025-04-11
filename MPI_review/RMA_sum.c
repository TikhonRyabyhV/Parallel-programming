#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[]) {

    int rank, size;
    double local_sum = 0.0, global_sum = 0.0;
    long long N;
    MPI_Win win;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc != 2) {
        if (rank == 0) {
            fprintf(stderr, "Usage: %s <N>\n", argv[0]);
        }
        MPI_Finalize();

        return 0;
    }

    N = atoll(argv[1]);
    if (N <= 0) {
        if (rank == 0) {
            fprintf(stderr, "N must be positive\n");
        }
        MPI_Finalize();

        return 0;
    }

    double *sums;
    if (rank == 0) {
        sums = (double *)calloc(size, sizeof(double));
    } 
    
    else {
        sums = NULL;
    }

    MPI_Win_create(sums, rank == 0 ? size * sizeof(double) : 0, sizeof(double), 
                   MPI_INFO_NULL, MPI_COMM_WORLD, &win);

    long long base_chunk = N / size;
    long long remainder = N % size;
    long long start, end;

    if (rank < remainder) {
        start = rank * (base_chunk + 1) + 1;
        end = start + base_chunk;
    } 
    
    else {
        start = rank * base_chunk + remainder + 1;
        end = start + base_chunk - 1;
    }


    if (end > N) {
	end = N;
    }

    for (long long i = start; i <= end; i++) {
        local_sum += 1.0 / i;
    }

    // for debug:
    //printf("Rank %d: start = %ld, end = %ld, local_sum = %f\n", rank, start, end, local_sum);

    MPI_Win_fence(0, win);

    if (rank != 0) {
        MPI_Put(&local_sum, 1, MPI_DOUBLE, 0, rank, 1, MPI_DOUBLE, win);
    } 
    
    else {
        sums[0] = local_sum;
    }

    MPI_Win_fence(0, win);

    if (rank == 0) {
        for (int i = 0; i < size; i++) {
            global_sum += sums[i];
        }
        
	printf("Sum from 1 to %Ld = %f\n", N, global_sum);
        
	free(sums);
    }

    MPI_Win_free(&win);
    MPI_Finalize();

    return 0;

}
