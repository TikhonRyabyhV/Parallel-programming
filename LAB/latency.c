#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>

#define NUM_MESSAGES 1000
#define MESSAGE_SIZE 100

int main(int argc, char *argv[]) {
    int rank, size;
    double start_time, end_time, total_time, avg_latency;
    char message[MESSAGE_SIZE];
    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (size != 2) {
        if (rank == 0) {
            fprintf(stderr, "This program requires exactly 2 processes\n");
        }
        MPI_Finalize();
        return 1;
    }

    for (int i = 0; i < MESSAGE_SIZE; i++) {
        message[i] = 0;
    }

    if (rank == 0) {
        start_time = MPI_Wtime();

        for (int i = 0; i < NUM_MESSAGES; i++) {
            MPI_Send(message, MESSAGE_SIZE, MPI_CHAR, 1, 0, MPI_COMM_WORLD);
            MPI_Recv(message, MESSAGE_SIZE, MPI_CHAR, 1, 0, MPI_COMM_WORLD, &status);
        }

        end_time = MPI_Wtime();
        total_time = end_time - start_time;

        avg_latency = (total_time / (2 * NUM_MESSAGES)) * 1e6; // В микросекундах

        printf("Average latency for %d-byte message: %.2f microseconds\n", 
               MESSAGE_SIZE, avg_latency);
    } else if (rank == 1) {
        for (int i = 0; i < NUM_MESSAGES; i++) {
            MPI_Recv(message, MESSAGE_SIZE, MPI_CHAR, 0, 0, MPI_COMM_WORLD, &status);
            MPI_Send(message, MESSAGE_SIZE, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
        }
    }

    MPI_Finalize();
    return 0;
}
