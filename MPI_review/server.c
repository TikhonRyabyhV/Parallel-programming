#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[]) {

    int world_size, rank;
    MPI_Comm client_comm;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (world_size != 1) {
        if (rank == 0) {
            fprintf(stderr, "Server should run with one process\n");
        }

        MPI_Finalize();

        return 0;
    }

    if (argc != 2) {
        fprintf(stderr, "Usage: %s <number_of_clients>\n", argv[0]);
        MPI_Finalize();

        return 0;
    }

    int num_clients = atoi(argv[1]);
    if (num_clients <= 0 || num_clients > 7) {
        fprintf(stderr, "Number of clients must be between 1 and 7\n");
        MPI_Finalize();

        return 0;
    }

    char *spawn_args[] = {NULL};
    MPI_Info info = MPI_INFO_NULL;
    int *errcodes = malloc(num_clients * sizeof(int));

    MPI_Comm_spawn("client", spawn_args, num_clients, info, 0, MPI_COMM_SELF, &client_comm, errcodes);

    for (int i = 0; i < num_clients; i++) {
        if (errcodes[i] != MPI_SUCCESS) {
            fprintf(stderr, "Error spawning client %d\n", i);
            MPI_Abort(MPI_COMM_WORLD, MPI_ERR_SPAWN);
        }
    }

    printf("Server: successfully spawned %d clients\n", num_clients);

    for (int i = 0; i < num_clients; i++) {
        int send_msg = 100 + i;
        MPI_Send(&send_msg, 1, MPI_INT, i, 0, client_comm);
        printf("Server: sent %d to client %d\n", send_msg, i);
    }

    for (int i = 0; i < num_clients; i++) {
        int recv_msg;
        MPI_Recv(&recv_msg, 1, MPI_INT, i, 0, client_comm, MPI_STATUS_IGNORE);
        printf("Server: received %d from client %d\n", recv_msg, i);
    }

    free(errcodes);
    MPI_Comm_free(&client_comm);
    MPI_Finalize();
    
    return 0;

}
