#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[]) {

    int rank;
    MPI_Comm server_comm;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Comm_get_parent(&server_comm);

    if (server_comm == MPI_COMM_NULL) {
        fprintf(stderr, "This is not a spawned process!\n");
        MPI_Finalize();

        return 0;
    }

    int recv_msg;
    MPI_Recv(&recv_msg, 1, MPI_INT, 0, 0, server_comm, MPI_STATUS_IGNORE);
    printf("Client %d: received %d from server\n", rank, recv_msg);

    int send_msg = recv_msg + 50;
    MPI_Send(&send_msg, 1, MPI_INT, 0, 0, server_comm);
    printf("Client %d: sent %d to server\n", rank, send_msg);

    MPI_Comm_free(&server_comm);
    MPI_Finalize();
    
    return 0;

}
