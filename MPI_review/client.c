#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#define MAX_PORT_NAME 256

int main(int argc, char *argv[]) {
    int rank, size;
    MPI_Comm server_comm;
    char port_name[MAX_PORT_NAME];
    int data_to_send = 123;
    int received_response;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    sleep(2); // Задержка для создания port.txt

    FILE *port_file = fopen("port.txt", "r");
    if (!port_file) {
        fprintf(stderr, "Client: Cannot open port file\n");
        MPI_Finalize();
        return 1;
    }
    if (fgets(port_name, MAX_PORT_NAME, port_file) == NULL) {
        fprintf(stderr, "Client: Cannot read port name\n");
        fclose(port_file);
        MPI_Finalize();
        return 1;
    }
    port_name[strcspn(port_name, "\n")] = 0;
    fclose(port_file);

    printf("Client: Connecting to port %s\n", port_name);
    if (MPI_Comm_connect(port_name, MPI_INFO_NULL, 0, MPI_COMM_WORLD, &server_comm) != MPI_SUCCESS) {
        fprintf(stderr, "Client: Failed to connect to server\n");
        MPI_Finalize();
        return 1;
    }
    printf("Client: Connected to server\n");

    MPI_Send(&data_to_send, 1, MPI_INT, 0, 0, server_comm);
    printf("Client: Sent %d to server\n", data_to_send);

    MPI_Recv(&received_response, 1, MPI_INT, 0, 0, server_comm, MPI_STATUS_IGNORE);
    printf("Client: Received %d from server\n", received_response);

    MPI_Comm_disconnect(&server_comm);
    printf("Client: Disconnected from server\n");

    MPI_Finalize();
    return 0;
}
