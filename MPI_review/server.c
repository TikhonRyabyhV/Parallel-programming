#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char *argv[]) {
    int rank, size;
    MPI_Comm client_comm;
    MPI_Status status;
    char port_name[MPI_MAX_PORT_NAME];
    int received_data;
    int response = 42;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (MPI_Open_port(MPI_INFO_NULL, port_name) != MPI_SUCCESS) {
        fprintf(stderr, "Server: Failed to open port\n");
        MPI_Finalize();
        return 1;
    }
    printf("Server opened port: %s\n", port_name);

    FILE *port_file = fopen("port.txt", "w");
    if (!port_file) {
        fprintf(stderr, "Server: Cannot open port file\n");
        MPI_Close_port(port_name);
        MPI_Finalize();
        return 1;
    }
    fprintf(port_file, "%s\n", port_name);
    fclose(port_file);

    printf("Server: Waiting for client...\n");
    if (MPI_Comm_accept(port_name, MPI_INFO_NULL, 0, MPI_COMM_WORLD, &client_comm) != MPI_SUCCESS) {
        fprintf(stderr, "Server: Failed to accept client connection\n");
        MPI_Close_port(port_name);
        MPI_Finalize();
        return 1;
    }
    printf("Server: Client connected\n");

    MPI_Recv(&received_data, 1, MPI_INT, 0, 0, client_comm, &status);
    printf("Server: Received %d from client\n", received_data);

    MPI_Send(&response, 1, MPI_INT, 0, 0, client_comm);
    printf("Server: Sent response %d to client\n", response);

    MPI_Comm_disconnect(&client_comm);
    printf("Server: Client disconnected\n");

    MPI_Close_port(port_name);
    MPI_Finalize();
    return 0;
}
