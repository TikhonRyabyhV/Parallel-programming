#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[]) {

    int parent_size, rank, size;
    MPI_Comm parent_comm;

    MPI_Init(&argc, &argv);
    
    MPI_Comm_get_parent(&parent_comm);
    
    if (parent_comm == MPI_COMM_NULL) {
        fprintf(stderr, "This is not a child process!\n");
        MPI_Finalize();

        return 0;
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    MPI_Comm_remote_size(parent_comm, &parent_size);

    int msg;
    if (rank == 0) {
        MPI_Recv(&msg, 1, MPI_INT, 0, 0, parent_comm, MPI_STATUS_IGNORE);
        printf("Slave %d: received %d from master\n", rank, msg);
    } 
    
    else {
        MPI_Recv(&msg, 1, MPI_INT, rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        printf("Slave %d: received %d from %d\n", rank, msg, rank-1);
    }

    msg += rank;
    
    if (rank == size-1) {
        MPI_Send(&msg, 1, MPI_INT, 0, 0, parent_comm);
        printf("Slave %d: sent %d to master\n", rank, msg);
    } 
    
    else {
        MPI_Send(&msg, 1, MPI_INT, rank+1, 0, MPI_COMM_WORLD);
        printf("Slave %d: sent %d to process %d\n", rank, msg, rank+1);
    }

    MPI_Comm_free(&parent_comm);
    MPI_Finalize();

    return 0;

}
