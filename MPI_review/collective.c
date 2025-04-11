#include <mpi.h>
#include <stdio.h>
#include <string.h>

int main(int argc, char *argv[]) {

    int rank, size;
    MPI_File fh;
    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    char buffer[16];
    snprintf(buffer, sizeof(buffer), "%d", rank);
    
    int len = strlen(buffer);

    memset(buffer + len, ' ', 16 - len - 1);
    buffer[15] = '\n';

    MPI_Offset offset = (size - 1 - rank) * 16;

    MPI_File_open(MPI_COMM_WORLD, "output.txt", 
                 MPI_MODE_CREATE | MPI_MODE_WRONLY, 
                 MPI_INFO_NULL, &fh);

    MPI_File_write_at_all(fh, offset, buffer, 16, MPI_CHAR, &status);

    MPI_File_close(&fh);

    if (rank == 0) {
        printf("Ranks written to output.txt in descending order\n");
    }

    MPI_Finalize();

    return 0;

}
