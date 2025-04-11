#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[]) {

    int world_size, rank, num_spawns = 3;
    MPI_Comm intercomm;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int total_processes = world_size + world_size * num_spawns;
    if (total_processes > 8) {
        if (rank == 0) {
            fprintf(stderr, "Total processes (%d) exceed limit of 8\n", total_processes);
        }

        MPI_Finalize();

        return 0;
    }

    char *spawn_args[] = {NULL};
    MPI_Info info = MPI_INFO_NULL;

    int errcodes[num_spawns];
    int total_requested = num_spawns;

    MPI_Comm_spawn("slave", spawn_args, num_spawns, info, 
                  0, MPI_COMM_SELF, &intercomm, errcodes);

    for (int i = 0; i < num_spawns; i++) {
        if (errcodes[i] != MPI_SUCCESS) {
            fprintf(stderr, "Master %d: error spawning process %d\n", rank, i);
            MPI_Abort(MPI_COMM_WORLD, MPI_ERR_SPAWN);
        }
    }

    printf("Master %d: successfully spawned %d processes\n", rank, total_requested);

    int start_msg = 42 + rank;
    MPI_Send(&start_msg, 1, MPI_INT, 0, 0, intercomm);

    int final_msg;
    MPI_Recv(&final_msg, 1, MPI_INT, num_spawns-1, 0, intercomm, MPI_STATUS_IGNORE);
    printf("Master %d: received final message %d\n", rank, final_msg);

    MPI_Comm_free(&intercomm);
    MPI_Finalize();
    
    return 0;

}
