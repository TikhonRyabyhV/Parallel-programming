#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#include "macro.h"

int main (int argc, char* argv[]) {

	int comm_size = 1, my_rank = 0;

	check(MPI_Init (&argc, &argv), "Cannot initialize MPI!")

	check(MPI_Comm_size (MPI_COMM_WORLD, &comm_size),
		       	"Cannot get communicator size!")
	check(MPI_Comm_rank (MPI_COMM_WORLD, &my_rank  ),
		  "Cannot get rank of current process!")

	int value = 0;

	if (my_rank == 0) {
		printf ("Initial value: %d.\n", value);

		check(MPI_Send (&value, 1, MPI_INT, 1, 0, MPI_COMM_WORLD),
				            "Cannot send first message!")
	}

	int prev_rank = 0;
	if (my_rank == 0) {
		prev_rank = comm_size - 1;
	}

	else {
		prev_rank = my_rank - 1;
	}	

	check(MPI_Recv (&value, 1, MPI_INT, prev_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE),
			                                          "Cannot receive message!")
	value += 1;

	printf ("Current value: %d (my_rank = %d).\n", value, my_rank);

	if (my_rank != 0) {
		check(MPI_Send (&value, 1, MPI_INT, ((my_rank + 1) % comm_size), 0, MPI_COMM_WORLD),
				                                            "Cannot send message!")
	}

	check(MPI_Finalize(), "Bad end of MPI!")
	
	return 0;

}
