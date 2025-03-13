#include <mpi.h>
#include <stdio.h>

#include "macro.h"

int main (int argc, char* argv[]) {

	int comm_size = 1, my_rank = 0;

	check(MPI_Init (&argc, &argv), "Cannot initialize MPI!")

	check(MPI_Comm_size (MPI_COMM_WORLD, &comm_size      ),
			"Cannot get communicator size!"      )
	check(MPI_Comm_rank (MPI_COMM_WORLD, &my_rank        ),
			"Cannot get rank of current process!")

	printf ("Communicator size = %d, my rank = %d\n", comm_size, my_rank);

	check(MPI_Finalize(), "Bad end of MPI!")
	
	return 0;

}	
