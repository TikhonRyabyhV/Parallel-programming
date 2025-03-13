#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>

#include "macro.h"

int main (int argc, char* argv[]) {

	int comm_size = 1, my_rank = 0;

	char* endptr = NULL;
	long long N = strtoll(argv[1], &endptr, 10);

	check(MPI_Init (&argc, &argv), "Cannot initialize MPI!")

	check(MPI_Comm_size (MPI_COMM_WORLD, &comm_size     ),
			"Cannot get communicator size!"     )
	check(MPI_Comm_rank (MPI_COMM_WORLD, &my_rank       ),
			"Cannot get rnk of current process!")

	long long comp_len = N / comm_size;

	long long start  = my_rank * comp_len + 1;
	long long finish = my_rank == (comm_size - 1) ? N + 1 : (my_rank + 1) * comp_len + 1; 

	double result = 0;

	for (long long i  = start; i < finish; ++i) {
		result += (1.0 / i);
	}

	if (my_rank != 0) {
		check(MPI_Send (&result, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD),
				                      "Cannot send message!")
	}

	else {
		double add = 0;

		for (int i = 1; i < comm_size; ++i) {
			check(MPI_Recv (&add, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE),
					                                   "Cannot receive message!")
			result += add;
		}

		printf ("Sum of number series 1/n (N = %lld): %.20lf.\n", N, result);
	}

	check(MPI_Finalize(), "Bad end of MPI!")

	return 0;

}
