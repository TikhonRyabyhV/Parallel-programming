#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "LongDouble.h"

#include "macro.h"


unsigned int comp_N (unsigned int accuracy) {

	double up   = 1.0           ;
	double down = ((double) accuracy) + 1.0;

	while (up - down > 0.5) {
        
		double middle = (up + down) / 2;
		double val    = middle * log10(middle);

		if (val < accuracy) {
			down = middle;
		}
		
		else {
			up = middle;
		}

	}

	return 2 * ((unsigned int) ((up + down) / 2));

}


int main (int argc, char* argv[]) {

	int comm_size = 1, my_rank = 0;

	char* endptr = NULL;
	unsigned int accuracy = atoi(argv[1]);

	check(MPI_Init (&argc, &argv), "Cannot initialize MPI!")

	check(MPI_Comm_size (MPI_COMM_WORLD, &comm_size     ),
			"Cannot get communicator size!"     )
	check(MPI_Comm_rank (MPI_COMM_WORLD, &my_rank       ),
			"Cannot get rank of current process!")


	unsigned int N = comp_N (accuracy);

	LongDouble data("1.0", accuracy + 1);
	LongDouble item("1.0", accuracy + 1);

	unsigned int comp_len = N / comm_size;

	unsigned int start  = my_rank * comp_len + 1;
	unsigned int finish = my_rank == (comm_size - 1) ? N + 1 : (my_rank + 1) * comp_len + 1;

	for (unsigned int i  = start; i < finish; ++i) {
		item /= i;
		
		data += item;
	}

	LongDouble add("1.0", accuracy + 1);
	LongDouble   e("0.0", accuracy + 1);
	

	if(my_rank < (comm_size - 1)) {
		check(MPI_Recv (&add, 1, MPI_LONG_DOUBLE, my_rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE),
					                                   "Cannot receive message!")

		data += (item * add);
	}

	if (my_rank != 0) {
		check(MPI_Send (&data, 1, MPI_LONG_DOUBLE, my_rank - 1, 0, MPI_COMM_WORLD),
				                      "Cannot send message!")
	}

	else {
		e += data;
	}

	if(my_rank == 0) {
		std::cout << e << std::endl;
	}

	check(MPI_Finalize(), "Bad end of MPI!")

	return 0;

}
