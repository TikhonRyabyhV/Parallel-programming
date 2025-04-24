#pragma once

#include <mpi.h>
#include <stdio.h>

#define check(func_res, message)                                          \
	if(func_res != MPI_SUCCESS) {                                     \
		printf ("%s\nError with code: %d.\n", message, func_res); \
		MPI_Abort (MPI_COMM_WORLD, func_res);                     \
	}
