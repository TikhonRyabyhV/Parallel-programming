#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <string.h>

#define OUTPUT_FILE "./results/parconveq.bin"

const double T = 100.0;    // 0 <= t <= T
const double X = 100.0;    // 0 <= x <= X
const double tau = 0.007;
const double h = 0.05;
const double c = 1.35;     // should be positive
const int K = (int)(T / tau);
const int M = (int)(X / h);

double phi(double x) {
    return sin(x / 5.0);
}

double psi(double t) {
    return 2.0 * sin(t);
}

double f(double t, double x) {
    int s = ((int)x) % 100;
    if (s < 40 || s > 60)
        return 0.0;
    return 5.0 * sin(x);
}

double solveByCornerScheme(double **u, int k, int m) {
    k--;
    return u[k][m] + tau * (f(k * tau, m * h) - c * (u[k][m] - u[k][m - 1]) / h);
}

double solveByCrossScheme(double **u, int k, int m) {
    k--;
    return u[k - 1][m] + 2.0 * tau * (f(k * tau, m * h) - c * (u[k][m + 1] - u[k][m - 1]) / (2.0 * h));
}

int main(int argc, char **argv) {
    int err = MPI_Init(&argc, &argv);
    if (err != MPI_SUCCESS) {
        printf("Bad MPI initialization\n");
        return 1;
    }

    int rank;
    err = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (err != MPI_SUCCESS) {
        printf("Can't get rank of this process\n");
        return 1;
    }

    int world_size;
    err = MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    if (err != MPI_SUCCESS) {
        printf("Can't get total count of processes\n");
        return 1;
    }

    int buffer_size = sizeof(double) * 4 + MPI_BSEND_OVERHEAD;
    char *buffer = (char *)malloc(buffer_size);
    err = MPI_Buffer_attach(buffer, buffer_size);
    if (err != MPI_SUCCESS) {
        printf("Can't attach buffer for Bsend\n");
        free(buffer);
        return 1;
    }

    int measurement = (argc > 1) && (strcmp(argv[1], "-m") == 0);

    // Allocate grid dynamically
    double **u = (double **)malloc((K + 1) * sizeof(double *));
    for (int i = 0; i <= K; i++) {
        u[i] = (double *)malloc((M + 1) * sizeof(double));
    }

    double start = MPI_Wtime();

    int partition_size = M / world_size;
    int last = world_size - 1;
    int min_m = (rank == 0) ? 1 : (rank * partition_size);
    int max_m = (rank == last) ? (M + 1) : ((rank + 1) * partition_size);

    // Initial conditions
    for (int m = min_m - 1; m < max_m; m++) {
        u[0][m] = phi(m * h);
    }

    // Boundary conditions on left side
    for (int k = 1; k <= K; k++) {
        if (rank == 0) {
            u[k][0] = psi(k * tau);
        }
    }

    // First layer by corner scheme
    for (int m = min_m; m < max_m; m++) {
        u[1][m] = solveByCornerScheme(u, 1, m);
    }

    if (rank > 0) {
        err = MPI_Bsend(&u[1][min_m - 1], 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
        if (err != MPI_SUCCESS) {
            printf("Can't send information to left process\n");
            return 1;
        }
    }
    if (rank < last) {
        err = MPI_Bsend(&u[1][max_m], 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
        if (err != MPI_SUCCESS) {
            printf("Can't send information to right process\n");
            return 1;
        }
    }

    MPI_Request left_request, right_request;
    for (int k = 2; k <= K; k++) {
        if (rank > 0) {
            err = MPI_Irecv(&u[k - 1][min_m - 1], 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &left_request);
            if (err != MPI_SUCCESS) {
                printf("Can't receive information from left process\n");
                return 1;
            }
        }
        if (rank < last) {
            err = MPI_Irecv(&u[k - 1][max_m], 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, &right_request);
            if (err != MPI_SUCCESS) {
                printf("Can't receive information from right process\n");
                return 1;
            }
        }

        for (int m = min_m + 1; m < max_m - 1; m++) {
            u[k][m] = solveByCrossScheme(u, k, m);
        }

        if (rank > 0) {
            err = MPI_Wait(&left_request, MPI_STATUS_IGNORE);
            if (err != MPI_SUCCESS) {
                printf("Can't receive information from left process\n");
                return 1;
            }
        }
        if (rank < last) {
            err = MPI_Wait(&right_request, MPI_STATUS_IGNORE);
            if (err != MPI_SUCCESS) {
                printf("Can't receive information from right process\n");
                return 1;
            }
        }

        u[k][min_m] = solveByCrossScheme(u, k, min_m);
        if (rank == last) {
            u[k][M] = solveByCornerScheme(u, k, M);
        } else {
            u[k][max_m - 1] = solveByCrossScheme(u, k, max_m - 1);
        }

        if (rank > 0) {
            err = MPI_Bsend(&u[k - 1][min_m], 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
            if (err != MPI_SUCCESS) {
                printf("Can't send information to left process\n");
                return 1;
            }
        }
        if (rank < last) {
            err = MPI_Bsend(&u[k - 1][max_m - 1], 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
            if (err != MPI_SUCCESS) {
                printf("Can't send information to right process\n");
                return 1;
            }
        }
    }

    double solving_time = MPI_Wtime() - start;

    if (!measurement) {
        MPI_File file;
        err = MPI_File_open(MPI_COMM_WORLD, OUTPUT_FILE, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &file);
        if (err != MPI_SUCCESS) {
            printf("Can't open file\n");
            return 1;
        }
        err = MPI_File_set_size(file, 0);
        if (err != MPI_SUCCESS) {
            printf("Can't clear file\n");
            return 1;
        }

        if (rank == 0) {
            err = MPI_File_write(file, &K, 1, MPI_INT, MPI_STATUS_IGNORE);
            if (err != MPI_SUCCESS) {
                printf("Can't write to the output file\n");
                return 1;
            }
            err = MPI_File_write(file, &M, 1, MPI_INT, MPI_STATUS_IGNORE);
            if (err != MPI_SUCCESS) {
                printf("Can't write to the output file\n");
                return 1;
            }
            err = MPI_File_write(file, &tau, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);
            if (err != MPI_SUCCESS) {
                printf("Can't write to the output file\n");
                return 1;
            }
            err = MPI_File_write(file, &h, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);
            if (err != MPI_SUCCESS) {
                printf("Can't write to the output file\n");
                return 1;
            }
            err = MPI_File_write(file, &c, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);
            if (err != MPI_SUCCESS) {
                printf("Can't write to the output file\n");
                return 1;
            }
            min_m = 0;
        }

        int header_size = 2 * sizeof(int) + 3 * sizeof(double);
        for (int k = 0; k <= K; k++) {
            int offset = header_size + (k * (M + 1) + min_m) * sizeof(double);
            err = MPI_File_write_at_all(file, offset, &u[k][min_m], max_m - min_m, MPI_DOUBLE, MPI_STATUS_IGNORE);
            if (err != MPI_SUCCESS) {
                printf("Can't save solution to the output file\n");
                return 1;
            }
        }
        MPI_File_close(&file);
    }

    double total_time = MPI_Wtime() - start;
    if (rank == 0) {
        printf("K = %d\n", K);
        printf("M = %d\n", M);
        printf("c = %.2f\n", c);
        printf("np = %d\n", world_size);
        printf("solving time: %.6f s\n", solving_time);
        printf("writing time: %.6f s\n", total_time - solving_time);
        printf("total   time: %.6f s\n", total_time);
    }

    // Free memory
    for (int i = 0; i <= K; i++) {
        free(u[i]);
    }
    free(u);
    MPI_Buffer_detach(&buffer, &buffer_size);
    free(buffer);

    err = MPI_Finalize();
    if (err != MPI_SUCCESS) {
        printf("Bad MPI finalization\n");
        return 1;
    }

    return 0;
}
