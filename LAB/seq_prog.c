#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <string.h>

#define OUTPUT_FILE "./results/conveq.bin"

const double T = 100.0;    // 0 <= t <= T
const double X = 100.0;    // 0 <= x <= X

const double tau = 0.007;
const double h   = 0.05 ;
const double c   = 1.35 ;

const int K = ((int) (T / tau));
const int M = ((int) (X / h  ));

double phi(double x) {

    return sin(x / 5.0);

}

double psi(double t) {

    return 2.0 * sin(t);

}

double f(double t, double x) {

    int s = ((int) x) % 100;

    if (s < 40 || s > 60) {
        return 0.0;
    }

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

    int measurement = (argc > 1) && (strcmp(argv[1], "-m") == 0);

    // Allocate grid dynamically
    double **u = (double **)malloc((K + 1) * sizeof(double *));
    for (int i = 0; i <= K; i++) {
        u[i] = (double *)malloc((M + 1) * sizeof(double));
    }

    double start = MPI_Wtime();

    // Initial conditions
    for (int m = 0; m <= M; m++) {
        u[0][m] = phi(m * h);
    }

    // Boundary conditions
    for (int k = 1; k <= K; k++) {
        u[k][0] = psi(k * tau);
    }

    // First layer by corner scheme
    for (int m = 1; m <= M; m++) {
        u[1][m] = solveByCornerScheme(u, 1, m);
    }

    // Remaining layers by cross scheme
    for (int k = 2; k <= K; k++) {
        for (int m = 1; m < M; m++) {
            u[k][m] = solveByCrossScheme(u, k, m);
        }
        u[k][M] = solveByCornerScheme(u, k, M);
    }

    double solving_time = MPI_Wtime() - start;

    // Write to binary file if not in measurement mode
    if (!measurement) {
        FILE *file = fopen(OUTPUT_FILE, "wb");
        if (!file) {
            printf("Error opening file %s\n", OUTPUT_FILE);
            for (int i = 0; i <= K; i++) {
                free(u[i]);
            }
            free(u);
            MPI_Finalize();
            return 1;
        }

        fwrite(&K, sizeof(int), 1, file);
        fwrite(&M, sizeof(int), 1, file);
        fwrite(&tau, sizeof(double), 1, file);
        fwrite(&h, sizeof(double), 1, file);
        fwrite(&c, sizeof(double), 1, file);
        for (int k = 0; k <= K; k++) {
            for (int m = 0; m <= M; m++) {
                fwrite(&u[k][m], sizeof(double), 1, file);
            }
        }
        fclose(file);
    }

    double total_time = MPI_Wtime() - start;

    // Output results
    printf("K = %d\n", K);
    printf("M = %d\n", M);
    printf("c = %.2f\n", c);
    printf("sequential\n");
    printf("solving time: %.6f s\n", solving_time);
    printf("writing time: %.6f s\n", total_time - solving_time);
    printf("total   time: %.6f s\n", total_time);

    // Free memory
    for (int i = 0; i <= K; i++) {
        free(u[i]);
    }
    free(u);

    err = MPI_Finalize();
    if (err != MPI_SUCCESS) {
        printf("Bad MPI finalization\n");
        return 1;
    }

    return 0;
}
