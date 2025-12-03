#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#define ISIZE 5000
#define JSIZE 5000

int main(int argc, char **argv) {
    double (*a)[JSIZE] = malloc(sizeof(double[ISIZE][JSIZE]));
    if (a == NULL) return 1;

    #pragma omp parallel for collapse(2)
    for (int i = 0; i < ISIZE; i++) {
        for (int j = 0; j < JSIZE; j++) {
            a[i][j] = 10 * i + j;
        }
    }

    // start timer
    double tstart = omp_get_wtime();

    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        int nthreads = omp_get_num_threads();

        // chunking external cycle
        int N = ISIZE - 4;
        int chunk = N / nthreads;
        int remainder = N % nthreads;

        int start, end;
        if (tid < remainder) {
            start = tid * (chunk + 1);
            end = start + chunk + 1;
        } else {
            start = tid * chunk + remainder;
            end = start + chunk;
        }

        // local copy (shadow) of 2 lines for each process
        double shadow[2][JSIZE];
        
        if (end < ISIZE) { 
             for (int j = 0; j < JSIZE; j++) {
                 shadow[0][j] = a[end][j];
             }
        }
        if (end + 1 < ISIZE) {
             for (int j = 0; j < JSIZE; j++) {
                 shadow[1][j] = a[end + 1][j];
             }
        }

        // wait for saving copies
        #pragma omp barrier

        for (int i = start; i < end; i++) {
            for (int j = 0; j < JSIZE - 4; j++) {
                double val;
                int row_idx = i + 2;

                if (row_idx < end) {
                    val = a[row_idx][j+4];
                } else {
                    int shadow_idx = row_idx - end;
                    val = shadow[shadow_idx][j+4];
                }

                a[i][j] = sin(0.1 * val);
            }
        }
    }
    
    // end timer
    double tend = omp_get_wtime();

    // exec time info
    printf("Task 2z finished. Time: %f\n", tend - tstart);

    FILE *ff = fopen("result_omp_opt.txt", "w");
    for (int i = 0; i < ISIZE; i++) {
        for (int j = 0; j < JSIZE; j++) {
            fprintf(ff, "%f ", a[i][j]);
        }
        fprintf(ff, "\n");
    }
    fclose(ff);
    free(a);

    return 0;
}

