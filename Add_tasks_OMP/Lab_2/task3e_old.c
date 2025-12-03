//
// D = (0, 0) => full-parallel
//
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#define ISIZE 5000
#define JSIZE 5000

double a[ISIZE][JSIZE];
double b[ISIZE][JSIZE];

int main(int argc, char **argv)
{
    int i, j;
    FILE *ff;

    for (i=0; i<ISIZE; i++){
        for (j=0; j<JSIZE; j++){
            a[i][j] = 10*i +j;
            b[i][j] = 0;
        }
    }

    double t_start = omp_get_wtime();

    #pragma omp parallel for private(j) collapse(2)
    for (i = 0; i < ISIZE; i++) {
        for (j = 0; j < JSIZE; j++) {
            a[i][j] = sin(0.002 * a[i][j]);
        }
    }

    #pragma omp parallel for private(j) collapse(2)
    for (i = 0; i < ISIZE - 4; i++) {
        for (j = 1; j < JSIZE; j++) {
            b[i][j] = a[i+4][j-1] * 1.5;
        }
    }

    double t_end = omp_get_wtime();
    
    printf("Time: %f\n", t_end - t_start);

    ff = fopen("result.txt", "w");
    for(i= 0; i < ISIZE; i++){
        for (j= 0; j < JSIZE; j++){
            fprintf(ff, "%f ", b[i][j]);
        }
        fprintf(ff, "\n");
    }
    fclose(ff);
    
    return 0;
}

