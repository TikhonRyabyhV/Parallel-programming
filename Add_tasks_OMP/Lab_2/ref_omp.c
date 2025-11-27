#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#define ISIZE 5000
#define JSIZE 5000

double a[ISIZE][JSIZE];

int main(int argc, char **argv)
{
    int i, j;
    FILE *ff;

    for (i=0; i<ISIZE; i++){
        for (j=0; j<JSIZE; j++){
            a[i][j] = 10*i +j;
        }
    }

    double t_start = omp_get_wtime();

    #pragma omp parallel for private(j) collapse(2)
    for (i=0; i<ISIZE; i++){
        for (j = 0; j < JSIZE; j++){
            a[i][j] = sin(2*a[i][j]);
        }
    }

    double t_end = omp_get_wtime();
    
    printf("Reference OpenMP Time: %f\n", t_end - t_start);

    ff = fopen("result_ref_omp.txt","w");
    for(i=0; i < ISIZE; i++){
        for (j=0; j < JSIZE; j++){
            fprintf(ff,"%f ",a[i][j]);
        }
        fprintf(ff,"\n");
    }
    fclose(ff);
    return 0;
}

