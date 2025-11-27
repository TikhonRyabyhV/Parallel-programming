#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

int main(int argc, char* argv[]) {
    if (argc < 2) {
        printf("Usage: %s <N>\n", argv[0]);
        return 1;
    }

    long long N = atoll(argv[1]);
    double sum = 0.0;

    #pragma omp parallel for reduction(+:sum)
    for (long long i = 1; i <= N; i++) {
        sum += 1.0 / (double)i;
    }

    printf("Sum: %.16f\n", sum);
    return 0;
}

