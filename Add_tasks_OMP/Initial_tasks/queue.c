#include <stdio.h>
#include <omp.h>

int main() {
    int shared_var = 0;      
    volatile int turn = 0;

    #pragma omp parallel shared(shared_var, turn)
    {
        int rank = omp_get_thread_num();

        while (turn != rank) {
            #pragma omp flush(turn)
        }

        
        shared_var++;
        printf("Thread %d accessed shared_var. New value: %d\n", rank, shared_var);


        turn++;
        #pragma omp flush(turn)
    }

    return 0;
}

