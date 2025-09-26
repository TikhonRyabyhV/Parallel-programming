#include <stdio.h>
#include <omp.h>

static int value = 0;
#pragma omp threadprivate(value)

int main() {

    value = 100;

    printf("value (main thread) before parallel part: %d\n", value);

    // Parallel part with copyin
    #pragma omp parallel num_threads(4) copyin(value)
    {
        int thread_id = omp_get_thread_num();
        printf("Thread %d: initial value = %d\n", thread_id, value);

        // One thread changes value in single block
        #pragma omp single // copyprivate(value)
        {
            value = 200 + thread_id;
            printf("Thread %d in single block: changed value = %d\n", thread_id, value);
        }

        // All threads print value after copyprivate
        #pragma omp barrier
        printf("Thread %d: value after copyprivate = %d\n", thread_id, value);
    }

    // value in main thread after parallel part
    printf("value in main thread after parallel part: %d\n", value);

    return 0;
}
