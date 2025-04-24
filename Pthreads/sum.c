#include <stdio.h>
#include <pthread.h>
#include <stdlib.h>

typedef struct {
    int start;
    int end;  
    double sum;
} thread_data_t;

void* compute_partial_sum(void* arg) {
    thread_data_t* data = (thread_data_t*)arg;
    data->sum = 0.0;
    for (int i = data->start; i <= data->end; i++) {
        data->sum += 1.0 / i;
    }
    return NULL;
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <N> <num_threads>\n", argv[0]);
        return 1;
    }

    int N = atoi(argv[1]);
    int num_threads = atoi(argv[2]);

    if (N <= 0 || num_threads <= 0) {
        fprintf(stderr, "N and num_threads must be positive\n");
        return 1;
    }

    if (num_threads > N) {
        num_threads = N;
    }

    pthread_t* threads = malloc(num_threads * sizeof(pthread_t));
    thread_data_t* thread_data = malloc(num_threads * sizeof(thread_data_t));

    int base_size = N / num_threads;
    int remainder = N % num_threads;
    int current_start = 1;

    for (int i = 0; i < num_threads; i++) {
        int extra = (i < remainder) ? 1 : 0;
        int block_size = base_size + extra;
        thread_data[i].start = current_start;
        thread_data[i].end = current_start + block_size - 1;
        thread_data[i].sum = 0.0;

        if (pthread_create(&threads[i], NULL, compute_partial_sum, &thread_data[i]) != 0) {
            fprintf(stderr, "Error creating thread %d\n", i + 1);
            free(threads);
            free(thread_data);
            return 1;
        }

        current_start += block_size;
    }

    for (int i = 0; i < num_threads; i++) {
        if (pthread_join(threads[i], NULL) != 0) {
            fprintf(stderr, "Error joining thread %d\n", i + 1);
        }
    }

    double total_sum = 0.0;
    for (int i = 0; i < num_threads; i++) {
        total_sum += thread_data[i].sum;
    }

    printf("Partial sum of harmonic series for N=%d: %.10f\n", N, total_sum);

    free(threads);
    free(thread_data);

    return 0;
}
