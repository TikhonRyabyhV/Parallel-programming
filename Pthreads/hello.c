#include <stdio.h>
#include <pthread.h>
#include <stdlib.h>

typedef struct {
    int thread_id;
    int total_threads;
} thread_data_t;

void* print_hello(void* arg) {
    thread_data_t* data = (thread_data_t*)arg;
    printf("Hello World from thread %d of %d\n", data->thread_id, data->total_threads);
    return NULL;
}

int main(int argc, char* argv[]) {
    int total_threads = 4;
    if (argc > 1) {
        total_threads = atoi(argv[1]);
    }

    if (total_threads <= 0) {
        fprintf(stderr, "Number of threads must be positive\n");
        return 1;
    }

    pthread_t* threads = malloc(total_threads * sizeof(pthread_t));
    thread_data_t* thread_data = malloc(total_threads * sizeof(thread_data_t));

    for (int i = 0; i < total_threads; i++) {
        thread_data[i].thread_id = i + 1;
        thread_data[i].total_threads = total_threads;
        if (pthread_create(&threads[i], NULL, print_hello, &thread_data[i]) != 0) {
            fprintf(stderr, "Error creating thread %d\n", i + 1);
            free(threads);
            free(thread_data);
            return 1;
        }
    }

    for (int i = 0; i < total_threads; i++) {
        if (pthread_join(threads[i], NULL) != 0) {
            fprintf(stderr, "Error joining thread %d\n", i + 1);
        }
    }

    free(threads);
    free(thread_data);

    return 0;
}
