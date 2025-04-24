#include <stdio.h>
#include <pthread.h>
#include <stdlib.h>
#include <unistd.h>

int shared_value = 0;

pthread_mutex_t mutex;

pthread_cond_t* cond_vars;

int current_thread = 0;

typedef struct {
    int thread_id;
    int total_threads;
} thread_data_t;

void* thread_function(void* arg) {
    thread_data_t* data = (thread_data_t*)arg;
    int my_id = data->thread_id - 1;

    pthread_mutex_lock(&mutex);

    while (current_thread != my_id) {
        pthread_cond_wait(&cond_vars[my_id], &mutex);
    }

    shared_value++;
    printf("Thread %d incremented shared_value to %d\n", data->thread_id, shared_value);

    current_thread = (current_thread + 1) % data->total_threads;
    if (current_thread == 0) {
        pthread_cond_signal(&cond_vars[0]);
    } else {
        pthread_cond_signal(&cond_vars[current_thread]);
    }

    pthread_mutex_unlock(&mutex);
    return NULL;
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        fprintf(stderr, "Usage: %s <num_threads>\n", argv[0]);
        return 1;
    }

    int num_threads = atoi(argv[1]);
    if (num_threads <= 0) {
        fprintf(stderr, "Number of threads must be positive\n");
        return 1;
    }

    if (pthread_mutex_init(&mutex, NULL) != 0) {
        fprintf(stderr, "Mutex initialization failed\n");
        return 1;
    }

    cond_vars = malloc(num_threads * sizeof(pthread_cond_t));
    for (int i = 0; i < num_threads; i++) {
        if (pthread_cond_init(&cond_vars[i], NULL) != 0) {
            fprintf(stderr, "Condition variable initialization failed\n");
            free(cond_vars);
            pthread_mutex_destroy(&mutex);
            return 1;
        }
    }

    pthread_t* threads = malloc(num_threads * sizeof(pthread_t));
    thread_data_t* thread_data = malloc(num_threads * sizeof(thread_data_t));

    for (int i = 0; i < num_threads; i++) {
        thread_data[i].thread_id = i + 1;
        thread_data[i].total_threads = num_threads;
        if (pthread_create(&threads[i], NULL, thread_function, &thread_data[i]) != 0) {
            fprintf(stderr, "Error creating thread %d\n", i + 1);
            free(threads);
            free(thread_data);
            free(cond_vars);
            pthread_mutex_destroy(&mutex);
            return 1;
        }
    }

    for (int i = 0; i < num_threads; i++) {
        if (pthread_join(threads[i], NULL) != 0) {
            fprintf(stderr, "Error joining thread %d\n", i + 1);
        }
    }

    for (int i = 0; i < num_threads; i++) {
        pthread_cond_destroy(&cond_vars[i]);
    }
    pthread_mutex_destroy(&mutex);
    free(threads);
    free(thread_data);
    free(cond_vars);

    printf("Final shared_value: %d\n", shared_value);
    return 0;
}
