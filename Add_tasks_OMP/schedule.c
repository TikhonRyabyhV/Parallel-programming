#include <stdio.h>
#include <omp.h>

#define N 20  // number of iterations
#define NUM_THREADS 4  // number of threads

int iteration_to_thread[N];

void init_array() {
    for (int i = 0; i < N; i++) {
        iteration_to_thread[i] = -1;
    }
}


void print_results(const char* schedule, int chunk) {
    printf("\nSchedule: %s, chunk: %d\n", schedule, chunk);
    printf("Iteration -> Thread:\n");
    for (int i = 0; i < N; i++) {
        printf("Iteration %2d: Thread %d\n", i, iteration_to_thread[i]);
    }
}

void run_loop_default() {
    init_array();
    #pragma omp parallel for num_threads(NUM_THREADS)
    for (int i = 0; i < N; i++) {
        iteration_to_thread[i] = omp_get_thread_num();
    }
    print_results("default", 0);
}

void run_loop_static(int chunk) {
    init_array();
    #pragma omp parallel for num_threads(NUM_THREADS) schedule(static, chunk)
    for (int i = 0; i < N; i++) {
        iteration_to_thread[i] = omp_get_thread_num();
    }
    print_results("static", chunk);
}

void run_loop_dynamic(int chunk) {
    init_array();
    #pragma omp parallel for num_threads(NUM_THREADS) schedule(dynamic, chunk)
    for (int i = 0; i < N; i++) {
        iteration_to_thread[i] = omp_get_thread_num();
    }
    print_results("dynamic", chunk);
}

void run_loop_guided(int chunk) {
    init_array();
    #pragma omp parallel for num_threads(NUM_THREADS) schedule(guided, chunk)
    for (int i = 0; i < N; i++) {
        iteration_to_thread[i] = omp_get_thread_num();
    }
    print_results("guided", chunk);
}

int main() {
    // Set number of threads
    omp_set_num_threads(NUM_THREADS);

    printf("=== Tests for balance methods ===\n");

    // 1. Default (without schedule)
    run_loop_default();

    // 2. Static, chunk=1
    run_loop_static(1);

    // 3. Static, chunk=4
    run_loop_static(4);

    // 4. Dynamic, chunk=1
    run_loop_dynamic(1);

    // 5. Dynamic, chunk=4
    run_loop_dynamic(4);

    // 6. Guided, chunk=1
    run_loop_guided(1);

    // 7. Guided, chunk=4
    run_loop_guided(4);

    return 0;
}
