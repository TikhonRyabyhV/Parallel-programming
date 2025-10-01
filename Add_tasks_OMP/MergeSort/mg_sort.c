#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>
#include <omp.h>
#include <sys/time.h>
#include <string.h>

// --- Sort function ---
void merge_into_temp(int arr[], int left, int mid, int right, int temp_arr[]) {
    int i = left;
    int j = mid + 1;
    int k = left;
    while (i <= mid && j <= right) {
        if (arr[i] <= arr[j]) {
            temp_arr[k++] = arr[i++];
        } else {
            temp_arr[k++] = arr[j++];
        }
    }
    while (i <= mid) {
        temp_arr[k++] = arr[i++];
    }
    while (j <= right) {
        temp_arr[k++] = arr[j++];
    }
}

void sequentialMergeSort(int arr[], int left, int right, int temp_arr[]) {
    if (left < right) {
        int mid = left + (right - left) / 2;
        sequentialMergeSort(arr, left, mid, temp_arr);
        sequentialMergeSort(arr, mid + 1, right, temp_arr);
        merge_into_temp(arr, left, mid, right, temp_arr);
        memcpy(arr + left, temp_arr + left, (right - left + 1) * sizeof(int));
    }
}

void parallelMergeSort(int arr[], int left, int right, int sequential_threshold, int temp_arr[]) {
    if (left < right) {
        if (right - left + 1 <= sequential_threshold) {
            sequentialMergeSort(arr, left, right, temp_arr);
        } else {
            int mid = left + (right - left) / 2;
            #pragma omp task shared(arr, temp_arr)
            parallelMergeSort(arr, left, mid, sequential_threshold, temp_arr);
            #pragma omp task shared(arr, temp_arr)
            parallelMergeSort(arr, mid + 1, right, sequential_threshold, temp_arr);
            #pragma omp taskwait
            merge_into_temp(arr, left, mid, right, temp_arr);
            memcpy(arr + left, temp_arr + left, (right - left + 1) * sizeof(int));
        }
    }
}

bool isSorted(const int arr[], int size) {
    for (int i = 0; i < size - 1; ++i) {
        if (arr[i] > arr[i+1]) {
            return false;
        }
    }
    return true;
}

double get_wall_time() {
    struct timeval time;
    if (gettimeofday(&time, NULL)) {
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}

int main() {
    srand(time(NULL)); // Random numbers generator initialization

    // Arrays sizes
    int N_sizes[] = {1 << 16, 1 << 18, 1 << 20, 1 << 22}; // 65536, 262144, 1048576, 4194304
    int num_N_sizes = sizeof(N_sizes) / sizeof(N_sizes[0]);

    // Threshold values
    int thresholds[] = {16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384}; // Уменьшим для простоты, если нужно больше, добавьте
    int num_thresholds = sizeof(thresholds) / sizeof(thresholds[0]);
    
    // Thread numbers
    int num_threads_arr[] = {1, 2, 3, 4, 5, 6, 7, 8};
    int num_num_threads = sizeof(num_threads_arr) / sizeof(num_threads_arr[0]);

    const int MAX_N = N_sizes[num_N_sizes - 1]; 
    int *arr = (int *)malloc(MAX_N * sizeof(int));
    int *arr_initial_copy = (int *)malloc(MAX_N * sizeof(int));
    int *global_temp_arr = (int *)malloc(MAX_N * sizeof(int));

    if (arr == NULL || arr_initial_copy == NULL || global_temp_arr == NULL) {
        perror("Memory allocation error");
        return 1;
    }

    // --- Main cycle ---
    // Output format: N_SIZE THREADS THRESHOLD SEQ_TIME PAR_TIME SPEEDUP EFFICIENCY IS_CORRECT
    // IS_CORRECT: 1 - correct, 0 - incorrect
    printf("N_SIZE THREADS THRESHOLD SEQ_TIME PAR_TIME SPEEDUP EFFICIENCY IS_CORRECT\n");

    for (int ni = 0; ni < num_N_sizes; ++ni) {
        int N = N_sizes[ni];
        
        // Fill arr_initial_copy with random data for current N
        for (int i = 0; i < N; ++i) {
            arr_initial_copy[i] = rand() % N;
        }

        // --- Measuring time for current N (sequential version) ---
        memcpy(arr, arr_initial_copy, N * sizeof(int));
        double start_seq = get_wall_time();
        sequentialMergeSort(arr, 0, N - 1, global_temp_arr);
        double end_seq = get_wall_time();
        double duration_seq = end_seq - start_seq;
        
        bool seq_correct = isSorted(arr, N);
        if (!seq_correct) {
            fprintf(stderr, "ERROR: Sequential sort failed for N=%d!\n", N);
        }
        printf("%d 0 0 %f %f %f %f %d\n", N, duration_seq, duration_seq, 1.0, 1.0, (int)seq_correct);


        // --- Parallel version for current N ---
        for (int t_idx = 0; t_idx < num_num_threads; ++t_idx) {
            int num_threads = num_threads_arr[t_idx];
            omp_set_num_threads(num_threads);

            for (int th_idx = 0; th_idx < num_thresholds; ++th_idx) {
                int threshold = thresholds[th_idx];
                
                memcpy(arr, arr_initial_copy, N * sizeof(int)); // array reset
                double start_par = get_wall_time();
                #pragma omp parallel
                #pragma omp single
                parallelMergeSort(arr, 0, N - 1, threshold, global_temp_arr);
                double end_par = get_wall_time();
                double duration_par = end_par - start_par;

                bool par_correct = isSorted(arr, N);
                if (!par_correct) {
                    fprintf(stderr, "ERROR: Parallel sort failed for N=%d, T=%d, P=%d!\n", N, threshold, num_threads);
                }

                double speedup = duration_seq / duration_par;
                double efficiency = speedup / num_threads;

                printf("%d %d %d %f %f %f %f %d\n", N, num_threads, threshold, duration_seq, duration_par, speedup, efficiency, (int)par_correct);
            }
        }
    }

    free(arr);
    free(arr_initial_copy);
    free(global_temp_arr);
    
    return 0;
}
