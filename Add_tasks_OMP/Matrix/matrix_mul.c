#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <immintrin.h>
#include <math.h>
#include <time.h>

#define ALIGN 32
#define BLOCK_SIZE 64
#define STRASSEN_THRESHOLD 128 

float get_rand_float() { return (float)rand() / RAND_MAX; }
void fill_matrix(float* A, int N) { for (int i=0; i<N*N; i++) A[i] = get_rand_float(); }
void zero_matrix(float* A, int N) { memset(A, 0, N*N*sizeof(float)); }

void transpose(const float* Src, float* Dst, int N) {
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < N; i += BLOCK_SIZE) {
        for (int j = 0; j < N; j += BLOCK_SIZE) {
            for (int ii = i; ii < i + BLOCK_SIZE && ii < N; ii++) {
                for (int jj = j; jj < j + BLOCK_SIZE && jj < N; jj++) {
                    Dst[jj * N + ii] = Src[ii * N + jj];
                }
            }
        }
    }
}

void matmul_naive(const float* A, const float* B, float* C, int N) {
    zero_matrix(C, N);
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            float sum = 0.0f;
            for (int k = 0; k < N; k++) sum += A[i * N + k] * B[k * N + j];
            C[i * N + j] = sum;
        }
    }
}

void matmul_transposed(const float* A, const float* B, float* C, int N) {
    float* B_T = (float*)aligned_alloc(ALIGN, N*N*sizeof(float));
    transpose(B, B_T, N);
    zero_matrix(C, N);
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            float sum = 0.0f;
            for (int k = 0; k < N; k++) sum += A[i * N + k] * B_T[j * N + k];
            C[i * N + j] = sum;
        }
    }
    free(B_T);
}

void matmul_blocked(const float* A, const float* B, float* C, int N) {
    zero_matrix(C, N);
    int bs = 64; 
    #pragma omp parallel for collapse(2)
    for (int ii = 0; ii < N; ii += bs) {
        for (int jj = 0; jj < N; jj += bs) {
            for (int kk = 0; kk < N; kk += bs) {
                for (int i = ii; i < ((ii + bs) > N ? N : (ii + bs)); i++) {
                    for (int j = jj; j < ((jj + bs) > N ? N : (jj + bs)); j++) {
                        float sum = C[i * N + j];
                        for (int k = kk; k < ((kk + bs) > N ? N : (kk + bs)); k++) {
                            sum += A[i * N + k] * B[k * N + j];
                        }
                        C[i * N + j] = sum;
                    }
                }
            }
        }
    }
}

void add_mat(const float* A, const float* B, float* C, int N) {
    #pragma omp parallel for
    for (int i = 0; i < N * N; i++) C[i] = A[i] + B[i];
}
void sub_mat(const float* A, const float* B, float* C, int N) {
    #pragma omp parallel for
    for (int i = 0; i < N * N; i++) C[i] = A[i] - B[i];
}
void matmul_naive_kernel(const float* A, const float* B, float* C, int N) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            float sum = 0.0f;
            for (int k = 0; k < N; k++) sum += A[i * N + k] * B[k * N + j];
            C[i * N + j] = sum;
        }
    }
}
void strassen_recursive(const float* A, const float* B, float* C, int N) {
    if (N <= STRASSEN_THRESHOLD) { matmul_naive_kernel(A, B, C, N); return; }
    int K = N / 2;
    int size = K * K * sizeof(float);
    float *A11 = (float*)aligned_alloc(ALIGN, size), *A12 = (float*)aligned_alloc(ALIGN, size);
    float *A21 = (float*)aligned_alloc(ALIGN, size), *A22 = (float*)aligned_alloc(ALIGN, size);
    float *B11 = (float*)aligned_alloc(ALIGN, size), *B12 = (float*)aligned_alloc(ALIGN, size);
    float *B21 = (float*)aligned_alloc(ALIGN, size), *B22 = (float*)aligned_alloc(ALIGN, size);
    float *C11 = (float*)aligned_alloc(ALIGN, size), *C12 = (float*)aligned_alloc(ALIGN, size);
    float *C21 = (float*)aligned_alloc(ALIGN, size), *C22 = (float*)aligned_alloc(ALIGN, size);
    float *P1 = (float*)aligned_alloc(ALIGN, size), *P2 = (float*)aligned_alloc(ALIGN, size);
    float *P3 = (float*)aligned_alloc(ALIGN, size), *P4 = (float*)aligned_alloc(ALIGN, size);
    float *P5 = (float*)aligned_alloc(ALIGN, size), *P6 = (float*)aligned_alloc(ALIGN, size);
    float *P7 = (float*)aligned_alloc(ALIGN, size);
    float *T1 = (float*)aligned_alloc(ALIGN, size), *T2 = (float*)aligned_alloc(ALIGN, size);

    #pragma omp parallel for collapse(2)
    for (int i = 0; i < K; i++) {
        for (int j = 0; j < K; j++) {
            A11[i*K+j] = A[i*N+j];         A12[i*K+j] = A[i*N+(j+K)];
            A21[i*K+j] = A[(i+K)*N+j];     A22[i*K+j] = A[(i+K)*N+(j+K)];
            B11[i*K+j] = B[i*N+j];         B12[i*K+j] = B[i*N+(j+K)];
            B21[i*K+j] = B[(i+K)*N+j];     B22[i*K+j] = B[(i+K)*N+(j+K)];
        }
    }

    add_mat(A11, A22, T1, K); add_mat(B11, B22, T2, K); strassen_recursive(T1, T2, P1, K);
    add_mat(A21, A22, T1, K); strassen_recursive(T1, B11, P2, K);
    sub_mat(B12, B22, T2, K); strassen_recursive(A11, T2, P3, K);
    sub_mat(B21, B11, T2, K); strassen_recursive(A22, T2, P4, K);
    add_mat(A11, A12, T1, K); strassen_recursive(T1, B22, P5, K);
    sub_mat(A21, A11, T1, K); add_mat(B11, B12, T2, K); strassen_recursive(T1, T2, P6, K);
    sub_mat(A12, A22, T1, K); add_mat(B21, B22, T2, K); strassen_recursive(T1, T2, P7, K);

    add_mat(P1, P4, T1, K); sub_mat(T1, P5, T2, K); add_mat(T2, P7, C11, K);
    add_mat(P3, P5, C12, K); add_mat(P2, P4, C21, K);
    sub_mat(P1, P2, T1, K); add_mat(T1, P3, T2, K); add_mat(T2, P6, C22, K);

    #pragma omp parallel for collapse(2)
    for (int i = 0; i < K; i++) {
        for (int j = 0; j < K; j++) {
            C[i*N+j] = C11[i*K+j]; C[i*N+(j+K)] = C12[i*K+j];
            C[(i+K)*N+j] = C21[i*K+j]; C[(i+K)*N+(j+K)] = C22[i*K+j];
        }
    }
    free(A11); free(A12); free(A21); free(A22); free(B11); free(B12); free(B21); free(B22);
    free(C11); free(C12); free(C21); free(C22); free(P1); free(P2); free(P3); free(P4); free(P5); free(P6); free(P7); free(T1); free(T2);
}

void matmul_simd(const float* A, const float* B, float* C, int N) {
    float* B_T = (float*)aligned_alloc(ALIGN, N*N*sizeof(float));
    transpose(B, B_T, N);
    zero_matrix(C, N);
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            __m256 sum_vec = _mm256_setzero_ps();
            int k = 0;
            for (; k <= N - 8; k += 8) {
                __m256 a = _mm256_loadu_ps(&A[i * N + k]);
                __m256 b = _mm256_loadu_ps(&B_T[j * N + k]);
                sum_vec = _mm256_fmadd_ps(a, b, sum_vec);
            }
            float temp[8];
            _mm256_storeu_ps(temp, sum_vec);
            float sum = temp[0] + temp[1] + temp[2] + temp[3] + temp[4] + temp[5] + temp[6] + temp[7];
            for (; k < N; k++) sum += A[i * N + k] * B_T[j * N + k];
            C[i * N + j] = sum;
        }
    }
    free(B_T);
}

void matmul_fast_combined(const float* A, const float* B, float* C, int N) {
    float* B_T = (float*)aligned_alloc(ALIGN, N*N*sizeof(float));
    transpose(B, B_T, N);
    #pragma omp parallel for collapse(2) schedule(dynamic)
    for (int i = 0; i < N; i += BLOCK_SIZE) {
        for (int j = 0; j < N; j += BLOCK_SIZE) {
            for (int ii = i; ii < i + BLOCK_SIZE && ii < N; ii++) {
                for (int jj = j; jj < j + BLOCK_SIZE && jj < N; jj++) {
                    __m256 sum_vec = _mm256_setzero_ps();
                    int k = 0;
                    for (; k <= N - 8; k += 8) {
                        __m256 a = _mm256_loadu_ps(&A[ii * N + k]);
                        __m256 b = _mm256_loadu_ps(&B_T[jj * N + k]);
                        sum_vec = _mm256_fmadd_ps(a, b, sum_vec);
                    }
                    float temp[8];
                    _mm256_storeu_ps(temp, sum_vec);
                    float sum = temp[0] + temp[1] + temp[2] + temp[3] + temp[4] + temp[5] + temp[6] + temp[7];
                    for (; k < N; k++) sum += A[ii * N + k] * B_T[jj * N + k];
                    C[ii * N + jj] = sum;
                }
            }
        }
    }
    free(B_T);
}

int main() {
    int sizes[] = {64, 128, 256, 512, 1024, 2048};
    int num_sizes = sizeof(sizes) / sizeof(sizes[0]);

    printf("Size\tNaive_1\tTransp_1\tBlock_1\tStrass_1\tSIMD_1\tFast_1\tNaive_8\tTransp_8\tBlock_8\tStrass_8\tSIMD_8\tFast_8\n");

    for (int s = 0; s < num_sizes; s++) {
        int N = sizes[s];
        size_t bytes = N * N * sizeof(float);
        float *A = (float*)aligned_alloc(ALIGN, bytes);
        float *B = (float*)aligned_alloc(ALIGN, bytes);
        float *C = (float*)aligned_alloc(ALIGN, bytes);
        fill_matrix(A, N); fill_matrix(B, N);

        printf("%d\t", N); fflush(stdout);

        for (int t = 0; t < 2; t++) {
            int num_threads = (t == 0) ? 1 : 8;
            omp_set_num_threads(num_threads);

            double start, end;

            if (N <= 1024) {
                start = omp_get_wtime(); matmul_naive(A, B, C, N); end = omp_get_wtime();
                printf("%.4f\t", end - start);
            } else { printf("Skip\t"); } fflush(stdout);

            start = omp_get_wtime(); matmul_transposed(A, B, C, N); end = omp_get_wtime();
            printf("%.4f\t", end - start); fflush(stdout);

            start = omp_get_wtime(); matmul_blocked(A, B, C, N); end = omp_get_wtime();
            printf("%.4f\t", end - start); fflush(stdout);

            if (N <= 2048) {
                start = omp_get_wtime(); strassen_recursive(A, B, C, N); end = omp_get_wtime();
                printf("%.4f\t", end - start);
            } else { printf("Skip\t"); } fflush(stdout);

            start = omp_get_wtime(); matmul_simd(A, B, C, N); end = omp_get_wtime();
            printf("%.4f\t", end - start); fflush(stdout);

            start = omp_get_wtime(); matmul_fast_combined(A, B, C, N); end = omp_get_wtime();
            if (t == 1) printf("%.4f\n", end - start);
            else printf("%.4f\t", end - start);
            fflush(stdout);
        }
        free(A); free(B); free(C);
    }
    return 0;
}

