#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include <time.h>

typedef double value_t;

typedef struct {
    void (*solver)(value_t*, value_t*, value_t*, value_t*, value_t*, int);
    int max_iterations;
    value_t epsilon;
    value_t mu;  // mu = a * h^2 / 12
} params_t;

void require(int condition, const char *message) {
    if (!condition) {
        fprintf(stderr, "Error: %s\n", message);
        exit(1);
    }
}

// Tridiagonal Matrix Algorithm
void solve_tridiagonal(value_t *delta, value_t *A, value_t *B, value_t *C, value_t *F, int N) {
    value_t *cp = (value_t*)malloc(N * sizeof(value_t));
    value_t *dp = (value_t*)malloc(N * sizeof(value_t));

    cp[0] = C[0] / B[0];
    dp[0] = -F[0] / B[0];
    
    for (int i = 1; i < N; i++) {
        value_t denom = B[i] - A[i] * cp[i-1];
        cp[i] = (i == N-1) ? 0.0 : (C[i] / denom);
        dp[i] = -(F[i] + A[i] * dp[i-1]) / denom;
    }

    delta[N-1] = dp[N-1];
    for (int i = N - 2; i >= 0; i--) {
        delta[i] = dp[i] - cp[i] * delta[i+1];
    }

    free(cp);
    free(dp);
}

// Recursive part of cyclic reduction
void solve_CR_impl(value_t *delta, value_t *a, value_t *b, value_t *c, value_t *f, int N, int shift) {
    int base = shift - 1;

    // base case
    if (N == 1) {
        delta[base] = -f[0] / b[0];
        return;
    }

    int m = N / 2;

    // Forward Reduction
    #pragma omp parallel for
    for (int i = 0; i < m; i++) {
        int idx = 2 * i + 1;
        
        value_t alpha = (idx - 1 >= 0) ? a[idx] / b[idx - 1] : 0.0;
        value_t gamma = (idx + 1 < N) ? c[idx] / b[idx + 1] : 0.0;

        a[i + N] = -alpha * a[idx - 1];
        b[i + N] = b[idx] - alpha * c[idx - 1] - gamma * ((idx + 1 < N) ? a[idx + 1] : 0.0);
        c[i + N] = -gamma * ((idx + 1 < N) ? c[idx + 1] : 0.0);
        f[i + N] = f[idx] - alpha * f[idx - 1] - gamma * ((idx + 1 < N) ? f[idx + 1] : 0.0);
    }

    solve_CR_impl(delta, a + N, b + N, c + N, f + N, m, 2 * shift);

    // Backward Substitution
    #pragma omp parallel for
    for (int i = 0; i < m; i++) {
        int k = 2 * i;
        value_t rhs = -f[k];
        
        if (i > 0) 
            rhs -= a[k] * delta[shift * (k - 1) + base];
        
        if (k + 1 < N) 
            rhs -= c[k] * delta[shift * (k + 1) + base];
            
        delta[shift * k + base] = rhs / b[k];
    }

    if (N % 2 == 1) {
        value_t rhs = -f[N-1] - a[N-1] * delta[shift * (N - 2) + base];
        delta[shift * (N-1) + base] = rhs / b[N-1];
    }
}

void solve_CR(value_t *delta, value_t *A, value_t *B, value_t *C, value_t *F, int N) {
    solve_CR_impl(delta, A, B, C, F, N, 1);
}

// f(y) = y^3 - y
value_t func_y(value_t y) {
    return y * y * y - y;
}

value_t calc_F(int m, const value_t *y, value_t mu) {
    return y[m+1] - 2*y[m] + y[m-1] - mu * (func_y(y[m+1]) + 10*func_y(y[m]) + func_y(y[m-1]));
}

// L2-norm 
value_t get_norm(const value_t *vec, int N) {
    value_t sum = 0.0;
    #pragma omp parallel for reduction(+:sum)
    for (int m = 0; m < N; m++) {
        sum += vec[m] * vec[m];
    }
    return sqrt(sum);
}

// Newton method
int solve_NODE(value_t *grid, int N_grid, params_t *params) {
    int N_sys = N_grid - 2;
    
    value_t *delta = (value_t*)calloc(N_sys, sizeof(value_t));
    value_t *A = (value_t*)calloc(2 * N_sys, sizeof(value_t));
    value_t *B = (value_t*)calloc(2 * N_sys, sizeof(value_t));
    value_t *C = (value_t*)calloc(2 * N_sys, sizeof(value_t));
    value_t *F_vals = (value_t*)calloc(2 * N_sys, sizeof(value_t));

    int iter = 0;
    for (; iter < params->max_iterations; iter++) {
        
        #pragma omp parallel for
        for (int m = 1; m < N_grid - 1; m++) {
            int i = m - 1;
            
            // f'(y) = 3y^2 - 1
            value_t df_minus = 3 * grid[m-1] * grid[m-1] - 1;
            value_t df_curr  = 3 * grid[m]   * grid[m]   - 1;
            value_t df_plus  = 3 * grid[m+1] * grid[m+1] - 1;

            A[i] = 1.0 - params->mu * df_minus;
            B[i] = -2.0 - 10.0 * params->mu * df_curr;
            C[i] = 1.0 - params->mu * df_plus;
            
            F_vals[i] = calc_F(m, grid, params->mu);
        }

        if (get_norm(F_vals, N_sys) < params->epsilon) {
            break;
        }

        // Solution of J * delta = -F
        params->solver(delta, A, B, C, F_vals, N_sys);

        #pragma omp parallel for
        for (int m = 1; m < N_grid - 1; m++) {
            grid[m] += delta[m - 1];
        }
    }

    free(delta); free(A); free(B); free(C); free(F_vals);
    return iter;
}

int main(int argc, char **argv) {
    params_t params;
    params.solver = solve_tridiagonal;
    params.epsilon = 1e-8;
    params.max_iterations = 1000;

    if (argc < 3) {
        fprintf(stderr, "Usage: %s [a] [h] {seq|CR[threads]}\n", argv[0]);
        return 1;
    }

    double a = strtod(argv[1], NULL);
    double h = strtod(argv[2], NULL);
    params.mu = a * h * h / 12.0;

    // Choose solver
    if (argc == 4 && strncmp(argv[3], "CR", 2) == 0) {
        int threads = atoi(argv[3] + 2);
        fprintf(stderr, "solver: cyclic reduction with %d threads\n", threads);
        omp_set_num_threads(threads);
        params.solver = solve_CR;
    } else {
        fprintf(stderr, "solver: sequential\n");
        omp_set_num_threads(1);
        params.solver = solve_tridiagonal;
    }

    int N = (int)(20.0 / h); 
    value_t *grid = (value_t*)malloc(N * sizeof(value_t));
    
    for (int i = 0; i < N; i++) {
        grid[i] = sqrt(2.0);
    }

    double start_time = omp_get_wtime();
    int iterations = solve_NODE(grid, N, &params);
    double end_time = omp_get_wtime();

    fprintf(stderr, "iterations: %d\n", iterations);
    fprintf(stderr, "time: %f s\n", end_time - start_time);

    // output info
    for (int i = 0; i < N; i++) {
        printf("%.10f ", grid[i]);
    }
    printf("\n");

    free(grid);
    return 0;
}

