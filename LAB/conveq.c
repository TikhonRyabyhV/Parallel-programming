#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <mpi.h>

// Параметры задачи
#define T 100.0    // Конечное время
#define X 100.0    // Конечное пространство
#define a 1.0      // Коэффициент a в уравнении
#define Nx 4000    // Число узлов по x
#define Nt 5000    // Число узлов по t
#define SYNC_INTERVAL 100 // Синхронизация каждые 100 шагов

// Функции для начальных и граничных условий, а также правой части
double phi(double x) {
    return sin(M_PI * x / X); // Начальное условие u(0, x) = sin(pi * x / X)
}

double psi(double t) {
    return exp(-t / T); // Граничное условие u(t, 0) = exp(-t / T)
}

double f(double x, double t) {
    return 0.0; // Правая часть f(x, t) = 0
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Проверка флага -nf
    int write_to_file = 1;
    if (argc > 1 && strcmp(argv[1], "-nf") == 0) {
        write_to_file = 0;
        if (rank == 0) {
            printf("Запись в файл отключена (флаг -nf)\n");
        }
    }

    // Шаги сетки
    double tau = T / Nt;
    double h = X / Nx;

    // Проверка условия устойчивости
    if (rank == 0) {
        double courant = tau / h;
        printf("Число Куранта: tau/h = %f (должно быть <= 1 для устойчивости схемы 'уголок')\n", courant);
        if (courant > 1.0 / a) {
            printf("Внимание: схема 'уголок' может быть неустойчивой (tau/h = %f, должно быть <= %f)\n", courant, 1.0 / a);
        }
    }

    // Статическое разбиение сетки
    int local_nx = Nx / size;
    int remainder = Nx % size;
    int local_start, local_size, local_end;

    if (rank < remainder) {
        local_size = local_nx + 1;
        local_start = rank * local_size;
    } else {
        local_size = local_nx;
        local_start = rank * local_nx + remainder;
    }
    local_end = local_start + local_size;
    if (local_start == 0) local_start = 1;
    local_size = local_end - local_start;

    // Массив для u
    double **u = (double **)malloc((Nt + 1) * sizeof(double *));
    for (int k = 0; k <= Nt; k++) {
        u[k] = (double *)calloc(Nx + 1, sizeof(double));
    }

    // Инициализация начального условия u(0, x) = phi(x)
    for (int m = 0; m <= Nx; m++) {
        u[0][m] = phi(m * h);
    }

    // Инициализация граничного условия u(t, 0) = psi(t)
    if (rank == 0) {
        for (int k = 1; k <= Nt; k++) {
            u[k][0] = psi(k * tau);
        }
    }

    // Синхронизация начального условия u[0][m] и u[k][0]
    for (int m = 0; m <= Nx; m++) {
        MPI_Bcast(&u[0][m], 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
    for (int k = 1; k <= Nt; k++) {
        MPI_Bcast(&u[k][0], 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    // Переменные для динамической балансировки
    int *next_m = (int *)calloc(size, sizeof(int));
    for (int i = 0; i < size; i++) {
        if (i < remainder) {
            next_m[i] = i * (local_nx + 1);
        } else {
            next_m[i] = i * local_nx + remainder;
        }
        if (next_m[i] == 0) next_m[i] = 1;
    }

    // Переменные для замера времени
    double start_time, end_time, total_time;
    double compute_start, compute_end, local_compute_time = 0.0, compute_time;

    // Синхронизация перед началом замера
    MPI_Barrier(MPI_COMM_WORLD);
    start_time = MPI_Wtime();

    // Основной цикл по времени
    for (int k = 0; k < Nt; k++) {
        // Обмен u[k][m-1] для схемы "уголок"
        double u_prev_k;
        if (rank == 0) {
            u_prev_k = u[k][0];
        } else {
            MPI_Recv(&u_prev_k, 1, MPI_DOUBLE, rank - 1, k, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        // Вычисления для локального поддомена схемой прямоугольник
        compute_start = MPI_Wtime();

        double u_prev_k1 = rank == 0 ? u[k][0] : u_prev_k;
        int m = local_start;
        while (m <= Nx) {
            if (m >= local_end) {
                int found = 0;
                for (int p = 0; p < size; p++) {
                    if (next_m[p] <= Nx) {
                        m = next_m[p];
                        next_m[p] += local_nx;
                        if (next_m[p] > Nx) next_m[p] = Nx + 1;
                        found = 1;
                        break;
                    }
                }
                if (!found) break;
            }

            double x_m = m * h;
            double t_kp1 = (k + 1) * tau;

            if (t_kp1 < x_m) {
                double x0 = x_m - t_kp1;
                u[k + 1][m] = phi(x0);
            } else {
                double x_m_half = (m - 0.5) * h;
                double t_k_half = (k + 0.5) * tau;
                double f_val = f(x_m_half, t_k_half);

                if (m == local_start || (m == next_m[rank] - local_nx && k > 0)) {
                    double u_prev_k_local = (m == 1) ? u[k][0] : u[k][m - 1];
                    u[k + 1][m] = u[k][m] - (tau / h) * (u[k][m] - u_prev_k_local) + tau * f_val;
                } else {
                    u[k + 1][m] = (h * (u[k][m] + u[k][m - 1]) - tau * (u[k][m] - u[k][m - 1]) +
                                   2 * tau * h * f_val - (h - tau) * u_prev_k1) / (h + tau);
                }
            }
            u_prev_k1 = u[k + 1][m];
            m++;
        }

        compute_end = MPI_Wtime();
        local_compute_time += (compute_end - compute_start);

        // Отправляем u[k+1][local_end-1] следующему процессору
        if (rank < size - 1 && local_end <= Nx) {
            MPI_Send(&u[k + 1][local_end - 1], 1, MPI_DOUBLE, rank + 1, k, MPI_COMM_WORLD);
        }

        // Обновляем u[k][m-1] для следующего шага
        if ((k + 1) % SYNC_INTERVAL == 0 || k == Nt - 1) {
            // Полная синхронизация
            if (rank == 0) {
                for (int i = 1; i < size; i++) {
                    int start = (i < remainder) ? i * (local_nx + 1) : i * local_nx + remainder;
                    int count = (i < remainder) ? local_nx + 1 : local_nx;
                    if (start < Nx) {
                        MPI_Recv(&u[k + 1][start], count, MPI_DOUBLE, i, k, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    }
                }
            } else {
                MPI_Send(&u[k + 1][local_start], local_size, MPI_DOUBLE, 0, k, MPI_COMM_WORLD);
            }
            MPI_Bcast(&u[k + 1][1], Nx, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        }
    }

    // Сбор данных на последнем шаге, если нужно записать в файл
    if (write_to_file) {
        for (int k = 0; k <= Nt; k++) {
            if (k % SYNC_INTERVAL != 0 && k != Nt) continue; // Уже синхронизировано
            if (rank == 0) {
                for (int i = 1; i < size; i++) {
                    int start = (i < remainder) ? i * (local_nx + 1) : i * local_nx + remainder;
                    int count = (i < remainder) ? local_nx + 1 : local_nx;
                    if (start < Nx) {
                        MPI_Recv(&u[k][start], count, MPI_DOUBLE, i, k, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    }
                }
            } else {
                MPI_Send(&u[k][local_start], local_size, MPI_DOUBLE, 0, k, MPI_COMM_WORLD);
            }
            MPI_Bcast(&u[k][1], Nx, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        }
    }

    // Синхронизация перед концом замера
    MPI_Barrier(MPI_COMM_WORLD);
    end_time = MPI_Wtime();
    total_time = end_time - start_time;

    // Сбор вычислительного времени
    MPI_Reduce(&local_compute_time, &compute_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    // Вывод времени на главном процессоре
    if (rank == 0) {
        printf("Общее время выполнения (включая коммуникации): %f секунд\n", total_time);
        printf("Примерное время коммуникаций: %f секунд\n", total_time - compute_time / size);
    }

    // Отладочный вывод
    /*if (rank == 0) {
        printf("u[0][0]=%f, u[0][1000]=%f, u[0][2000]=%f\n", u[0][0], u[0][1000], u[0][2000]);
        printf("u[100][0]=%f, u[100][1000]=%f, u[100][2000]=%f\n", u[100][0], u[100][1000], u[100][2000]);
    }*/

    // Запись результатов в файл, если не указан флаг -nf
    if (write_to_file && rank == 0) {
        FILE *output = fopen("output.txt", "w");
        for (int k = 0; k <= Nt; k++) {
            for (int m = 0; m <= Nx; m++) {
                fprintf(output, "%f %f %f\n", k * tau, m * h, u[k][m]);
            }
        }
        fclose(output);
        printf("Результаты записаны в output.txt\n");
    }

    // Освобождение памяти
    for (int k = 0; k <= Nt; k++) {
        free(u[k]);
    }
    free(u);
    free(next_m);

    MPI_Finalize();
    return 0;
}
