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
#define BLOCK_SIZE 100 // Размер блока для динамической балансировки

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

    // Переменные для замера времени
    double start_time, end_time, total_time;
    double compute_start, compute_end, local_compute_time = 0.0, compute_time;

    // Переменная для текущего индекса m (глобальная для всех процессоров)
    int global_next_m = 1; // Начинаем с m=1, так как m=0 — граничное условие

    // Открываем файл для записи с помощью MPI-IO
    MPI_File fh;
    if (write_to_file) {
        MPI_File_open(MPI_COMM_WORLD, "output.txt", MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
        // Очистка файла на процессоре 0
        if (rank == 0) {
            MPI_File_set_size(fh, 0);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    // Синхронизация перед началом замера
    MPI_Barrier(MPI_COMM_WORLD);
    start_time = MPI_Wtime();

    // Основной цикл по времени
    for (int k = 0; k < Nt; k++) {
        // Каждый процессор получает блок индексов m
        int local_start_m = -1;
        int local_block_size = 0;
        double u_prev_k1 = u[k][0]; // Начальное значение для m=1

        // Динамическая балансировка
        while (global_next_m <= Nx) {
            // Критическая секция: атомарно получаем следующий блок
            if (rank == 0) {
                local_start_m = global_next_m;
                local_block_size = (local_start_m + BLOCK_SIZE - 1 <= Nx) ? BLOCK_SIZE : (Nx - local_start_m + 1);
                global_next_m += local_block_size;

                // Рассылаем информацию о новом global_next_m
                for (int i = 1; i < size; i++) {
                    MPI_Send(&global_next_m, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
                }
            } else {
                MPI_Recv(&global_next_m, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                local_start_m = global_next_m - local_block_size;
                local_block_size = (local_start_m + BLOCK_SIZE - 1 <= Nx) ? BLOCK_SIZE : (Nx - local_start_m + 1);
            }

            if (local_start_m > Nx) break;

            // Получаем u[k][m-1] для первого индекса блока
            if (local_start_m > 1) {
                if (rank == 0) {
                    u_prev_k1 = u[k][local_start_m - 1];
                } else {
                    MPI_Recv(&u_prev_k1, 1, MPI_DOUBLE, rank - 1, k, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }

            // Вычисления для текущего блока
            compute_start = MPI_Wtime();

            for (int m = local_start_m; m < local_start_m + local_block_size; m++) {
                double x_m = m * h;
                double t_kp1 = (k + 1) * tau;

                if (t_kp1 < x_m) {
                    double x0 = x_m - t_kp1;
                    u[k + 1][m] = phi(x0);
                } else {
                    double x_m_half = (m - 0.5) * h;
                    double t_k_half = (k + 0.5) * tau;
                    double f_val = f(x_m_half, t_k_half);

                    if (m == local_start_m) {
                        double u_prev_k_local = (m == 1) ? u[k][0] : u[k][m - 1];
                        u[k + 1][m] = u[k][m] - (tau / h) * (u[k][m] - u_prev_k_local) + tau * f_val;
                    } else {
                        u[k + 1][m] = (h * (u[k][m] + u[k][m - 1]) - tau * (u[k][m] - u[k][m - 1]) +
                                       2 * tau * h * f_val - (h - tau) * u_prev_k1) / (h + tau);
                    }
                }
                u_prev_k1 = u[k + 1][m];
            }

            compute_end = MPI_Wtime();
            local_compute_time += (compute_end - compute_start);

            // Отправляем u[k+1][m] последнему процессору в цепочке
            if (rank < size - 1 && local_start_m + local_block_size - 1 <= Nx) {
                MPI_Send(&u[k + 1][local_start_m + local_block_size - 1], 1, MPI_DOUBLE, rank + 1, k, MPI_COMM_WORLD);
            }

            // Записываем результаты в файл с соответствующим смещением
            if (write_to_file) {
                for (int m = local_start_m; m < local_start_m + local_block_size; m++) {
                    char line[256];
                    int len = snprintf(line, sizeof(line), "%f %f %f\n", (k + 1) * tau, m * h, u[k + 1][m]);
                    MPI_Offset offset = ((k + 1) * (Nx + 1) + m) * len;
                    MPI_File_write_at(fh, offset, line, len, MPI_CHAR, MPI_STATUS_IGNORE);
                }
            }

            // Обновляем u[k][m] для следующего шага k
            if (rank == 0) {
                for (int i = 1; i < size; i++) {
                    int start_m, block_size;
                    MPI_Recv(&start_m, 1, MPI_INT, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Recv(&block_size, 1, MPI_INT, i, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    double *recv_buf = (double *)malloc(block_size * sizeof(double));
                    MPI_Recv(recv_buf, block_size, MPI_DOUBLE, i, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    for (int m = 0; m < block_size; m++) {
                        u[k + 1][start_m + m] = recv_buf[m];
                    }
                    free(recv_buf);
                }
                // Рассылаем только необходимые данные для следующего шага
                for (int i = 1; i < size; i++) {
                    int start_m = i * (Nx / size);
                    int end_m = (i + 1) * (Nx / size);
                    if (end_m > Nx) end_m = Nx;
                    int count = end_m - start_m + 1;
                    MPI_Send(&u[k + 1][start_m], count, MPI_DOUBLE, i, 4, MPI_COMM_WORLD);
                }
            } else {
                MPI_Send(&local_start_m, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
                MPI_Send(&local_block_size, 1, MPI_INT, 0, 2, MPI_COMM_WORLD);
                MPI_Send(&u[k + 1][local_start_m], local_block_size, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD);
                // Получаем данные, необходимые для следующего шага
                int start_m = rank * (Nx / size);
                int end_m = (rank + 1) * (Nx / size);
                if (end_m > Nx) end_m = Nx;
                int count = end_m - start_m + 1;
                MPI_Recv(&u[k + 1][start_m], count, MPI_DOUBLE, 0, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }

        // Сбрасываем global_next_m для следующего шага k
        if (rank == 0) {
            global_next_m = 1;
            for (int i = 1; i < size; i++) {
                MPI_Send(&global_next_m, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            }
        } else {
            MPI_Recv(&global_next_m, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }

    // Закрываем файл
    if (write_to_file) {
        MPI_File_close(&fh);
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
        printf("Суммарное вычислительное время (по всем процессорам): %f секунд\n", compute_time);
        printf("Примерное время коммуникаций: %f секунд\n", total_time - compute_time / size);
    }

    // Отладочный вывод
    if (rank == 0) {
        printf("u[0][0]=%f, u[0][1000]=%f, u[0][2000]=%f\n", u[0][0], u[0][1000], u[0][2000]);
        printf("u[100][0]=%f, u[100][1000]=%f, u[100][2000]=%f\n", u[100][0], u[100][1000], u[100][2000]);
    }

    // Освобождение памяти
    for (int k = 0; k <= Nt; k++) {
        free(u[k]);
    }
    free(u);

    MPI_Finalize();
    return 0;
}
