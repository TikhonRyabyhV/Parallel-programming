#include <stdio.h>
#include <stdlib.h>
#include <math.h>
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

    // Массив для u (хранится на всех процессорах)
    double **u = (double **)malloc((Nt + 1) * sizeof(double *));
    for (int k = 0; k <= Nt; k++) {
        u[k] = (double *)calloc(Nx + 1, sizeof(double)); // Инициализация нулями
    }

    // Инициализация начального условия u(0, x) = phi(x) на всех процессорах
    for (int m = 0; m <= Nx; m++) {
        u[0][m] = phi(m * h);
    }

    // Инициализация граничного условия u(t, 0) = psi(t) на главном процессоре
    if (rank == 0) {
        for (int k = 1; k <= Nt; k++) { // Не трогаем k=0, чтобы сохранить u[0][0]
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
    double start_time, end_time, compute_time;
    double local_compute_time = 0.0;

    // Основной цикл по времени
    for (int k = 0; k < Nt; k++) {
        if (size == 1) {
            // Специальный случай для одного процесса: процессор 0 выполняет все вычисления
            start_time = MPI_Wtime(); // Начало замера времени

            double u_prev_k1 = u[k][0]; // Начальное значение для схемы
            for (int m = 1; m <= Nx; m++) {
                double x_m = m * h;
                double t_kp1 = (k + 1) * tau;

                if (t_kp1 < x_m) {
                    // Используем начальные данные: u(t_{k+1}, x_m) = u(0, x_m - t_{k+1})
                    double x0 = x_m - t_kp1;
                    u[k + 1][m] = phi(x0);
                } else {
                    // Вычисление f в точке (k+0.5, m-0.5)
                    double x_m_half = (m - 0.5) * h;
                    double t_k_half = (k + 0.5) * tau;
                    double f_val = f(x_m_half, t_k_half);

                    if (m == 1) {
                        // Схема "уголок" для первого узла
                        double u_prev_k = u[k][0];
                        u[k + 1][m] = u[k][m] - (tau / h) * (u[k][m] - u_prev_k) + tau * f_val;
                    } else {
                        // Неявная схема для внутренних узлов
                        u[k + 1][m] = (h * (u[k][m] + u[k][m - 1]) - tau * (u[k][m] - u[k][m - 1]) +
                                       2 * tau * h * f_val - (h - tau) * u_prev_k1) / (h + tau);
                    }
                }
                u_prev_k1 = u[k + 1][m]; // Обновляем для следующей итерации
            }

            end_time = MPI_Wtime(); // Конец замера времени
            local_compute_time += (end_time - start_time);
        } else {
            // Динамическая балансировка для нескольких процессов
            if (rank == 0) {
                // Менеджер: распределяет задачи
                start_time = MPI_Wtime(); // Начало замера времени (только управление задачами)

                int next_m = 1; // Начинаем с m=1, так как m=0 — граничное условие
                int active_workers = size - 1; // Число активных рабочих процессов

                while (active_workers > 0) {
                    // Ожидаем запрос от рабочего процесса
                    MPI_Status status;
                    int worker_rank;
                    MPI_Recv(&worker_rank, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);

                    if (next_m <= Nx) {
                        // Отправляем задачу: начальный индекс и размер блока
                        int block_size = (next_m + BLOCK_SIZE - 1 <= Nx) ? BLOCK_SIZE : (Nx - next_m + 1);
                        int task[2] = {next_m, block_size};
                        MPI_Send(task, 2, MPI_INT, worker_rank, 1, MPI_COMM_WORLD);

                        // Получаем результаты от рабочего
                        double *recv_buf = (double *)malloc(block_size * sizeof(double));
                        MPI_Recv(recv_buf, block_size, MPI_DOUBLE, worker_rank, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                        // Сохраняем результаты в u[k+1]
                        for (int m = 0; m < block_size; m++) {
                            u[k + 1][next_m + m] = recv_buf[m];
                        }
                        free(recv_buf);

                        next_m += block_size;
                    } else {
                        // Больше задач нет, отправляем сигнал завершения
                        int task[2] = {-1, 0};
                        MPI_Send(task, 2, MPI_INT, worker_rank, 1, MPI_COMM_WORLD);
                        active_workers--;
                    }
                }

                end_time = MPI_Wtime(); // Конец замера времени
                local_compute_time += (end_time - start_time);
            } else {
                // Рабочий процесс
                while (1) {
                    // Запрашиваем задачу у менеджера
                    MPI_Send(&rank, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);

                    // Получаем задачу
                    int task[2];
                    MPI_Recv(task, 2, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    int start_m = task[0];
                    int block_size = task[1];

                    if (start_m == -1) {
                        // Сигнал завершения
                        break;
                    }

                    // Замеряем время только для вычислений
                    start_time = MPI_Wtime();

                    // Вычисляем u[k+1][m] для данного блока
                    double *local_u = (double *)malloc(block_size * sizeof(double));
                    double u_prev_k1 = (start_m == 1) ? u[k][0] : u[k][start_m - 1]; // Начальное значение для схемы

                    for (int i = 0; i < block_size; i++) {
                        int m = start_m + i;
                        double x_m = m * h;
                        double t_kp1 = (k + 1) * tau;

                        if (t_kp1 < x_m) {
                            // Используем начальные данные: u(t_{k+1}, x_m) = u(0, x_m - t_{k+1})
                            double x0 = x_m - t_kp1;
                            local_u[i] = phi(x0);
                        } else {
                            // Вычисление f в точке (k+0.5, m-0.5)
                            double x_m_half = (m - 0.5) * h;
                            double t_k_half = (k + 0.5) * tau;
                            double f_val = f(x_m_half, t_k_half);

                            if (i == 0) {
                                // Схема "уголок" для первого узла блока
                                double u_prev_k = (m == 1) ? u[k][0] : u[k][m - 1];
                                local_u[i] = u[k][m] - (tau / h) * (u[k][m] - u_prev_k) + tau * f_val;
                            } else {
                                // Неявная схема для внутренних узлов блока
                                local_u[i] = (h * (u[k][m] + u[k][m - 1]) - tau * (u[k][m] - u[k][m - 1]) +
                                              2 * tau * h * f_val - (h - tau) * u_prev_k1) / (h + tau);
                            }
                        }
                        u_prev_k1 = local_u[i]; // Обновляем для следующей итерации
                    }

                    end_time = MPI_Wtime();
                    local_compute_time += (end_time - start_time);

                    // Отправляем результаты менеджеру
                    MPI_Send(local_u, block_size, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD);
                    free(local_u);
                }
            }
        }

        // Синхронизация u[k+1][m] между всеми процессорами после шага k
        for (int m = 1; m <= Nx; m++) {
            MPI_Bcast(&u[k + 1][m], 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        }
    }

    // Сбор времени вычислений со всех процессов
    MPI_Reduce(&local_compute_time, &compute_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    // Вывод времени на главном процессоре
    if (rank == 0) {
        printf("Время вычислений: %f секунд\n", compute_time);
    }

    // Отладочный вывод
    if (rank == 0) {
        printf("u[0][0]=%f, u[0][1000]=%f, u[0][2000]=%f\n", u[0][0], u[0][1000], u[0][2000]);
        printf("u[100][0]=%f, u[100][1000]=%f, u[100][2000]=%f\n", u[100][0], u[100][1000], u[100][2000]);
    }

    // Запись результатов в файл на главном процессоре
    if (rank == 0) {
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

    MPI_Finalize();
    return 0;
}
