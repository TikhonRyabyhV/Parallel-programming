#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "config.h"

// Структура для отрезка [a, b] и его вклада в интеграл
typedef struct {
    double a;      // левая граница
    double b;      // правая граница
    double fa;     // значение функции в a
    double fb;     // значение функции в b
    double sum;    // приближение интеграла
} Segment;

// Структура для стека
typedef struct {
    Segment *data; // массив отрезков
    int top;       // указатель вершины стека
    int max_size;  // максимальный размер стека
} Stack;

// Макросы для работы со стеком
#define IS_EMPTY(stack) ((stack).top == 0)
#define IS_NOT_EMPTY(stack) ((stack).top > 0)

#define PUSH(stack, nA, nB, nfA, nfB, nsAB) \
    do { \
        (stack).data[(stack).top].a = (nA); \
        (stack).data[(stack).top].b = (nB); \
        (stack).data[(stack).top].fa = (nfA); \
        (stack).data[(stack).top].fb = (nfB); \
        (stack).data[(stack).top].sum = (nsAB); \
        (stack).top++; \
    } while (0)

#define POP(stack, nA, nB, nfA, nfB, nsAB) \
    do { \
        (stack).top--; \
        nA = (stack).data[(stack).top].a; \
        nB = (stack).data[(stack).top].b; \
        nfA = (stack).data[(stack).top].fa; \
        nfB = (stack).data[(stack).top].fb; \
        nsAB = (stack).data[(stack).top].sum; \
    } while (0)

// Глобальные переменные
pthread_mutex_t global_stack_mutex;
pthread_mutex_t global_sum_mutex;
pthread_mutex_t task_present_mutex;
Stack global_stack = {0};
double global_sum = 0.0;
int active_handlers_count = 0;

// Параметры потока
typedef struct {
    int threads_count; // общее число потоков
    double eps;       // погрешность
    double total_length; // длина интервала интегрирования
} ThreadArgs;

// Инициализация стека
void init_stack(Stack *stack, int max_size) {
    stack->data = (Segment *)malloc(max_size * sizeof(Segment));
    stack->top = 0;
    stack->max_size = max_size;
}

// Освобождение стека
void free_stack(Stack *stack) {
    free(stack->data);
}

// Функция обработки потока
void *handle_global_stack(void *arg) {
    ThreadArgs *args = (ThreadArgs *)arg;
    int threads_count = args->threads_count;
    double eps = args->eps;
    double total_length = args->total_length;

    Stack local_stack;
    init_stack(&local_stack, max_stack_size);
    double s = 0.0;
    double a, b, fa, fb, sum;

    while (1) {
        pthread_mutex_lock(&task_present_mutex);
        pthread_mutex_lock(&global_stack_mutex);

        if (IS_EMPTY(global_stack)) {
            pthread_mutex_unlock(&global_stack_mutex);
            pthread_mutex_unlock(&task_present_mutex);
            break;
        }

        POP(global_stack, a, b, fa, fb, sum);
        if (IS_NOT_EMPTY(global_stack)) {
            pthread_mutex_unlock(&task_present_mutex);
        }
        if (a < b) active_handlers_count++;
        pthread_mutex_unlock(&global_stack_mutex);

        if (a > b) { // терминальный отрезок
            pthread_mutex_unlock(&task_present_mutex);
            break;
        }

        while (1) {
            double mid = (a + b) / 2.0;
            double fm = f(mid);
            double h = (mid - a) / 2.0;
            double sum_left = (fa + fm) * h;
            double sum_right = (fm + fb) * h;
            double sum_new = sum_left + sum_right;

            if (fabs(sum_new - sum) < eps * fabs(sum_new)) {
                s += sum_new;
                if (IS_EMPTY(local_stack)) break;
                POP(local_stack, a, b, fa, fb, sum);
            } else {
                if (local_stack.top >= local_stack.max_size) {
                    printf("Ошибка: локальный стек переполнен!\n");
                    break;
                }
                PUSH(local_stack, a, mid, fa, fm, sum_left);
                a = mid;
                fa = fm;
                sum = sum_right;
            }

            if (local_stack.top >= local_stack_critical_line) {
                pthread_mutex_lock(&global_stack_mutex);
                if (IS_NOT_EMPTY(global_stack)) {
                    pthread_mutex_unlock(&global_stack_mutex);
                    continue;
                }
                if (global_stack.top + local_stack_critical_line > global_stack.max_size) {
                    printf("Ошибка: глобальный стек переполнен!\n");
                    pthread_mutex_unlock(&global_stack_mutex);
                    break;
                }
                memcpy(global_stack.data + global_stack.top,
                       local_stack.data,
                       sizeof(Segment) * local_stack_critical_line);
                global_stack.top += local_stack_critical_line;
                local_stack.top = 0;
                pthread_mutex_unlock(&task_present_mutex);
                pthread_mutex_unlock(&global_stack_mutex);
            }
        }

        pthread_mutex_lock(&global_stack_mutex);
        active_handlers_count--;
        if (active_handlers_count == 0 && IS_EMPTY(global_stack)) {
            for (int i = 0; i < threads_count; i++) {
                PUSH(global_stack, 1.0, 0.0, 0.0, 0.0, 0.0); // терминальный отрезок
            }
            pthread_mutex_unlock(&task_present_mutex);
        }
        pthread_mutex_unlock(&global_stack_mutex);
    }

    pthread_mutex_lock(&global_sum_mutex);
    global_sum += s;
    pthread_mutex_unlock(&global_sum_mutex);

    free_stack(&local_stack);
    return NULL;
}

int main(int argc, char *argv[]) {
    if (argc != 3) {
        printf("Использование: %s <число_потоков> <погрешность>\n", argv[0]);
        return 1;
    }

    int num_threads = atoi(argv[1]);
    epsilon = atof(argv[2]);

    if (num_threads < 1) {
        printf("Ошибка: число потоков должно быть положительным\n");
        return 1;
    }
    if (epsilon <= 0) {
        printf("Ошибка: погрешность должна быть положительной\n");
        return 1;
    }
    if (to >= 7.0) {
        printf("Ошибка: функция не определена в точке x=7\n");
        return 1;
    }

    pthread_mutex_init(&global_stack_mutex, NULL);
    pthread_mutex_init(&global_sum_mutex, NULL);
    pthread_mutex_init(&task_present_mutex, NULL);

    init_stack(&global_stack, max_stack_size);
    PUSH(global_stack, from, to, f(from), f(to), (f(from) + f(to)) * (to - from) / 2.0);

    pthread_t *threads = (pthread_t *)malloc(num_threads * sizeof(pthread_t));
    ThreadArgs *args = (ThreadArgs *)malloc(num_threads * sizeof(ThreadArgs));

    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);

    for (int i = 0; i < num_threads; i++) {
        args[i].threads_count = num_threads;
        args[i].eps = epsilon;
        args[i].total_length = to - from;
        pthread_create(&threads[i], NULL, handle_global_stack, &args[i]);
    }

    for (int i = 0; i < num_threads; i++) {
        pthread_join(threads[i], NULL);
    }

    clock_gettime(CLOCK_MONOTONIC, &end);
    double elapsed = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("I = %.10f\n", global_sum);
    printf("Threads: %d, time: %.6f s\n", num_threads, elapsed);

    free_stack(&global_stack);
    free(threads);
    free(args);
    pthread_mutex_destroy(&global_stack_mutex);
    pthread_mutex_destroy(&global_sum_mutex);
    pthread_mutex_destroy(&task_present_mutex);

    return 0;
}
