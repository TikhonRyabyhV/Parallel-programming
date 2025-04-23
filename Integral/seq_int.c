#include <stdio.h>
#include <stdlib.h>
#include <math.h>
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

// Интегрирование методом трапеций с адаптивным разбиением
double integrate(double a, double b, double eps, int num_segments) {
    Stack stack;
    init_stack(&stack, max_stack_size);
    double total_sum = 0.0;
    double total_length = b - a;

    // Начальное разбиение интервала
    double h = (b - a) / num_segments;
    for (int i = 0; i < num_segments; i++) {
        double start = a + i * h;
        double end = (i == num_segments - 1) ? b : start + h;
        double fa = f(start);
        double fb = f(end);
        double sum = (end - start) * (fa + fb) / 2.0;
        PUSH(stack, start, end, fa, fb, sum);
    }

    while (IS_NOT_EMPTY(stack)) {
        double a, b, fa, fb, sum;
        POP(stack, a, b, fa, fb, sum);

        double mid = (a + b) / 2.0;
        double h = (b - a) / 2.0;
        double fm = f(mid);
        double sum_left = h * (fa + fm) / 2.0;
        double sum_right = h * (fm + fb) / 2.0;
        double sum_new = sum_left + sum_right;

        if (fabs(sum_new - sum) < eps * fabs(sum_new)) {
            total_sum += sum_new;
        } else {
            if (stack.top + 2 >= stack.max_size) {
                printf("Ошибка: стек переполнен!\n");
                free_stack(&stack);
                return 0.0;
            }
            PUSH(stack, mid, b, fm, fb, sum_right);
            PUSH(stack, a, mid, fa, fm, sum_left);
        }
    }

    free_stack(&stack);
    return total_sum;
}

int main(int argc, char *argv[]) {
    if (argc != 3) {
        printf("Использование: %s <число_отрезков> <погрешность>\n", argv[0]);
        return 1;
    }

    int num_segments = atoi(argv[1]);
    epsilon = atof(argv[2]);

    if (num_segments < 1) {
        printf("Ошибка: число отрезков должно быть положительным\n");
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

    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);

    double result = integrate(from, to, epsilon, num_segments);
    
    clock_gettime(CLOCK_MONOTONIC, &end);
    double elapsed = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("I = %.10f\n", result);
    printf("Segments: %d, time: %.6f s\n", num_segments, elapsed);

    return 0;
}
