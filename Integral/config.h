#pragma once

// Параметры интегрирования
const double from = 1e-3;
const double to = 1;
double epsilon = 1e-6; // Погрешность по умолчанию
const int max_stack_size = 1000;
const int local_stack_critical_line = 10;

// Функция для интегрирования:
//double f(double x) {
//    return cos(1.0 / (7.0 - x));
//}
//
double f(double x) {
	return 1 / (x * x) * sin(1 / (x * x));
}
