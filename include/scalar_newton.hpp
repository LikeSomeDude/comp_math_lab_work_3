#pragma once

#include <cstddef>

struct ScalarNewtonOptions { // параметры решения скалярного уравнения методом
                             // Ньютона
  double eps = 1e-4;         // точность
  double search_left = -20.0; // левая граница поиска корня
  double search_right = 20.0; // правая граница поиска корня
  double search_step = 0.01; // шаг поиска корня
  std::size_t max_iterations = 100; // максимальное количество итераций
};

struct ScalarNewtonResult { // результат решения скалярного уравнения методом
                            // Ньютона
  double root = 0.0;        // корень уравнения
  double residual = 0.0;    // невязка
  double segment_left = 0.0;  // левая граница сегмента
  double segment_right = 0.0; // правая граница сегмента
  std::size_t iterations = 0; // количество итераций
  bool converged = false;     // сошелся ли метод
};

// уравнение варианта 20: (x + 3)cos(x) = 1
double scalar_function(double x); // приводим к виду f(x) = (x + 3)cos(x) - 1

// производная функции f(x) = (x + 3)cos(x) - 1
double scalar_function_derivative(double x);

// находит минимальный по модулю ненулевой корень методом Ньютона
ScalarNewtonResult solve_scalar_newton(const ScalarNewtonOptions &options);
