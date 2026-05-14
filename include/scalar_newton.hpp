#pragma once

#include <cstddef>

struct ScalarNewtonOptions {
    double eps = 1e-4;
    double search_left = -20.0;
    double search_right = 20.0;
    double search_step = 0.01;
    std::size_t max_iterations = 100;
};

struct ScalarNewtonResult {
    double root = 0.0;
    double residual = 0.0;
    double segment_left = 0.0;
    double segment_right = 0.0;
    std::size_t iterations = 0;
    bool converged = false;
};

// Уравнение варианта 20: (x + 3)cos(x) = 1.
double scalar_function(double x);

// Производная функции f(x) = (x + 3)cos(x) - 1.
double scalar_function_derivative(double x);

// Находит минимальный по модулю ненулевой корень методом Ньютона.
ScalarNewtonResult solve_scalar_newton(const ScalarNewtonOptions& options);
