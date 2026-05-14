#pragma once

#include "types.hpp"

#include <cstddef>

// Параметры единой реализации метода Ньютона.
struct NewtonOptions {
    double eps = 1e-10;
    std::size_t max_iterations = 500;

    // Первые k итераций выполняются как полный метод Ньютона:
    // матрица Якоби пересчитывается на каждой итерации.
    // k = 1, m = infinity соответствует модифицированному методу.
    std::size_t full_jacobian_iterations = 1;

    // После первых k итераций матрица Якоби пересчитывается каждые m итераций.
    // m = 1, k = infinity соответствует полному методу Ньютона.
    // std::numeric_limits<std::size_t>::max() используется как infinity.
    std::size_t jacobian_recalculation_period = 1;
};

// Результат одного запуска метода Ньютона.
struct NewtonResult {
    Vector x;
    std::size_t iterations = 0;
    std::size_t jacobian_recalculations = 0;
    bool converged = false;
    double residual_norm = 0.0;
    OperationStats linear_solve_stats;
};

// Решает систему F(x) = 0 методом Ньютона с параметрами m и k.
NewtonResult solve_newton(const Vector& x0, const NewtonOptions& options);

// Бесконечная норма вектора: max_i |v_i|.
double norm_inf(const Vector& v);
