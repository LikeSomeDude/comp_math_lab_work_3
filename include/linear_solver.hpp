#pragma once

#include "types.hpp"

// Результат LU-разложения с частичным выбором главного элемента.
struct LuFactorization {
    Matrix l;
    Matrix u;
    std::vector<std::size_t> permutation;
    OperationStats stats;
};

// Решение одной СЛАУ и число операций, потраченных на это решение.
struct LinearSolveResult {
    Vector x;
    OperationStats stats;
};

// Строит LU-разложение матрицы A = P^(-1)LU.
LuFactorization decompose_lu(Matrix a);

// Решает СЛАУ по уже готовому LU-разложению.
// Это нужно для модифицированного метода Ньютона, где матрица Якоби не меняется.
LinearSolveResult solve_with_lu(const LuFactorization& lu, const Vector& b);

// Удобная обертка: сразу строит LU-разложение и решает СЛАУ.
LinearSolveResult solve_linear_system_lu(Matrix a, Vector b);
