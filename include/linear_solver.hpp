#pragma once

#include "types.hpp"

// результат LU-разложения с частичным выбором главного элемента
struct LuFactorization {
  Matrix l;
  Matrix u;
  std::vector<std::size_t> permutation; // перестановка
  OperationStats stats;                 // статистика
};

// решение одной СЛАУ и число операций
struct LinearSolveResult {
  Vector x;
  OperationStats stats;
};

// строит LU-разложение матрицы
LuFactorization decompose_lu(Matrix a);

// решает СЛАУ по уже готовому LU-разложению
// это нужно для модифицированного метода Ньютона, где матрица Якоби не
// меняется
LinearSolveResult solve_with_lu(const LuFactorization &lu, const Vector &b);

// сразу строит LU-разложение и решает СЛАУ
LinearSolveResult solve_linear_system_lu(Matrix a, Vector b);
