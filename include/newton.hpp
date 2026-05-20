#pragma once

#include "types.hpp"

#include <cstddef>
#include <functional>

using SystemFunction = std::function<Vector(const Vector &)>;
using JacobianFunction = std::function<Matrix(const Vector &)>;

// параметры единой реализации метода Ньютона
struct NewtonOptions {
  double eps = 1e-10;
  std::size_t max_iterations = 500;

  // первые k итераций выполняются как полный метод Ньютона
  // матрица Якоби пересчитывается на каждой итерации
  // k = 1, m = infinity соответствует модифицированному методу
  std::size_t full_jacobian_iterations = 1;

  // после первых k итераций матрица Якоби пересчитывается каждые m итераций
  // m = 1, k = infinity соответствует полному методу Ньютона
  // std::numeric_limits<std::size_t>::max() используется как infinity
  std::size_t jacobian_recalculation_period = 1;
};

// результат одного запуска метода Ньютона
struct NewtonResult {
  Vector x;                   // решение
  std::size_t iterations = 0; // количество итераций
  std::size_t jacobian_recalculations = 0; // счетчик пересчетов Якоби
  bool converged = false;                  // сошелся ли метод
  double residual_norm = 0.0;              // норма невязки
  OperationStats linear_solve_stats; // статистика решения СЛАУ
};

// решает произвольную систему F(x) = 0 методом Ньютона с параметрами m и k
NewtonResult solve_newton(const Vector &x0, const SystemFunction &values,
                          const JacobianFunction &jacobian,
                          const NewtonOptions &options);

// решает основную 10-мерную систему из задания
NewtonResult solve_newton(const Vector &x0, const NewtonOptions &options);

// бесконечная норма вектора max_i |v_i|
double norm_inf(const Vector &v);
