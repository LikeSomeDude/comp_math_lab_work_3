#include "newton.hpp"

#include "linear_solver.hpp"
#include "nonlinear_system.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace {
// возвращает вектор с противоположными знаками
Vector negate(const Vector &v) {
  Vector result(v.size());
  for (std::size_t i = 0; i < v.size(); ++i) {
    result[i] = -v[i];
  }
  return result;
}

// прибавляет rhs к lhs
void add_in_place(Vector &lhs, const Vector &rhs) {
  for (std::size_t i = 0; i < lhs.size(); ++i) {
    lhs[i] += rhs[i];
  }
}

// возвращает true, если нужно пересчитать матрицу Якоби
bool should_recalculate_jacobian(std::size_t iteration,
                                 const NewtonOptions &options) {
  if (iteration == 0) {
    return true;
  }

  // первые k раз пересчитываем каждую итерацию
  if (iteration < options.full_jacobian_iterations) {
    return true;
  }

  // не пересчитываем если m = infinity
  if (options.jacobian_recalculation_period == 0 ||
      options.jacobian_recalculation_period ==
          std::numeric_limits<std::size_t>::max()) {
    return false;
  }

  // после первых k итераций пересчитываем каждые m
  return (iteration - options.full_jacobian_iterations + 1) %
             options.jacobian_recalculation_period ==
         0;
}
} // namespace

// бесконечная норма вектора
double norm_inf(const Vector &v) {
  double result = 0.0;
  for (double value : v) {
    result = std::max(result, std::abs(value));
  }
  return result;
}

// обертка для короткого вызова метода на 10-мерной системе
NewtonResult solve_newton(const Vector &x0, const NewtonOptions &options) {
  return solve_newton(x0, system_values, jacobian, options);
}

// универсальная функция решения системы методом Ньютона
NewtonResult solve_newton(const Vector &x0, const SystemFunction &values,
                          const JacobianFunction &jacobian_function,
                          const NewtonOptions &options) {
  NewtonResult result;
  result.x = x0; // начальное приближение

  LuFactorization current_lu; // LU-разложение текущей матрицы Якоби

  for (std::size_t iteration = 0; iteration < options.max_iterations;
       ++iteration) {
    const Vector f = values(result.x);  // вектор функций
    result.residual_norm = norm_inf(f); // норма невязки

    if (result.residual_norm < options.eps) { // проверка точности
      result.converged = true;
      result.iterations = iteration;
      return result;
    }

    // нужно ли пересчитать Якоби
    if (iteration == 0 || should_recalculate_jacobian(iteration, options)) {
      // при пересчете Якоби сразу строим LU-разложение и потом
      // переиспользуем его, пока матрица Якоби остается той же
      current_lu = decompose_lu(jacobian_function(result.x));

      ++result.jacobian_recalculations;

      result.linear_solve_stats.additions += current_lu.stats.additions;
      result.linear_solve_stats.multiplications +=
          current_lu.stats.multiplications;
      result.linear_solve_stats.divisions += current_lu.stats.divisions;
    }

    // шаг Ньютона находится из системы J(x_k) * dx = -F(x_k)
    LinearSolveResult step = solve_with_lu(current_lu, negate(f));
    result.linear_solve_stats.additions += step.stats.additions;
    result.linear_solve_stats.multiplications += step.stats.multiplications;
    result.linear_solve_stats.divisions += step.stats.divisions;

    add_in_place(result.x, step.x);
    result.iterations = iteration + 1;

    if (norm_inf(step.x) < options.eps) {
      result.residual_norm = norm_inf(values(result.x));
      result.converged = result.residual_norm < options.eps;
      return result;
    }
  }

  result.residual_norm = norm_inf(values(result.x));
  return result;
}
