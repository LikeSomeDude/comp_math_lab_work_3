#include "small_system.hpp"

#include <cmath>
#include <stdexcept>

namespace {
void validate_xy_size(const Vector &x) {
  if (x.size() != 2) {
    throw std::invalid_argument("Expected vector of size 2.");
  }
}
} // namespace

Vector small_system_initial_guess() { // начальное приближение
  return {-1.0, -1.0};
}

Vector small_system_values(const Vector &x) { // возвращает F(x)
  validate_xy_size(x);

  const double x_value = x[0];
  const double y_value = x[1];

  return {
      std::sin(y_value + 0.5) - x_value - 1.0,
      y_value + std::cos(x_value - 2.0),
  };
}

Matrix small_system_jacobian(const Vector &x) { // возвращает J(x)
  validate_xy_size(x);

  const double x_value = x[0];
  const double y_value = x[1];

  return {
      {-1.0, std::cos(y_value + 0.5)},
      {-std::sin(x_value - 2.0), 1.0},
  };
}
