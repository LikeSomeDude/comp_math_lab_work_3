#include "nonlinear_system.hpp"

#include <cmath>
#include <stdexcept>

namespace {
void validate_x_size(const Vector &x) {
  if (x.size() != 10) {
    throw std::invalid_argument("Expected vector of size 10.");
  }
}
} // namespace

Vector initial_guess() { // начальное приближение
  return {0.5, 0.5, 1.5, -1.0, -0.5, 1.5, 0.5, -0.5, 1.5, -1.5};
}

Vector changed_initial_guess() { // измененное начальное приближение
  Vector x = initial_guess();
  x[4] = -0.2;
  return x;
}

Vector system_values(const Vector &x) { // считает значения функций системы
  validate_x_size(x);

  const double x1 = x[0];
  const double x2 = x[1];
  const double x3 = x[2];
  const double x4 = x[3];
  const double x5 = x[4];
  const double x6 = x[5];
  const double x7 = x[6];
  const double x8 = x[7];
  const double x9 = x[8];
  const double x10 = x[9];

  Vector f(10);
  f[0] = std::cos(x2 * x1) - std::exp(-(3.0 * x3)) + x4 * x5 * x5 - x6 -
         std::sinh(2.0 * x8) * x9 + 2.0 * x10 + 2.000433974165385440;
  f[1] = std::sin(x2 * x1) + x3 * x9 * x7 - std::exp(-x10 + x6) +
         3.0 * x5 * x5 - x6 * (x8 + 1.0) + 10.886272036407019994;
  f[2] =
      x1 - x2 + x3 - x4 + x5 - x6 + x7 - x8 + x9 - x10 - 3.1361904761904761904;
  f[3] = 2.0 * std::cos(-x9 + x4) + x5 / (x3 + x1) - std::sin(x2 * x2) +
         std::pow(std::cos(x7 * x10), 2.0) - x8 - 0.1707472705022304757;
  f[4] = std::sin(x5) + 2.0 * x8 * (x3 + x1) - std::exp(-x7 * (-x10 + x6)) +
         2.0 * std::cos(x2) - 1.0 / (-x9 + x4) - 0.3685896273101277862;
  f[5] = std::exp(x1 - x4 - x9) + x5 * x5 / x8 +
         std::cos(3.0 * x10 * x2) / 2.0 - x6 * x3 + 2.0491086016771875115;
  f[6] = std::pow(x2, 3.0) * x7 - std::sin(x10 / x5 + x8) +
         (x1 - x6) * std::cos(x4) + x3 - 0.7380430076202798014;
  f[7] = x5 * std::pow(x1 - 2.0 * x6, 2.0) - 2.0 * std::sin(-x9 + x3) +
         1.5 * x4 - std::exp(x2 * x7 + x10) + 3.5668321989693809040;
  f[8] = 7.0 / x6 + std::exp(x5 + x4) - 2.0 * x2 * x8 * x10 * x7 + 3.0 * x9 -
         3.0 * x1 - 8.4394734508383257499;
  f[9] = x10 * x1 + x9 * x2 - x8 * x3 + std::sin(x4 + x5 + x6) * x7 -
         0.78238095238095238096;

  return f;
}

Matrix jacobian(const Vector &x) { // якобиан
  validate_x_size(x);

  const double x1 = x[0];
  const double x2 = x[1];
  const double x3 = x[2];
  const double x4 = x[3];
  const double x5 = x[4];
  const double x6 = x[5];
  const double x7 = x[6];
  const double x8 = x[7];
  const double x9 = x[8];
  const double x10 = x[9];

  Matrix j(10, Vector(10));

  j[0][0] = -x2 * std::sin(x2 * x1);
  j[0][1] = -x1 * std::sin(x2 * x1);
  j[0][2] = 3.0 * std::exp(-(3.0 * x3));
  j[0][3] = x5 * x5;
  j[0][4] = 2.0 * x4 * x5;
  j[0][5] = -1.0;
  j[0][6] = 0.0;
  j[0][7] = -2.0 * std::cosh(2.0 * x8) * x9;
  j[0][8] = -std::sinh(2.0 * x8);
  j[0][9] = 2.0;
  j[1][0] = x2 * std::cos(x2 * x1);
  j[1][1] = x1 * std::cos(x2 * x1);
  j[1][2] = x9 * x7;
  j[1][3] = 0.0;
  j[1][4] = 6.0 * x5;
  j[1][5] = -std::exp(-x10 + x6) - x8 - 1.0;
  j[1][6] = x3 * x9;
  j[1][7] = -x6;
  j[1][8] = x3 * x7;
  j[1][9] = std::exp(-x10 + x6);
  j[2][0] = 1.0;
  j[2][1] = -1.0;
  j[2][2] = 1.0;
  j[2][3] = -1.0;
  j[2][4] = 1.0;
  j[2][5] = -1.0;
  j[2][6] = 1.0;
  j[2][7] = -1.0;
  j[2][8] = 1.0;
  j[2][9] = -1.0;
  j[3][0] = -x5 * std::pow(x3 + x1, -2.0);
  j[3][1] = -2.0 * x2 * std::cos(x2 * x2);
  j[3][2] = -x5 * std::pow(x3 + x1, -2.0);
  j[3][3] = -2.0 * std::sin(-x9 + x4);
  j[3][4] = 1.0 / (x3 + x1);
  j[3][5] = 0.0;
  j[3][6] = -2.0 * std::cos(x7 * x10) * x10 * std::sin(x7 * x10);
  j[3][7] = -1.0;
  j[3][8] = 2.0 * std::sin(-x9 + x4);
  j[3][9] = -2.0 * std::cos(x7 * x10) * x7 * std::sin(x7 * x10);
  j[4][0] = 2.0 * x8;
  j[4][1] = -2.0 * std::sin(x2);
  j[4][2] = 2.0 * x8;
  j[4][3] = std::pow(-x9 + x4, -2.0);
  j[4][4] = std::cos(x5);
  j[4][5] = x7 * std::exp(-x7 * (-x10 + x6));
  j[4][6] = -(x10 - x6) * std::exp(-x7 * (-x10 + x6));
  j[4][7] = 2.0 * x3 + 2.0 * x1;
  j[4][8] = -std::pow(-x9 + x4, -2.0);
  j[4][9] = -x7 * std::exp(-x7 * (-x10 + x6));
  j[5][0] = std::exp(x1 - x4 - x9);
  j[5][1] = -3.0 / 2.0 * x10 * std::sin(3.0 * x10 * x2);
  j[5][2] = -x6;
  j[5][3] = -std::exp(x1 - x4 - x9);
  j[5][4] = 2.0 * x5 / x8;
  j[5][5] = -x3;
  j[5][6] = 0.0;
  j[5][7] = -x5 * x5 * std::pow(x8, -2.0);
  j[5][8] = -std::exp(x1 - x4 - x9);
  j[5][9] = -3.0 / 2.0 * x2 * std::sin(3.0 * x10 * x2);
  j[6][0] = std::cos(x4);
  j[6][1] = 3.0 * x2 * x2 * x7;
  j[6][2] = 1.0;
  j[6][3] = -(x1 - x6) * std::sin(x4);
  j[6][4] = x10 * std::pow(x5, -2.0) * std::cos(x10 / x5 + x8);
  j[6][5] = -std::cos(x4);
  j[6][6] = std::pow(x2, 3.0);
  j[6][7] = -std::cos(x10 / x5 + x8);
  j[6][8] = 0.0;
  j[6][9] = -1.0 / x5 * std::cos(x10 / x5 + x8);
  j[7][0] = 2.0 * x5 * (x1 - 2.0 * x6);
  j[7][1] = -x7 * std::exp(x2 * x7 + x10);
  j[7][2] = -2.0 * std::cos(-x9 + x3);
  j[7][3] = 1.5;
  j[7][4] = std::pow(x1 - 2.0 * x6, 2.0);
  j[7][5] = -4.0 * x5 * (x1 - 2.0 * x6);
  j[7][6] = -x2 * std::exp(x2 * x7 + x10);
  j[7][7] = 0.0;
  j[7][8] = 2.0 * std::cos(-x9 + x3);
  j[7][9] = -std::exp(x2 * x7 + x10);
  j[8][0] = -3.0;
  j[8][1] = -2.0 * x8 * x10 * x7;
  j[8][2] = 0.0;
  j[8][3] = std::exp(x5 + x4);
  j[8][4] = std::exp(x5 + x4);
  j[8][5] = -7.0 * std::pow(x6, -2.0);
  j[8][6] = -2.0 * x2 * x8 * x10;
  j[8][7] = -2.0 * x2 * x10 * x7;
  j[8][8] = 3.0;
  j[8][9] = -2.0 * x2 * x8 * x7;
  j[9][0] = x10;
  j[9][1] = x9;
  j[9][2] = -x8;
  j[9][3] = std::cos(x4 + x5 + x6) * x7;
  j[9][4] = std::cos(x4 + x5 + x6) * x7;
  j[9][5] = std::cos(x4 + x5 + x6) * x7;
  j[9][6] = std::sin(x4 + x5 + x6);
  j[9][7] = -x3;
  j[9][8] = x2;
  j[9][9] = x1;

  return j;
}
