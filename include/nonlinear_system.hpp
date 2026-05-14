#pragma once

#include "types.hpp"

// Считает значения функций системы: F(x).
Vector system_values(const Vector& x);

// Считает матрицу Якоби J(x) для текущего приближения.
Matrix jacobian(const Vector& x);

// Возвращает начальное приближение из условия.
Vector initial_guess();

// Возвращает начальное приближение для пункта f): x5 = -0.2.
Vector changed_initial_guess();
