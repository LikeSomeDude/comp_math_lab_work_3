#pragma once

#include "types.hpp"

// считает значения функций системы F(x)
Vector system_values(const Vector &x);

// считает матрицу Якоби J(x) для текущего приближения
Matrix jacobian(const Vector &x);

// возвращает начальное приближение из условия
Vector initial_guess();

// возвращает начальное приближение для пункта f, где x5 = -0.2
Vector changed_initial_guess();
