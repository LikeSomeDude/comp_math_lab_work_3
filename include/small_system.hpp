#pragma once

#include "types.hpp"

// система варианта 6 из методички
// sin(y + 0.5) - x = 1,
// y + cos(x - 2) = 0
Vector small_system_values(const Vector &x); // возвращает F(x)
Matrix small_system_jacobian(const Vector &x); // возвращает якобиан
Vector small_system_initial_guess(); // начальное приближение
