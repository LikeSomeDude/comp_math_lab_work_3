#pragma once

#include <vector>

using Vector = std::vector<double>;
using Matrix = std::vector<std::vector<double>>;

// счетчик арифметических операций при решении линейных систем
struct OperationStats {
  long long additions = 0;       // сложение
  long long multiplications = 0; // умножение
  long long divisions = 0;       // деление

  long long total() const;
  void reset();
};
