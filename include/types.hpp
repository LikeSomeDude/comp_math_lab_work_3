#pragma once

#include <vector>

using Vector = std::vector<double>;
using Matrix = std::vector<std::vector<double>>;

// Счетчик арифметических операций при решении линейных систем.
struct OperationStats {
  long long additions = 0;
  long long multiplications = 0;
  long long divisions = 0;

  long long total() const;
  void reset();
};
