#include "linear_solver.hpp"

#include <cmath>
#include <stdexcept>
#include <utility>

namespace {
void validate_square_matrix(const Matrix& a) {
    if (a.empty()) {
        throw std::invalid_argument("Matrix must not be empty.");
    }

    const std::size_t n = a.size();
    for (const Vector& row : a) {
        if (row.size() != n) {
            throw std::invalid_argument("Matrix must be square.");
        }
    }
}

bool is_square_size_compatible(const LuFactorization& lu) {
    return lu.l.size() == lu.u.size()
        && lu.l.size() == lu.permutation.size();
}
}

long long OperationStats::total() const {
    return additions + multiplications + divisions;
}

void OperationStats::reset() {
    additions = 0;
    multiplications = 0;
    divisions = 0;
}

LuFactorization decompose_lu(Matrix a) {
    validate_square_matrix(a);

    LuFactorization result;
    const std::size_t n = a.size();
    result.l.assign(n, Vector(n, 0.0));
    result.u = std::move(a);
    result.permutation.resize(n);

    for (std::size_t i = 0; i < n; ++i) {
        result.l[i][i] = 1.0;
        result.permutation[i] = i;
    }

    for (std::size_t k = 0; k < n; ++k) {
        // Частичный выбор главного элемента: ищем максимум в текущем столбце.
        std::size_t pivot = k;
        double max_value = std::abs(result.u[k][k]);

        for (std::size_t i = k + 1; i < n; ++i) {
            const double value = std::abs(result.u[i][k]);
            if (value > max_value) {
                max_value = value;
                pivot = i;
            }
        }

        if (std::abs(result.u[pivot][k]) < 1e-10) {
            throw std::runtime_error("System cannot be solved: singular matrix.");
        }

        if (pivot != k) {
            std::swap(result.u[pivot], result.u[k]);
            std::swap(result.permutation[pivot], result.permutation[k]);

            for (std::size_t j = 0; j < k; ++j) {
                std::swap(result.l[pivot][j], result.l[k][j]);
            }
        }

        for (std::size_t i = k + 1; i < n; ++i) {
            // Зануляем элементы ниже диагонали; множители сохраняются в L.
            result.l[i][k] = result.u[i][k] / result.u[k][k];
            ++result.stats.divisions;

            for (std::size_t j = k; j < n; ++j) {
                result.u[i][j] -= result.l[i][k] * result.u[k][j];
                ++result.stats.multiplications;
                ++result.stats.additions;
            }
        }
    }

    return result;
}

LinearSolveResult solve_with_lu(const LuFactorization& lu, const Vector& b) {
    if (!is_square_size_compatible(lu)) {
        throw std::invalid_argument("Invalid LU factorization.");
    }

    if (b.size() != lu.l.size()) {
        throw std::invalid_argument("Right-hand side has incompatible size.");
    }

    LinearSolveResult result;
    const std::size_t n = lu.l.size();

    Vector pb(n);
    for (std::size_t i = 0; i < n; ++i) {
        pb[i] = b[lu.permutation[i]];
    }

    // Прямой ход: решаем Ly = Pb.
    Vector y(n);
    for (std::size_t i = 0; i < n; ++i) {
        double sum = pb[i];
        for (std::size_t k = 0; k < i; ++k) {
            sum -= lu.l[i][k] * y[k];
            ++result.stats.multiplications;
            ++result.stats.additions;
        }
        y[i] = sum;
    }

    // Обратный ход: решаем Ux = y.
    result.x.assign(n, 0.0);
    for (std::size_t reverse_i = 0; reverse_i < n; ++reverse_i) {
        const std::size_t i = n - 1 - reverse_i;
        double sum = y[i];

        for (std::size_t k = i + 1; k < n; ++k) {
            sum -= lu.u[i][k] * result.x[k];
            ++result.stats.multiplications;
            ++result.stats.additions;
        }

        if (std::abs(lu.u[i][i]) < 1e-10) {
            throw std::runtime_error("System cannot be solved: singular matrix.");
        }

        result.x[i] = sum / lu.u[i][i];
        ++result.stats.divisions;
    }

    return result;
}

LinearSolveResult solve_linear_system_lu(Matrix a, Vector b) {
    LuFactorization lu = decompose_lu(std::move(a));
    LinearSolveResult result = solve_with_lu(lu, b);
    result.stats.additions += lu.stats.additions;
    result.stats.multiplications += lu.stats.multiplications;
    result.stats.divisions += lu.stats.divisions;
    return result;
}
