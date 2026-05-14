#include "newton.hpp"
#include "nonlinear_system.hpp"

#include <chrono>
#include <exception>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace {
void print_vector(const Vector& x) {
    std::cout << std::fixed << std::setprecision(12);
    for (std::size_t i = 0; i < x.size(); ++i) {
        std::cout << "x" << i + 1 << " = " << x[i] << '\n';
    }
}

struct TimedResult {
    NewtonResult result;
    long long microseconds = 0;
};

struct ExperimentRow {
    std::size_t m = 1;
    std::size_t k = 1;
    TimedResult timed;
};

TimedResult run_newton(const Vector& x0, const NewtonOptions& options) {
    const auto started_at = std::chrono::steady_clock::now();
    const NewtonResult result = solve_newton(x0, options);
    const auto finished_at = std::chrono::steady_clock::now();

    const auto microseconds = std::chrono::duration_cast<std::chrono::microseconds>(
        finished_at - started_at
    ).count();

    return {result, microseconds};
}

void print_case(const std::string& name, const Vector& x0, const NewtonOptions& options) {
    const TimedResult timed = run_newton(x0, options);
    const NewtonResult& result = timed.result;

    std::cout << "\n" << name << "\n";
    std::cout << "converged: " << std::boolalpha << result.converged << '\n';
    std::cout << "iterations: " << result.iterations << '\n';
    std::cout << "jacobian recalculations: " << result.jacobian_recalculations << '\n';
    std::cout << "||F(x)||_inf: " << result.residual_norm << '\n';
    std::cout << "linear solve operations: " << result.linear_solve_stats.total() << '\n';
    std::cout << "time, us: " << timed.microseconds << '\n';
    print_vector(result.x);
}

NewtonOptions options_for(std::size_t m, std::size_t k) {
    NewtonOptions options;
    options.jacobian_recalculation_period = m;
    options.full_jacobian_iterations = k;
    return options;
}

void print_table_header(const std::string& title) {
    std::cout << "\n" << title << "\n";
    std::cout
        << std::setw(10) << "m"
        << std::setw(10) << "k"
        << std::setw(12) << "conv"
        << std::setw(12) << "iters"
        << std::setw(12) << "J count"
        << std::setw(16) << "||F(x)||"
        << std::setw(16) << "operations"
        << std::setw(12) << "time us"
        << '\n';
}

void print_table_row(std::size_t m, std::size_t k, const TimedResult& timed) {
    const NewtonResult& result = timed.result;
    const std::string m_text = m == std::numeric_limits<std::size_t>::max() ? "inf" : std::to_string(m);
    const std::string k_text = k == std::numeric_limits<std::size_t>::max() ? "inf" : std::to_string(k);

    std::cout
        << std::setw(10) << m_text
        << std::setw(10) << k_text
        << std::setw(12) << std::boolalpha << result.converged
        << std::setw(12) << result.iterations
        << std::setw(12) << result.jacobian_recalculations
        << std::setw(16) << std::scientific << std::setprecision(3) << result.residual_norm
        << std::setw(16) << std::fixed << result.linear_solve_stats.total()
        << std::setw(12) << timed.microseconds
        << '\n';
}

void print_best_results(const std::vector<ExperimentRow>& rows) {
    const ExperimentRow* best_time = nullptr;
    const ExperimentRow* best_operations = nullptr;

    for (const ExperimentRow& row : rows) {
        if (!row.timed.result.converged) {
            continue;
        }

        if (best_time == nullptr || row.timed.microseconds < best_time->timed.microseconds) {
            best_time = &row;
        }

        if (best_operations == nullptr
            || row.timed.result.linear_solve_stats.total()
                < best_operations->timed.result.linear_solve_stats.total()) {
            best_operations = &row;
        }
    }

    if (best_time == nullptr || best_operations == nullptr) {
        std::cout << "No converged runs in this table.\n";
        return;
    }

    std::cout << "Best by time:\n";
    print_table_row(best_time->m, best_time->k, best_time->timed);
    std::cout << "Best by operations:\n";
    print_table_row(best_operations->m, best_operations->k, best_operations->timed);
}

void run_k_experiments(const Vector& x0) {
    const std::size_t infinity = std::numeric_limits<std::size_t>::max();
    std::vector<ExperimentRow> rows;

    print_table_header("Experiments with k, then modified Newton");

    // Пункт c): первые k итераций пересчитываем Якоби, потом фиксируем ее.
    for (std::size_t k = 1; k <= 10; ++k) {
        rows.push_back({infinity, k, run_newton(x0, options_for(infinity, k))});
        print_table_row(rows.back().m, rows.back().k, rows.back().timed);
    }

    print_best_results(rows);
}

void run_m_experiments(const Vector& x0) {
    const std::size_t infinity = std::numeric_limits<std::size_t>::max();
    const std::vector<std::size_t> periods = {1, 2, 3, 4, 5, 7, 10, infinity};
    std::vector<ExperimentRow> rows;

    print_table_header("Experiments with cyclic m");
    // Пункт d): после начального шага пересчитываем Якоби каждые m итераций.
    for (std::size_t m : periods) {
        rows.push_back({m, 1, run_newton(x0, options_for(m, 1))});
        print_table_row(rows.back().m, rows.back().k, rows.back().timed);
    }

    print_best_results(rows);
}

void run_changed_initial_guess_experiments() {
    const std::size_t infinity = std::numeric_limits<std::size_t>::max();
    const Vector x0 = changed_initial_guess();
    std::vector<ExperimentRow> rows;

    print_table_header("Experiments with x5(0) = -0.2");
    // Пункт f): проверяем поведение при другом x5 и разных k.
    for (std::size_t k = 1; k <= 10; ++k) {
        rows.push_back({infinity, k, run_newton(x0, options_for(infinity, k))});
        print_table_row(rows.back().m, rows.back().k, rows.back().timed);
    }

    print_best_results(rows);
}
}

int main() {
    try {
        const std::size_t infinity = std::numeric_limits<std::size_t>::max();

        print_case("Full Newton", initial_guess(), options_for(1, infinity));
        print_case("Modified Newton", initial_guess(), options_for(infinity, 1));
        print_case("Cyclic Newton, m = 3", initial_guess(), options_for(3, 1));

        run_k_experiments(initial_guess());
        run_m_experiments(initial_guess());
        run_changed_initial_guess_experiments();
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << '\n';
        return 1;
    }

    return 0;
}
