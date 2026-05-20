#include "newton.hpp"
#include "nonlinear_system.hpp"
#include "scalar_newton.hpp"
#include "small_system.hpp"

#include <chrono>
#include <exception>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace {
void print_section_title(const std::string& title) {
    std::cout << "\n==== " << title << " ====\n";
}

void print_scalar_newton_result() {
    const ScalarNewtonResult result = solve_scalar_newton({});

    print_section_title("Scalar Newton method");
    std::cout << "Equation: (x + 3)cos(x) = 1\n";
    std::cout << "localized segment: ["
              << result.segment_left << ", " << result.segment_right << "]\n";
    std::cout << "converged: " << std::boolalpha << result.converged << '\n';
    std::cout << "iterations: " << result.iterations << '\n';
    std::cout << "root: " << std::fixed << std::setprecision(8) << result.root << '\n';
    std::cout << "|f(x)|: " << std::scientific << result.residual << '\n';
    std::cout << std::defaultfloat;
}

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

TimedResult run_newton(
    const Vector& x0,
    const SystemFunction& values,
    const JacobianFunction& jacobian,
    const NewtonOptions& options
) {
    const auto started_at = std::chrono::steady_clock::now();
    const NewtonResult result = solve_newton(x0, values, jacobian, options);
    const auto finished_at = std::chrono::steady_clock::now();

    const auto microseconds = std::chrono::duration_cast<std::chrono::microseconds>(
        finished_at - started_at
    ).count();

    return {result, microseconds};
}

void print_case(const std::string& name, const Vector& x0, const NewtonOptions& options, bool print_solution) {
    const TimedResult timed = run_newton(x0, options);
    const NewtonResult& result = timed.result;

    std::cout << std::defaultfloat;
    std::cout << "\n" << name << "\n";
    std::cout << "converged: " << std::boolalpha << result.converged << '\n';
    std::cout << "iterations: " << result.iterations << '\n';
    std::cout << "jacobian recalculations: " << result.jacobian_recalculations << '\n';
    std::cout << "||F(x)||_inf: " << result.residual_norm << '\n';
    std::cout << "linear solve operations: " << result.linear_solve_stats.total() << '\n';
    std::cout << "time, us: " << timed.microseconds << '\n';

    if (print_solution) {
        print_vector(result.x);
    }
}

void print_small_system_case() {
    NewtonOptions options;
    options.jacobian_recalculation_period = 1;
    options.full_jacobian_iterations = std::numeric_limits<std::size_t>::max();

    const TimedResult timed = run_newton(
        small_system_initial_guess(),
        small_system_values,
        small_system_jacobian,
        options
    );

    const NewtonResult& result = timed.result;

    print_section_title("Two-dimensional system, variant 6");
    std::cout << "System:\n";
    std::cout << "sin(y + 0.5) - x = 1\n";
    std::cout << "y + cos(x - 2) = 0\n";
    std::cout << "converged: " << std::boolalpha << result.converged << '\n';
    std::cout << "iterations: " << result.iterations << '\n';
    std::cout << "jacobian recalculations: " << result.jacobian_recalculations << '\n';
    std::cout << "||F(x)||_inf: " << result.residual_norm << '\n';
    std::cout << "linear solve operations: " << result.linear_solve_stats.total() << '\n';
    std::cout << "time, us: " << timed.microseconds << '\n';
    std::cout << std::fixed << std::setprecision(12);
    std::cout << "x = " << result.x[0] << '\n';
    std::cout << "y = " << result.x[1] << '\n';
    std::cout << std::defaultfloat;
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

    print_table_header("Parameter k: full Newton first, modified Newton after that");

    // пункт c: первые k итераций пересчитываем Якоби, потом фиксируем ее
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

    print_table_header("Parameter m: cyclic Jacobian recalculation");
    // пункт d: после начального шага пересчитываем Якоби каждые m итераций
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

    print_table_header("Changed initial guess: x5(0) = -0.2");
    // пункт f: проверяем поведение при другом x5 и разных k
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

        print_scalar_newton_result();
        print_small_system_case();

        print_section_title("Ten-dimensional system Newton methods");
        print_case("Full Newton", initial_guess(), options_for(1, infinity), true);
        print_case("Modified Newton", initial_guess(), options_for(infinity, 1), false);
        print_case("Cyclic Newton, m = 3", initial_guess(), options_for(3, 1), false);

        run_k_experiments(initial_guess());
        run_m_experiments(initial_guess());
        run_changed_initial_guess_experiments();
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << '\n';
        return 1;
    }

    return 0;
}
