#include "scalar_newton.hpp"

#include <cmath>
#include <stdexcept>
#include <vector>

namespace {
struct Segment {
    double left = 0.0;
    double right = 0.0;
};

bool is_nonzero_root(double x, double eps) {
    return std::abs(x) > eps;
}

bool is_new_root(const std::vector<ScalarNewtonResult>& roots, double root, double eps) {
    for (const ScalarNewtonResult& result : roots) {
        if (std::abs(result.root - root) < 10.0 * eps) {
            return false;
        }
    }

    return true;
}

ScalarNewtonResult solve_on_segment(const Segment& segment, const ScalarNewtonOptions& options) {
    ScalarNewtonResult result;
    result.segment_left = segment.left;
    result.segment_right = segment.right;
    result.root = (segment.left + segment.right) / 2.0;

    for (std::size_t iteration = 0; iteration < options.max_iterations; ++iteration) {
        const double value = scalar_function(result.root);
        const double derivative = scalar_function_derivative(result.root);

        if (std::abs(value) < options.eps) {
            result.residual = std::abs(value);
            result.iterations = iteration;
            result.converged = true;
            return result;
        }

        if (std::abs(derivative) < 1e-12) {
            break;
        }

        const double next = result.root - value / derivative;
        result.iterations = iteration + 1;

        if (std::abs(next - result.root) < options.eps) {
            result.root = next;
            result.residual = std::abs(scalar_function(result.root));
            result.converged = result.residual < options.eps;
            return result;
        }

        result.root = next;
    }

    result.residual = std::abs(scalar_function(result.root));
    return result;
}

std::vector<Segment> localize_roots(const ScalarNewtonOptions& options) {
    if (options.search_step <= 0.0 || options.search_left >= options.search_right) {
        throw std::invalid_argument("Invalid scalar root localization interval.");
    }

    std::vector<Segment> segments;
    double left = options.search_left;
    double f_left = scalar_function(left);

    for (double right = left + options.search_step; right <= options.search_right; right += options.search_step) {
        const double f_right = scalar_function(right);

        if (f_left == 0.0 || f_left * f_right < 0.0) {
            segments.push_back({left, right});
        }

        left = right;
        f_left = f_right;
    }

    return segments;
}
}

double scalar_function(double x) {
    return (x + 3.0) * std::cos(x) - 1.0;
}

double scalar_function_derivative(double x) {
    return std::cos(x) - (x + 3.0) * std::sin(x);
}

ScalarNewtonResult solve_scalar_newton(const ScalarNewtonOptions& options) {
    const std::vector<Segment> segments = localize_roots(options);
    std::vector<ScalarNewtonResult> roots;

    for (const Segment& segment : segments) {
        ScalarNewtonResult result = solve_on_segment(segment, options);

        if (result.converged
            && is_nonzero_root(result.root, options.eps)
            && is_new_root(roots, result.root, options.eps)) {
            roots.push_back(result);
        }
    }

    if (roots.empty()) {
        throw std::runtime_error("Scalar Newton method did not find a nonzero root.");
    }

    ScalarNewtonResult best = roots.front();
    for (const ScalarNewtonResult& result : roots) {
        if (std::abs(result.root) < std::abs(best.root)) {
            best = result;
        }
    }

    return best;
}
