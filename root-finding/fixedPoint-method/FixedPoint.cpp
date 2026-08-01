#include <cmath>
#include <functional>
#include <iostream>
#include <stdexcept>

double fixedPoint(
    const std::function<double(double)>& g,
    double initialGuess,
    const double tolerance,
    const int maxIterations
    ) {

    if (tolerance <= 0.0) {
        throw std::invalid_argument(
            "Tolerance must be positive."
        );
    }

    if (maxIterations <= 0) {
        throw std::invalid_argument(
            "Maximum iterations must be positive."
        );
    }

    double xCurrent = initialGuess;

    std::cout << "Iteration    Approximation          g(x)\n";
    std::cout << "----------------------------------------\n";

    for (int iteration = 0; iteration < maxIterations; iteration++) {
        double xNext = g(xCurrent);

        if (!std::isfinite(xNext)) {
            throw std::runtime_error(
                "The iteration produced a non-finite value."
            );
        }

        double error = std::abs(xNext - xCurrent);

        std::cout << iteration + 1 << "            "
                  << xCurrent << "            "
                  << xNext << '\n';

        if (error < tolerance) {
            return xNext;
        }

        xCurrent = xNext;
    }

    throw std::runtime_error(
        "The method did not converge within the iteration limit."
    );
}

int main() {
    const auto f = [](double x) {
        return x * x - x - 1.0;
    };

    const auto g = [](double x) {
        return std::sqrt(x + 1.0);
    };

    try {
        const double root = fixedPoint(g, 1.0, 1e-8, 100);

        std::cout << "Approximate root: " << root << '\n';
        std::cout << "f(root): " << f(root) << '\n';
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
    }
}