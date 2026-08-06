#include <cmath>
#include <functional>
#include <iostream>
#include <stdexcept>

double secant(
    const std::function<double(double)>& f,
    double x0,
    double x1,
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

    double f0 = f(x0);
    double f1 = f(x1);

    std::cout << "Iteration    Approximation          f(x)\n";
    std::cout << "----------------------------------------\n";

    for (int iteration = 0; iteration < maxIterations; iteration++) {
        const double denominator = f1 - f0;
        if (std::abs(denominator) < 1e-15) {
            throw std::runtime_error(
                "The denominator is too small."
            );
        }

        const double x2 = x1 - f1 * (x1 - x0) / denominator;
        const double f2 = f(x2);

        std::cout << iteration + 1 << "            "
                  << x2 << "            "
                  << f2 << '\n';

        if (std::abs(f2) < tolerance) {
            return x2;
        }

        x0 = x1;
        f0 = f1;

        x1 = x2;
        f1 = f2;
    }

    throw std::runtime_error(
        "Secant Method did not converge "
        "within the maximum number of iterations."
    );
}

int main() {
    const auto f = [](double x) {
        return x * x - 2.0;
    };

    try {
        const double root = secant(f, 1.0, 2.0, 1e-8, 100);

        std::cout << "Approximate root: " << root << '\n';
        std::cout << "f(root): " << f(root) << '\n';
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
    }
}