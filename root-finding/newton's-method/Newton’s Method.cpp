#include <cmath>
#include <functional>
#include <iostream>
#include <stdexcept>

double newtonMethod(
    const std::function<double(double)>& f,
    const std::function<double(double)>& derivative,
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

    std::cout << "Iteration    Approximation          f(x)\n";
    std::cout << "----------------------------------------\n";

    for (int iteration = 1; iteration <= maxIterations; iteration++) {
        const double fx = f(xCurrent);
        const double dfx = derivative(xCurrent);

        if (std::abs(dfx) < 1e-12) {
            throw std::runtime_error(
                "The derivative is too close to zero."
            );
        }

        const double xNext = xCurrent - fx / dfx;
        const double fxNext = f(xNext);

        std::cout
            << iteration << "            "
            << xNext << "            "
            << fxNext << '\n';

        if (std::abs(fxNext) < tolerance) {
            return xNext;
        }

        xCurrent = xNext;
    }

    throw std::runtime_error(
        "Newton's Method did not converge "
        "within the maximum number of iterations."
    );
}

int main() {
    const auto f = [](double x) {
        return x * x - 2.0;
    };

    const auto derivative = [](double x) {
        return 2.0 * x;
    };

    try {
        const double root = newtonMethod(
            f,
            derivative,
            1.0,
            1e-8,
            100
        );

        std::cout << "\nApproximate root: " << root << '\n';
    }
    catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << '\n';
    }

    return 0;
}