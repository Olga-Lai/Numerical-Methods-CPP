#include <cmath>
#include <functional>
#include <iostream>
#include <stdexcept>

double falsePosition(
    const std::function<double(double)>& f,
    double a,
    double b,
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

    if (a >= b) {
        throw std::invalid_argument(
            "The left endpoint must be smaller than the right endpoint."
        );
    }

    double fa = f(a);
    double fb = f(b);

    if ((fa > 0.0) == (fb > 0.0)) {
        throw std::invalid_argument(
            "The interval must contain a sign change.");
    }

    std::cout << "Iteration    xr                    f(xr)\n";
    std::cout << "----------------------------------------\n";

    double xr = a;
    for (int iteration = 0; iteration < maxIterations; iteration++) {
        double denominator = fb - fa;
        if (std::abs(denominator) < 1e-15) {
            throw std::runtime_error(
                "The denominator is too small."
            );
        }

        //xr = (a * fb - b * fa) / denominator;
        xr = a - fa * (b-a) / denominator;
        double fxr = f(xr);

        std::cout << iteration + 1 << "            "
                  << xr << "            "
                  << fxr << '\n';

        if (std::abs(fxr) < tolerance) {
            return xr;
        }

        if ((fa > 0.0) != (fxr > 0.0)) {
            b = xr;
            fb = fxr;
        } else {
            a = xr;
            fa = fxr;
        }
    }

    throw std::runtime_error(
        "The method did not converge within the iteration limit."
    );
}

int main() {
    const auto f = [](double x) {
        return x * x - 2.0;
    };

    try {
        const double root = falsePosition(f, 1.0, 2.0, 1e-8, 100);

        std::cout << "Approximate root: " << root << '\n';
        std::cout << "f(root): " << f(root) << '\n';
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
    }
}