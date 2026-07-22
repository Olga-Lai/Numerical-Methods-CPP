#include <cmath>
#include <functional>
#include <iostream>
#include <stdexcept>

double bisection(
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

    std::cout << "Iteration    Midpoint       f(midpoint)\n";
    std::cout << "----------------------------------------\n";


    for (int iteration = 0; iteration < maxIterations; iteration++) {
        const double midpoint = a + (b-a) / 2.0;
        const double fm = f(midpoint);

        std::cout << iteration + 1 << "            "
                  << midpoint << "            "
                  << fm << '\n';

        if (std::abs(fm) < tolerance ||
            (b-a) < tolerance) {
            return midpoint;
            }

        if ((fa > 0.0) != (fm > 0.0)) {
            b = midpoint;
        } else {
            a = midpoint;
            fa = fm;
        }
    }

    return a + (b-a) / 2.0;
}

int main() {
    const auto f = [](double x) {
        return x * x - 2.0;
    };

    try {
        const double root = bisection(f, 1.0, 2.0, 1e-8, 100);

        std::cout << "Approximate root: " << root << '\n';
        std::cout << "f(root): " << f(root) << '\n';
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
    }
}