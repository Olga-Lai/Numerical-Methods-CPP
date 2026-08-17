#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

std::vector<double> gaussianElimination(
    std::vector<std::vector<double>> A,
    std::vector<double> b
) {
    const int n = A.size();
    const double epsilon = 1e-12;

    if (n == 0 || b.size() != n) {
        throw std::invalid_argument(
            "Invalid matrix or vector size."
        );
    }

    for (const auto& row : A) {
        if (row.size() != n) {
            throw std::invalid_argument(
                "Matrix must be square."
            );
        }
    }

    // Forward elimination
    for (int k = 0; k < n-1; k++) {

        // 1. Select pivot row
        int pivotRow = k;

        for (int i = k + 1; i < n; i++) {
            if (std::abs(A[i][k]) >
                std::abs(A[pivotRow][k])) {
                pivotRow = i;
            }
        }

        // 2. Check pivot
        if (std::abs(A[pivotRow][k]) < epsilon) {
            throw std::runtime_error(
                "Matrix is singular or nearly singular."
            );
        }

        // 3. Swap rows
        if (pivotRow != k) {
            std::swap(A[k], A[pivotRow]);
            std::swap(b[k], b[pivotRow]);
        }

        // 4. Eliminate entries below pivot
        for (int i = k + 1; i < n; i++) {

            const double multiplier =
                A[i][k] / A[k][k];

            for (int j = k; j < n; j++) {
                A[i][j] -= multiplier * A[k][j];
            }

            b[i] -= multiplier * b[k];
        }
    }

    // Check final pivot
    if (std::abs(A[n - 1][n - 1]) < epsilon) {
        throw std::runtime_error(
            "Matrix is singular or nearly singular."
        );
    }

    // Back substitution
    std::vector<double> x(n);

    for (int i = n - 1; i >= 0; i--) {
        double sum = b[i];

        for (int j = i + 1; j < n; j++) {
            sum -= A[i][j] * x[j];
        }

        x[i] = sum / A[i][i];
    }

    return x;
}

int main() {
    std::vector<std::vector<double>> A = {
        { 2.0,  1.0, -1.0},
        {-3.0, -1.0,  2.0},
        {-2.0,  1.0,  2.0}
    };

    std::vector<double> b = {
         8.0,
       -11.0,
        -3.0
    };

    try {
        const std::vector<double> x =
            gaussianElimination(A, b);

        std::cout << "Solution:\n";

        for (int i = 0; i < x.size(); ++i) {
            std::cout
                << "x" << i + 1
                << " = " << x[i]
                << '\n';
        }

    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
    }

    return 0;
}