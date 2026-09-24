#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <numeric>

struct LUResult {
    std::vector<std::vector<double>> LU;
    std::vector<int> P;
};

LUResult luDecomposition(std::vector<std::vector<double>> A) {
    const int n = A.size();
    const double epsilon = 1e-12;

    if (n == 0) {
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

    std::vector<int> P(n);
    std::iota(P.begin(), P.end(), 0);

    for (int k = 0; k < n; k++) {

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
                "Matrix is singular or nearly singular.");
        }

        // 3. Swap rows
        if (pivotRow != k) {
            std::swap(A[k], A[pivotRow]);
            std::swap(P[k], P[pivotRow]);
        }

        // 4. Store multipliers and eliminate
        for (int i = k + 1; i < n; i++) {

            A[i][k] /= A[k][k]; // Store L multiplier

            for (int j = k + 1; j < n; j++) {
                A[i][j] -= A[i][k] * A[k][j];
            }
        }
    }

    return {A, P};
}

void printLU(const LUResult& lu) {
    const int n = lu.LU.size();

    std::cout << "L:\n";
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == j)
                std::cout << 1.0 << "\t";
            else if (i > j)
                std::cout << lu.LU[i][j] << "\t";
            else
                std::cout << 0.0 << "\t";
        }
        std::cout << '\n';
    }

    std::cout << "\nU:\n";
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i <= j)
                std::cout << lu.LU[i][j] << "\t";
            else
                std::cout << 0.0 << "\t";
        }
        std::cout << '\n';
    }

    std::cout << '\n';
}

std::vector<double> luSolve(const LUResult& lu, std::vector<double> b) {
    const int n = lu.LU.size();

    if (b.size() != n) {
        throw std::invalid_argument(
            "Vector size matches incorrectly.");
    }

    std::vector<double> Pb(n);
    for (int i = 0; i < n; i++) {
        Pb[i] = b[lu.P[i]];
    }

    // forward substitution
    std::vector<double> y(n);
    for (int i = 0; i < n; i++) {
        double sum = Pb[i];
        for (int j = 0; j < i; j++) {
            sum -= lu.LU[i][j] * y[j];
        }
        y[i] = sum;
    }

    // backward substitution
    std::vector<double> x(n);
    for (int i = n - 1; i >= 0; i--) {
        double sum = y[i];
        for (int j = i + 1; j < n; j++) {
            sum -= lu.LU[i][j] * x[j];
        }
        x[i] = sum / lu.LU[i][i];
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
        LUResult lu = luDecomposition(A);

        printLU(lu);

        const std::vector<double> x = luSolve(lu, b);

        std::cout << "Solution:\n";
        for (size_t i = 0; i < x.size(); ++i) {
            std::cout << "x" << i + 1 << " = " << x[i] << '\n';
        }

    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
    }

    return 0;
}