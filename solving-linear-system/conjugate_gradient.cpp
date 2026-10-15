
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

using namespace std;

struct CGResult {
    vector<double> x;
    int iterations;
    double residualNorm;
    bool converged;
};

double dot(const vector<double>& a, const vector<double>& b) {
    if (a.size() != b.size()) {
        throw invalid_argument("Vector dimensions do not match.");
    }

    double result = 0.0;

    for (size_t i = 0; i < a.size(); ++i) {
        result += a[i] * b[i];
    }

    return result;
}

vector<double> matVec(
    const vector<vector<double>>& A,
    const vector<double>& x
) {
    const size_t n = A.size();
    vector<double> result(n, 0.0);

    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            result[i] += A[i][j] * x[j];
        }
    }

    return result;
}

CGResult conjugateGradient(
    const vector<vector<double>>& A,
    const vector<double>& b,
    const vector<double>& x0,
    int maxIterations = 100,
    double tolerance = 1e-8
) {
    const size_t n = A.size();

    if (n == 0 || b.size() != n || x0.size() != n) {
        throw invalid_argument("Invalid matrix or vector dimensions.");
    }

    for (const auto& row : A) {
        if (row.size() != n) {
            throw invalid_argument("Matrix must be square.");
        }
    }

    if (maxIterations < 0 || tolerance <= 0.0 ||
        !isfinite(tolerance)) {
        throw invalid_argument("Invalid iteration limit or tolerance.");
    }

    // Initial guess and residual
    vector<double> x = x0;
    vector<double> Ax = matVec(A, x);
    vector<double> r(n);

    for (size_t i = 0; i < n; ++i) {
        r[i] = b[i] - Ax[i];
    }

    double rTr = dot(r, r);
    double residualNorm = sqrt(rTr);

    if (residualNorm <= tolerance) {
        return {x, 0, residualNorm, true};
    }

    // Initial search direction
    vector<double> p = r;

    for (int k = 0; k < maxIterations; ++k) {
        // Compute A * p
        vector<double> Ap = matVec(A, p);

        // Compute alpha_k
        double pTAp = dot(p, Ap);

        if (!(pTAp > 0.0) || !isfinite(pTAp)) {
            throw runtime_error(
                "CG breakdown: A may not be SPD."
            );
        }

        double alpha = rTr / pTAp;

        // Update x and r
        for (size_t i = 0; i < n; ++i) {
            x[i] += alpha * p[i];
            r[i] -= alpha * Ap[i];
        }

        // Check convergence
        double rTrNew = dot(r, r);
        residualNorm = sqrt(rTrNew);

        if (!isfinite(residualNorm)) {
            throw runtime_error(
                "CG breakdown: non-finite residual."
            );
        }

        if (residualNorm <= tolerance) {
            return {x, k + 1, residualNorm, true};
        }

        // Compute beta_k
        double beta = rTrNew / rTr;

        // Update search direction
        for (size_t i = 0; i < n; ++i) {
            p[i] = r[i] + beta * p[i];
        }

        rTr = rTrNew;
    }

    return {x, maxIterations, residualNorm, false};
}



int main() {
    vector<vector<double>> A = {
        {4.0, 1.0, 0.0},
        {1.0, 3.0, 1.0},
        {0.0, 1.0, 2.0}
    };

    vector<double> b = {1.0, 2.0, 3.0};
    vector<double> x0 = {0.0, 0.0, 0.0};

    CGResult result = conjugateGradient(A, b, x0);

    cout << "Solution:\n";

    for (double value : result.x) {
        cout << value << '\n';
    }

    cout << "Iterations: " << result.iterations << '\n';
    cout << "Residual norm: " << result.residualNorm << '\n';
    cout << "Converged: " << (result.converged ? "Yes" : "No") << '\n';

    return 0;
}

