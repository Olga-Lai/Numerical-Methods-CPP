#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

vector<double> gaussSeidel(
    const vector<vector<double>>& A,
    const vector<double>& b,
    vector<double> x,
    int maxIterations = 100,
    double tolerance = 1e-8
) {
    int n = A.size();

    for (int k = 0; k < maxIterations; ++k) {
        vector<double> x_old = x;

        for (int i = 0; i < n; ++i) {
            if (fabs(A[i][i]) < 1e-12) {
                throw runtime_error("Zero diagonal element.");
            }

            double sum = 0.0;

            // Use the newest values of x for j < i
            for (int j = 0; j < i; j++) {
                sum += A[i][j] * x[j];
            }

            // Use the old values of x for j > i
            for (int j = i + 1; j < n; ++j) {
                sum += A[i][j] * x_old[j];
            }

            x[i] = (b[i] - sum) / A[i][i];
        }

        // Check convergence
        double error = 0.0;
        for (int i = 0; i < n; i++) {
            error = max(error, fabs(x[i] - x_old[i]));
        }

        if (error < tolerance) {
            return x;
        }
    }

    return x;
}

int main() {
    vector<vector<double>> A = {
        {4, 1, 1},
        {2, 7, 1},
        {1, -3, 12}
    };

    vector<double> b = {7, -4, 15};

    vector<double> x = {0, 0, 0};

    vector<double> solution =
        gaussSeidel(A, b, x);

    cout << "Solution:" << endl;

    for (double value : solution) {
        cout << value << endl;
    }

    return 0;
}