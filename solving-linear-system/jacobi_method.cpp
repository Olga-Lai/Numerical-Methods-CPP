#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

std::vector<double> jacobiMethod(
    const std::vector<std::vector<double>>& A,
    const std::vector<double>& b,
    int maxIterations = 100,
    double tolerance = 1e-8)
{
    int n = A.size();

    if (n == 0 || b.size() != n)
        throw std::invalid_argument(
            "Invalid matrix dimensions.");

    std::vector<double> x(n, 0.0);
    std::vector<double> xNew(n);

    for (int iter = 0; iter < maxIterations; iter++)
    {
        for (int i = 0; i < n; i++)
        {
            if (std::abs(A[i][i]) < 1e-12)
                throw std::runtime_error(
                    "Zero diagonal element.");

            double sum = 0.0;

            for (int j = 0; j < n; j++)
            {
                if (i != j)
                    sum += A[i][j] * x[j];
            }

            xNew[i] = (b[i] - sum) / A[i][i];
        }

        // ||x^(k+1)-x^k||
        double error = 0.0;
        for (int i = 0; i < n; i++)
            error = std::max(error, std::abs(xNew[i] - x[i]));

        x = xNew;

        if (error < tolerance){
            std::cout << "Iteration: " << n << "\n";
            break;
        }

        // Residual Norm ||Ax-b||
        /*x = xNew;
        double residual = 0.0;
        for (int i = 0; i < n; i++)
        {
            double Ax = 0.0;

            for (int j = 0; j < n; j++)
                Ax += A[i][j] * x[j];

            residual = std::max(residual, std::abs(Ax - b[i]));
        }

        if (residual < tolerance) {
            std::cout << "Iteration: " << n << "\n";
            break;
        }*/
    }

    return x;
}

int main()
{
    std::vector<std::vector<double>> A = {
        {10, -1,  2},
        {-2,  3,  2},
        { 2, -1, 10}
    };

    std::vector<double> b = {6, 2, -10};

    auto solution = jacobiMethod(A, b);

    std::cout << std::fixed << std::setprecision(6);

    std::cout << "Solution (Jacobi Method):\n";
    for (double x : solution)
        std::cout << x << "\n";

    return 0;
}