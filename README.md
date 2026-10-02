# Numerical-Methods-CPP
Modern C++ implementations of classical numerical methods with mathematical derivations, convergence analysis, and numerical examples.

## What's Covered

Each method includes not only a C++ implementation but also the mathematical concepts behind the algorithm.

### Root-Finding Methods
| Method | Key Concepts |
|--------|---------------------|
| Bisection Method | Intermediate Value Theorem (IVT), error bound, linear convergence |
| False Position Method | Linear interpolation, secant line, linear convergence |
| Secant Method | Taylor expansion, order of convergence, superlinear convergence |
| Fixed Point Method | Limits, continuity, differentiability, local convergence |
| Newton's Method | Taylor's theorem with Lagrange remainder, quadratic convergence |

### Solving Linear Systems
| Method | Key Concepts |
|---|---|
| Gaussian Elimination | Forward elimination, back substitution, partial pivoting, numerical stability |
| LU Decomposition | Connection to Gaussian Elimination, PA=LU, LU solving |
| Jacobi Method | Iteration formula, Triangle inequality proof, Spectral radius, Convergence condition, Strict diagonal dominance |

## Directory Structure
```text
Numerical-Methods-CPP
│
├── root-finding/
│   ├── bisection-method/
│   ├── falsePosition-method/
│   ├── secant-method/
│   ├── fixedPoint-method/
│   └── newton's-method/
│
└── solving-linear-system/
    ├── Solving_Linear_Systems.pdf
    ├── gaussian_elimination.cpp
    ├── lu_decomposition.cpp
    └── jacobi_method.cpp
```

## License
This project is released under the MIT License.
