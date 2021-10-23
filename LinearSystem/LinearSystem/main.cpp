#include "Math.h"

#include <iostream>
#include <iomanip>

#define DOCTEST_CONFIG_IMPLEMENT
#include "doctest.h"

#include "LinearSystemSolver.h"

int main()
{
    doctest::Context context;
    context.run();

    //using namespace NumericalAnalysis;

    using namespace NumericalAnalysis;

    std::cout << "========================================================================" << std::endl
        << "Gaussian elimination method:" << std::endl;
    {
        std::vector<std::vector<double>> A {
            {5.0, 2.0, 1.0, 0.0},
            {1.0, 3.0, 2.0, 8.0},
            {4.0, -6.0, 1.0, 0.0},
            {5.0, 0.0, 3.0, 2.0},
        };
        std::vector<double> b{ 14.0, 65.0, -3.0, 32.0 };
        auto result = LinearSystemSolver::GaussianElimination(A, b);
        std::cout << "(" << result.at(0);
        for (size_t i = 1; i < result.size(); ++i) {
            std::cout << ", " << result[i];
        }
        std::cout << ")" << std::endl;
    }
    std::cout << "========================================================================" << std::endl
        << "Tridiagonal matrix method:" << std::endl;
    {
        std::vector<std::vector<double>> A{
            {2.0, 4.0, 0.0},
            {4.0, 1.0, 5.0},
            {0.0, 5.0, 2.0},
        };
        std::vector<double> b{ 20.0, 37.0, 30.0 };
        auto result = LinearSystemSolver::TridiagonalMatrix(A, b);
        std::cout << "(" << result.at(0);
        for (size_t i = 1; i < result.size(); ++i) {
            std::cout << ", " << result[i];
        }
        std::cout << ")" << std::endl;
    }
    std::cout << "========================================================================" << std::endl
        << "Jacobi method:" << std::endl;
    {
        std::vector<std::vector<double>> A{
            {6.0, 3.0, 1.0, 0.0},
            {3.0, 5.0, 0.0, 2.0},
            {1.0, 0.0, 3.0, 1.0},
            {0.0, 2.0, 1.0, 5.0},
        };
        std::vector<double> b{ 25.0, 31.0, 19.0, 35.0 };
        auto result = LinearSystemSolver::Jacobi(A, b, 0.001);
        for (size_t i = 0; i < result.size(); ++i) {
            std::cout << "Iteration " << i << ": (" << result.at(i).at(0);
            for (size_t j = 1; j < result.at(i).size(); ++j) {
                std::cout << ", " << result.at(i).at(j);
            }
            std::cout << ")" << std::endl;
        }
    }
    std::cout << "========================================================================" << std::endl
        << "Upper relaxation method:" << std::endl;
    {
        std::vector<std::vector<double>> A{
            {6.0, 3.0, 1.0, 0.0},
            {3.0, 5.0, 0.0, 2.0},
            {1.0, 0.0, 3.0, 1.0},
            {0.0, 2.0, 1.0, 5.0},
        };
        std::vector<double> b{ 25.0, 31.0, 19.0, 35.0 };
        auto result = LinearSystemSolver::UpperRelaxation(A, b, 0.001);
        for (size_t i = 0; i < result.size(); ++i) {
            std::cout << "Iteration " << i << ": (" << result.at(i).at(0);
            for (size_t j = 1; j < result.at(i).size(); ++j) {
                std::cout << ", " << result.at(i).at(j);
            }
            std::cout << ")" << std::endl;
        }
    }
}