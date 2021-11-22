#include "Math.h"

#include <iostream>
#include <iomanip>

#define DOCTEST_CONFIG_IMPLEMENT
#include "doctest.h"

#include "LinearSystemSolver.h"

#include <cmath>

double calculateF1(double x, double y) {
    return std::pow(x, 3) - 2*x + 5 - y;
}

double calculateF2(double x, double y) {
    return -std::pow(y, -2) - 3*std::pow(y, -2) - x;
}

double getJacobi(int i, double arg) {
    switch (i) {
    case 1:
        return 3 * std::pow(arg,2) - 2;
    case 2:
        return -1;
    case 3:
        return -1;
    case 4:
        return (3*arg + 2) / std::pow(arg, 3);
    }
}

std::vector<double> multiplyMatrixVector(const std::vector<std::vector<double>>& A, const std::vector<double> b) {
    std::vector<double> res(b.size(), 0);
    for (size_t i = 0; i < A.size(); i++) {
        for (size_t j = 0; j < A.at(i).size(); j++) {
            res[i] += A[i][j] * b[j];
        }
    }
    return res;
}

void solveNewtonNonLinear(double precision) {
    double x = -2;
    double y = 0.1;
    int iter = 0;
    while (std::abs(calculateF1(x, y)) > precision || std::abs(calculateF2(x, y)) > precision) {
        std::cout << std::endl << "Iteration #" << iter++ << std::endl;
        std::cout << "x: " << x << std::endl << "y: " << y << std::endl;
        std::cout << "f1 (x, y): " << calculateF1(x, y) << std::endl;
        std::cout << "f2 (x, y): " << calculateF2(x, y) << std::endl;
        std::vector<double> b{ calculateF1(x, y), calculateF2(x, y) };
        std::vector<std::vector<double>> A{
            {getJacobi(1, x), getJacobi(2, y)},
            {getJacobi(3, x), getJacobi(4, y)}
        };
        std::cout << "Jacobi: " << std::endl;
        std::cout << A[0][0] << "\t" << A[0][1] << std::endl << A[1][0] << "\t\t" << A[1][1] << std::endl;
        auto res = NumericalAnalysis::LinearSystemSolver::GaussianElimination(A, b);
        std::cout << "z:(" << res[0] << ", " << res[1] << ")" << std::endl;
        x -= res[0];
        y -= res[1];
    }
    std::cout << std::endl << "Result: " << std::endl;
    std::cout << "x: " << x << " y: " << y << std::endl;
    std::cout << "f1 (x, y): " << calculateF1(x, y) << std::endl;
    std::cout << "f2 (x, y): " << calculateF2(x, y);
}

void printVec(const std::vector<double>& vec) {
    std::cout << "(";
    for (size_t i = 0; i < vec.size() - 1; i++) {
        std::cout << vec[i] << ", ";
    }
    std::cout << vec.back() << ")";
}

double scalarMult(const std::vector<double>& vec1, const std::vector<double>& vec2) {
    double result = 0.0;
    for (size_t i = 0; i < vec1.size(); i++) {
        result += vec1[i] * vec2[i];
    }
    return result;
}

double solvePowerIteration(const std::vector<std::vector<double>>& A, double precision, bool approach) {
    /*double a = -1;
    for (size_t i = 0; i < A.size(); i++) {
        double tryA = 0;
        for (size_t j = 0; j < A.at(i).size(); j++) {
            tryA += std::abs(A[i][j]);
        }
        a = std::max(a, tryA);
    }*/
    //auto B = A;
    /*for (auto& row : A) {
        for (auto& b : row) {
            b *= a;
        }
    }*/
    std::vector<double> x0 = { 1, 1, 1 };
    std::vector<double> x1 = multiplyMatrixVector(A, x0);
    const size_t idx = 0;
    int counter = 0;
    double m0 = x0[idx];
    std::cout << "x" << counter << ": ";
    printVec(x0);
    std::cout << std::endl;
    std::cout << "m" << counter++ << ": " << m0 << std::endl;
    double m1 = approach ? x1[idx] / x0[idx] : scalarMult(x1, x0) / scalarMult(x0, x0);
    while (std::abs(m1 - m0) > precision) {
        std::cout << "x" << counter << ": ";
        printVec(x1);
        std::cout << std::endl;
        std::cout << "m" << counter++ << ": " << m1 << std::endl;
        m0 = m1;
        x1 = multiplyMatrixVector(A, x1);
        m1 = approach ? x1[idx] / x0[idx] : scalarMult(x1, x0) / scalarMult(x0, x0);
        x0 = x1;
    }
    std::cout << "x" << counter << ": ";
    printVec(x1);
    std::cout << std::endl;
    std::cout << "Found max eigenvalue " << counter << ": " << m1 << std::endl;
    return m1;
    //std::cout << "Min eigenvalue A = " << a - m1 << std::endl;
}

int main()
{
    /*doctest::Context context;
    context.run();

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
        auto result = LinearSystemSolver::UpperRelaxation(A, b, 0.001, 1.5);
        for (size_t i = 0; i < result.size(); ++i) {
            std::cout << "Iteration " << i << ": (" << result.at(i).at(0);
            for (size_t j = 1; j < result.at(i).size(); ++j) {
                std::cout << ", " << result.at(i).at(j);
            }
            std::cout << ")" << std::endl;
        }
    }*/

    std::cout << "===============================================================\nNewton:" << std::endl << std::endl;
    solveNewtonNonLinear(0.001);

}