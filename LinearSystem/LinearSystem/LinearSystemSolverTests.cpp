#include "LinearSystemSolver.h"

#include "doctest.h"

#include <cmath>

using namespace NumericalAnalysis;

TEST_CASE("GaussianElimination test") {
    std::vector<std::vector<double>> A;
    std::vector<double> b;
    std::vector<double> expectedResult;
    SUBCASE("1") {
        A = {
            {1.0, 0.0, 0.0},
            {0.0, 1.0, 0.0},
            {0.0, 0.0, 1.0},
        };
        b = { 5.0, 4.0, 3.0 };
        expectedResult = b;
    }
    SUBCASE("2") {
        A = {
            {1.0, 0.0, 0.0},
            {1.0, 1.0, 0.0},
            {0.0, 0.0, 1.0},
        };
        b = { 5.0, 4.0, 3.0 };
        expectedResult = { 5.0, -1.0, 3.0 };
    }
    SUBCASE("3") {
        A = {
            {5.0, 2.0, 1.0, 0.0},
            {1.0, 3.0, 2.0, 8.0},
            {4.0, -6.0, 1.0, 0.0},
            {5.0, 0.0, 3.0, 2.0},
        };
        b = { 14.0, 65.0, -3.0, 32.0 };
        expectedResult = { 1.0, 2.0, 5.0, 6.0 };
    }
    auto result = LinearSystemSolver::GaussianElimination(A, b);
    REQUIRE(result.size() == expectedResult.size());
    for (size_t i = 0; i < expectedResult.size(); ++i) {
        CHECK(std::abs(expectedResult.at(i) - result.at(i)) < 0.00001);
    }
}

TEST_CASE("TridiagonalMatrix test") {
    std::vector<std::vector<double>> A;
    std::vector<double> b;
    std::vector<double> expectedResult;
    SUBCASE("1") {
        A = {
            {1.0, 0.0, 0.0},
            {0.0, 1.0, 0.0},
            {0.0, 0.0, 1.0},
        };
        b = { 5.0, 4.0, 3.0 };
        expectedResult = b;
    }
    SUBCASE("2") {
        A = {
            {1.0, 0.0, 0.0},
            {1.0, 1.0, 0.0},
            {0.0, 0.0, 1.0},
        };
        b = { 5.0, 4.0, 3.0 };
        expectedResult = { 5.0, -1.0, 3.0 };
    }
    SUBCASE("3") {
        A = {
            {2.0, 4.0, 0.0},
            {4.0, 1.0, 5.0},
            {0.0, 5.0, 2.0},
        };
        b = { 20.0, 37.0, 30.0 };
        expectedResult = { 2.0, 4.0, 5.0 };
    }
    SUBCASE("3") {
        A = {
            {7.0, 4.0, 0.0},
            {1.0, 3.0, 1.0},
            {0.0, 2.0, 5.0},
        };
        b = { 30.0, 19.0, 33.0 };
        expectedResult = { 2.0, 4.0, 5.0 };
    }
    auto result = LinearSystemSolver::TridiagonalMatrix(A, b);
    REQUIRE(result.size() == expectedResult.size());
    for (size_t i = 0; i < expectedResult.size(); ++i) {
        CHECK(std::abs(expectedResult.at(i) - result.at(i)) < 0.00001);
    }
}

TEST_CASE("Jacobi test") {
    std::vector<std::vector<double>> A;
    std::vector<double> b;
    std::vector<double> expectedResult;
    double precision = 0.001;
    SUBCASE("1") {
        A = {
            {1.0, 0.0, 0.0},
            {0.0, 1.0, 0.0},
            {0.0, 0.0, 1.0},
        };
        b = { 5.0, 4.0, 3.0 };
        expectedResult = b;
    }
    SUBCASE("2") {
        A = {
            {1.0, 0.0, 0.0},
            {1.0, 1.0, 0.0},
            {0.0, 0.0, 1.0},
        };
        b = { 5.0, 4.0, 3.0 };
        expectedResult = { 5.0, -1.0, 3.0 };
    }
    SUBCASE("3") {
        A = {
            {6.0, 3.0, 1.0, 0.0},
            {3.0, 5.0, 0.0, 2.0},
            {1.0, 0.0, 3.0, 1.0},
            {0.0, 2.0, 1.0, 5.0},
        };
        b = { 25.0, 31.0, 19.0, 35.0 };
        expectedResult = { 2.0, 3.0, 4.0, 5.0 };
    }
    auto result = LinearSystemSolver::Jacobi(A, b, precision);
    REQUIRE(result.back().size() == expectedResult.size());
    for (size_t i = 0; i < expectedResult.size(); ++i) {
        CHECK(std::abs(expectedResult.at(i) - result.back().at(i)) < precision);
    }
}

TEST_CASE("UpperRelaxation test") {
    std::vector<std::vector<double>> A;
    std::vector<double> b;
    std::vector<double> expectedResult;
    double precision = 0.001;
    SUBCASE("1") {
        A = {
            {1.0, 0.0, 0.0},
            {0.0, 1.0, 0.0},
            {0.0, 0.0, 1.0},
        };
        b = { 5.0, 4.0, 3.0 };
        expectedResult = b;
    }
    SUBCASE("2") {
        A = {
            {1.0, 0.0, 0.0},
            {1.0, 1.0, 0.0},
            {0.0, 0.0, 1.0},
        };
        b = { 5.0, 4.0, 3.0 };
        expectedResult = { 5.0, -1.0, 3.0 };
    }
    SUBCASE("3") {
        A = {
            {5.0, 1.0, 1.0, 0.0},
            {1.0, 2.0, 0.0, 0.0},
            {1.0, 0.0, 4.0, 2.0},
            {0.0, 0.0, 2.0, 3.0},
        };
        b = { 17.0, 8.0, 28.0, 23.0 };
        expectedResult = { 2.0, 3.0, 4.0, 5.0 };
    }
    auto result = LinearSystemSolver::UpperRelaxation(A, b, precision, 1.5);
    REQUIRE(result.back().size() == expectedResult.size());
    for (size_t i = 0; i < expectedResult.size(); ++i) {
        CHECK(std::abs(expectedResult.at(i) - result.back().at(i)) < precision);
    }
}

TEST_CASE("Seidel test") {
    std::vector<std::vector<double>> A;
    std::vector<double> b;
    std::vector<double> expectedResult;
    double precision = 0.001;
    SUBCASE("1") {
        A = {
            {1.0, 0.0, 0.0},
            {0.0, 1.0, 0.0},
            {0.0, 0.0, 1.0},
        };
        b = { 5.0, 4.0, 3.0 };
        expectedResult = b;
    }
    SUBCASE("2") {
        A = {
            {1.0, 0.0, 0.0},
            {1.0, 1.0, 0.0},
            {0.0, 0.0, 1.0},
        };
        b = { 5.0, 4.0, 3.0 };
        expectedResult = { 5.0, -1.0, 3.0 };
    }
    SUBCASE("3") {
        A = {
            {5.0, 1.0, 1.0, 0.0},
            {1.0, 2.0, 0.0, 0.0},
            {1.0, 0.0, 4.0, 2.0},
            {0.0, 0.0, 2.0, 3.0},
        };
        b = { 17.0, 8.0, 28.0, 23.0 };
        expectedResult = { 2.0, 3.0, 4.0, 5.0 };
    }
    auto result = LinearSystemSolver::Seidel(A, b, precision);
    REQUIRE(result.back().size() == expectedResult.size());
    for (size_t i = 0; i < expectedResult.size(); ++i) {
        CHECK(std::abs(expectedResult.at(i) - result.back().at(i)) < precision);
    }
}