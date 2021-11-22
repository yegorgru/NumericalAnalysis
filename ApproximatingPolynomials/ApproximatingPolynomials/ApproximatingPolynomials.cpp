#include "ApproximatingPolynomials.h"

#include "LinearSystemSolver.h"

#include <map>
#include <utility>
#include <vector>
#include <iostream>
#include <iomanip>

namespace NumericalAnalysis
{

Polynomial NumericalAnalysis::ApproximatingPolynomials::NewtonPolynomial(const Polynomial& function, int nodes, bool debug)
{
    std::vector<std::pair<std::pair<double, double>, double>> dividedDifferences;
    double range = 1.0 / ((nodes + 1) / 2);
    if (debug) {
        std::cout << "function values and divided differences:" << std::endl;
    }
    std::vector<double> ran{ 1.0, 2.0, 3.0, 4.0, 5.0 };
    for (double node : ran) {
        dividedDifferences.emplace_back(std::make_pair(node, node), function.getValue(node));
        if (debug) {
            std::cout << "f(" << node << ") = " << dividedDifferences.back().second << std::endl;
        }
    }
    size_t begin = 0;
    for (int iteration = 0; iteration < nodes; iteration++) {
        size_t sizebuf = dividedDifferences.size();
        for (size_t i = begin; i < begin + nodes - iteration - 1; i++) {
            dividedDifferences.emplace_back(std::make_pair(dividedDifferences[i].first.first, dividedDifferences[i+1].first.second),
                                            (dividedDifferences[i + 1].second - dividedDifferences[i].second)
                                                /
                                            (dividedDifferences[i+1].first.second - dividedDifferences[i].first.first)
            );
            if (debug) {
                std::cout << "f(" << dividedDifferences[i].first.first << ", " <<
                    dividedDifferences[i + 1].first.second << ") = " << dividedDifferences.back().second << std::endl;
            }
        }
        begin = sizebuf;
    }
    Polynomial result(Polynomial::Monomials{ {0, dividedDifferences[0].second} });
    size_t left = 0;
    size_t right = 0;
    if (debug) {
        std::cout << "\nPolynomial:\nP = f(" << dividedDifferences[0].first.first << ")";
    }
    for (int iteration = 0; iteration < nodes - 1; iteration++) {
        right += nodes - iteration;
        Polynomial nextPolynomial({ {0, dividedDifferences[right].second} });
        if (debug) {
            std::cout << " + f(" << dividedDifferences[right].first.first << ", " << dividedDifferences[right].first.second << ")";
        }
        for (int i = 0; i < iteration + 1; i++) {
            Polynomial multiplier({ 
                {1, 1},
                {0, -dividedDifferences[i].first.first}
            });
            if (debug) {
                std::cout << "(x - " << dividedDifferences[i].first.first << ")";
            }
            nextPolynomial = Polynomial::multiply(nextPolynomial, multiplier);
        }
        result = Polynomial::add(result, nextPolynomial);
        left = right;
    }
    if (debug) {
        std::cout << std::endl;
    }
    return result;
}

Polynomial NumericalAnalysis::ApproximatingPolynomials::DegreePolynomial(const Polynomial& function, int nodes, bool debug)
{
    if (nodes == 0) return Polynomial();
    std::vector<std::vector<double>> A(nodes, std::vector<double>(nodes, 0.0));
    std::vector<double> b(nodes);
    A[0][0] = 1;
    b[0] = function.getValue(0);
    for (int i = 1; i < nodes; i++) {
        for (int j = 0; j < nodes; j++) {
            A[i][j] = std::pow(i, j);
        }
        b[i] = function.getValue(i);
    }
    if (debug) {
        for (size_t i = 0; i < nodes; i++) {
            for (size_t j = 0; j < nodes; j++) {
                std::cout << std::setw(14) << A[i][j];
            }
            std::cout << std::setw(5) << "|" << std::setw(5) << b[i] << std::endl;
        }
    }
    auto result = LinearSystemSolver::GaussianElimination(A, b);
    Polynomial::Monomials monomials;
    for (int i = 0; i < nodes; i++) {
        monomials[i] = result[i];
    }
    return Polynomial(monomials);
}

}