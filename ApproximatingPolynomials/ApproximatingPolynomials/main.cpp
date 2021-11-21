#include <iostream>
#include <iomanip>

#include "ApproximatingPolynomials.h"

using namespace NumericalAnalysis;

void printPolynomial(Polynomial polynomial) {
    if (polynomial.getMonomials().size() == 0) return;
    std::cout << polynomial.getMonomials().begin()->second << "x^" << polynomial.getMonomials().begin()->first;
    for (auto it = std::next(polynomial.getMonomials().begin()); it != polynomial.getMonomials().end(); ++it) {
        std::cout << " + " << it->second << "x^" << it->first;
    }
    std::cout << std::endl;
}

int main()
{
    std::cout << "Newton Polynomial method: " << std::endl;
    Polynomial p({
        {15, 5},
        {13, 2},
        {12, 7}
    });
    auto result = ApproximatingPolynomials::NewtonPolynomial(p, 12, true);
    std::cout << "\nNewton polynomial result: " << std::endl;
    printPolynomial(result);
    std::cout << "\n\n" << std::setw(20) << "x" << std::setw(20) << "true value" << std::setw(25) << "approximation value" << std::endl;
    double range = 1.0 / 12;
    for (double node = -1.0; node < 0.999; node += range) {
        std::cout << std::setw(20) << node << std::setw(20) << p.getValue(node) << std::setw(25) << result.getValue(node) << std::endl;
    }
}