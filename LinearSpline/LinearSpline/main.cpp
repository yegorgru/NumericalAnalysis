#include <iostream>

#include "SplineFactory.h"

#include <iomanip>

int main()
{
    using namespace NumericalAnalysis;
    Polynomial function({
        {8, 1},
        {3, -1},
        {2, 2},
        {0, -2}
    });
    auto spline = SplineFactory::createLinearSpline(function, { 1.0, 6.0 }, 6, true);
    std::cout << "\n\n" << std::setw(20) << "x" << std::setw(20) << "true value" << std::setw(25) << "approximation value" << std::endl;
    double x = 1.0;
    std::cout << std::setw(20) << x << std::setw(20) << function.getValue(x) << std::setw(25) << spline.begin()->second.getValue(x) << std::endl;
    for (const auto& p : spline) {
        x += 0.5;
        std::cout << std::setw(20) << x << std::setw(20) << function.getValue(x) << std::setw(25) << p.second.getValue(x) << std::endl;
        x += 0.5;
        std::cout << std::setw(20) << x << std::setw(20) << function.getValue(x) << std::setw(25) << p.second.getValue(x) << std::endl;
    }
}