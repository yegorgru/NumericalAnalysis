#include "Math.h"

#include <iostream>
#include <iomanip>

#define DOCTEST_CONFIG_IMPLEMENT
#include "doctest.h"

int main()
{
    doctest::Context context;
    context.run();

    using namespace NumericalCalculus;

    std::cout << "========================================================================" << std::endl
        << "Relaxation:" << std::endl;
    auto result = Math::Relaxation(Polynomial({
        {3, 1},
        {2, -6},
        {1, 5},
        {0, 12},
    }), 0.001, { -2, 0 });
    std::cout << "\t" << std::left << std::setw(20) << "n" << std::setw(20) << "xn" << std::setw(20) << "f(xn)" << std::endl;
    for (size_t i = 0; i < result.size(); i++) {
        std::cout << "\t" << std::left << std::setw(20) << i << std::setw(20) << result[i].first << std::setw(20) << result[i].second << std::endl;
    }
}