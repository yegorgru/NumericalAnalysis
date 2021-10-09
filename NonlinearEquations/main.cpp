#include "Math.h"

#define DOCTEST_CONFIG_IMPLEMENT
#include "doctest.h"

int main()
{
    doctest::Context context;
    context.run();

    using namespace NumericalCalculus;
    Point x = Math::Relaxation(Polynomial({
        {3, 1},
        {2, -6},
        {1, 5},
        {0, 12},
    }), 0.001, { -2, 0 });
    std::cout << "Answer: " << x.first << " " << x.second;
}