#include "Math.h"

#include <sstream>

#define DOCTEST_CONFIG_IMPLEMENT
#include "doctest.h"

int main()
{
    doctest::Context context;
    context.run();

    using namespace NumericalCalculus;

    std::ostringstream os;
    os << "========================================================================" << std::endl
        << "Relaxation:" << std::endl;
    Point x = Math::Relaxation(Polynomial({
        {3, 1},
        {2, -6},
        {1, 5},
        {0, 12},
    }), 0.001, { -2, 0 }, os);
    os << "Answer: " << x.first << " " << x.second;
    std::cout << os.str();
}