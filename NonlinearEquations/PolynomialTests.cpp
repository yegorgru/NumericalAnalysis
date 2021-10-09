#include "Polynomial.h"

#include "doctest.h"

using namespace NumericalCalculus;

TEST_CASE("take derivative test") {
    Polynomial polynomial;
    Polynomial expectedDerivative;
    SUBCASE("1 degree") {
        polynomial = Polynomial({
            {1, 2.0}
        });
        expectedDerivative = Polynomial({
            {0, 2.0}
        });
    }
    SUBCASE("2 degree") {
        polynomial = Polynomial({
            {1, 2.0},
            {2, 5.0}
        });
        expectedDerivative = Polynomial({
            {0, 2.0},
            {1, 10.0}
        });
    }
    SUBCASE("3 degree") {
        polynomial = Polynomial({
            {1, 2.0},
            {2, 5.5},
            {3, 9.0}
        });
        expectedDerivative = Polynomial({
            {0, 2.0},
            {1, 11.0},
            {2, 27.0}
        });
    }
    SUBCASE("3 degree, more complicated") {
        polynomial = Polynomial({
            {0, 12.0},
            {1, 5.0},
            {2, -6.0},
            {3, 1.0}
        });
        expectedDerivative = Polynomial({
            {0, 5.0},
            {1, -12.0},
            {2, 3.0}
        });
    }
    SUBCASE("gaps") {
        polynomial = Polynomial({
            {0, 12.0},
            {1, 5.0},
            {3, 3.0}
        });
        expectedDerivative = Polynomial({
            {0, 5.0},
            {2, 9.0}
        });
    }
    SUBCASE("negative degree") {
        polynomial = Polynomial({
            {-2, -3.5},
            {-1, 5.0},
            {0, 12.0},
            {1, 5.0},
            {3, 3.0}
        });
        expectedDerivative = Polynomial({
            {-3, 7.0},
            {-2, -5.0},
            {0, 5.0},
            {2, 9.0}
        });
    }
    CHECK(polynomial.takeDerivative().getMonomials() == expectedDerivative.getMonomials());
}