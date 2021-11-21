#include "Polynomial.h"

#include "doctest.h"

using namespace NumericalAnalysis;

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

TEST_CASE("is increasing, is positive, is decreasing, is negative test") {
    Polynomial pol({
        { 2, 1 },
        { 0, -1 }
    });
    CHECK(pol.isIncreasing({ 0, 2 }, 0.05));
    CHECK_FALSE(pol.isIncreasing({ -2, 0 }, 0.05));
    CHECK_FALSE(pol.isDecreasing({ 0, 2 }, 0.05));
    CHECK(pol.isDecreasing({ -2, 0 }, 0.05));
    CHECK_FALSE(pol.isPositive({ 0, 1 }, 0.05));
    CHECK_FALSE(pol.isPositive({ -1, 0 }, 0.05));
    CHECK(pol.isPositive({ 1, 2 }, 0.05));
    CHECK(pol.isPositive({ -2, -1 }, 0.05));
    CHECK(pol.isNegative({ 0, 1 }, 0.05));
    CHECK(pol.isNegative({ -1, 0 }, 0.05));
    CHECK_FALSE(pol.isNegative({ 1, 2 }, 0.05));
    CHECK_FALSE(pol.isNegative({ -2, -1 }, 0.05));
}

TEST_CASE("change sign test") {
    Polynomial pol({
        { 2, 1 },
        { 0, -1 }
        });
    CHECK(pol.changeSign({ 0, 2 }, 0.05));
    CHECK_FALSE(pol.changeSign({ 2, 3 }, 0.05));
    CHECK(pol.changeSign({ -2, 0 }, 0.05));
    CHECK_FALSE(pol.changeSign({ -3, -2 }, 0.05));
    CHECK_FALSE(pol.changeSign({ -0.5, 0.5 }, 0.05));
}

TEST_CASE("get value test") {
    Polynomial pol({
        { 2, 1 },
        { 0, -1 }
    });
    CHECK_EQ(pol.getValue(0), -1.0);
    CHECK_EQ(pol.getValue(1), 0.0);
    CHECK_EQ(pol.getValue(2), 3.0);
    CHECK_EQ(pol.getValue(-1), 0.0);
    CHECK_EQ(pol.getValue(-2), 3.0);
    CHECK_EQ(pol.getValue(0.5), -0.75);
}

TEST_CASE("min max test") {
    Polynomial pol({
        { 2, 1 },
        { 0, -1 }
    });
    auto min = pol.findAbsMin({ -0.5, 2 }, 0.001);
    CHECK(std::abs(min.first - 1.0) < 0.000001);
    CHECK(min.second < 0.000001);
    auto max = pol.findAbsMax({ 0, 2 }, 0.001);
    CHECK(std::abs(max.first - 2.0) < 0.000001);
    CHECK(std::abs(max.second - 3.0) < 0.000001);
}

TEST_CASE("multiply add test") {
    {
        Polynomial pol1({
            { 3, 10 },
            { 2, 1 },
            { 0, -1 }
        });
        Polynomial pol2(Polynomial::Monomials{ { 0, 1 } });
        CHECK(pol1.equal(Polynomial::multiply(pol1, pol2), 0.0001));
    }
    {
        Polynomial pol1({
            { 10, 2 },
            { 4, -3 },
            { 3, 10 },
            { 2, 1 },
            { 0, -1 }
        });
        Polynomial pol2({
            { 13, 5 },
            { 8, 2 },
            { 6, -30 },
            { 5, 1 },
            { 2, 3 },
            { 0, 5 }
        });
        Polynomial expected({
            { 23, 10 },
            { 18, 4 },
            { 17, -15 },
            { 16, -10 },
            { 15, 7 },
            { 13, -5 },
            { 11, 20 },
            { 10, 102 },
            { 9, -303 },
            { 8, -22 },
            { 7, 1 },
            { 6, 21 },
            { 5, 29 },
            { 4, -12 },
            { 3, 50 },
            { 2, 2 },
            { 0, -5 }
        });
        Polynomial result = Polynomial::multiply(pol1, pol2);
        CHECK(expected.equal(result, 0.0001));
    }
    {
        Polynomial pol1({
            { 3, 10 },
            { 2, 1 },
            { 0, -1 }
            });
        Polynomial pol2(Polynomial::Monomials{ { 0, 1 } });
        Polynomial expected({
            { 3, 10 },
            { 2, 1 }
            });
        CHECK(expected.equal(Polynomial::add(pol1, pol2), 0.0001));
    }
    {
        Polynomial pol1({
            { 10, 2 },
            { 4, -3 },
            { 3, 10 },
            { 2, 1 },
            { 0, -1 }
            });
        Polynomial pol2({
            { 13, 5 },
            { 10, 2 },
            { 6, -30 },
            { 5, 1 },
            { 2, 3 },
            { 0, 5 }
            });
        Polynomial expected({
            { 13, 5 },
            { 10, 4 },
            { 6, -30 },
            { 5, 1 },
            { 4, -3 },
            { 3, 10 },
            { 2, 4 },
            { 0, 4 }
            });
        Polynomial result = Polynomial::add(pol1, pol2);
        CHECK(expected.equal(result, 0.0001));
    }
}