#include "Math.h"

#include "doctest.h"

using namespace NumericalCalculus;

TEST_CASE("own example") {
    double precision = 0.001;
    Iterations result;
    Polynomial function = Polynomial({
        {2, 1},
        {0, -1},
    });
    SUBCASE("relaxation") {
        result = Math::Relaxation(function, precision, { 0.5, 1.5 });
    }
    SUBCASE("newton") {
        result = Math::Newton(function, precision, { 0.5, 1.5 });
    }
    SUBCASE("secant") {
        result = Math::Secant(function, precision, { 0.5, 1.5 });
    }
    REQUIRE(!result.empty());
    CHECK(std::abs(result.back().first - 1) < 10 * precision);
    CHECK(result.back().second < precision);
}

TEST_CASE("own example negative") {
    double precision = 0.001;
    Iterations result;
    Polynomial function = Polynomial({
        {2, 1},
        {0, -1},
        });
    SUBCASE("relaxation") {
        result = Math::Relaxation(function, precision, { -1.5, -0.5 });
    }
    SUBCASE("newton") {
        result = Math::Newton(function, precision, { -1.5, -0.5 });
    }
    SUBCASE("secant") {
        result = Math::Secant(function, precision, { -1.5, -0.5 });
    }
    REQUIRE(!result.empty());
    CHECK(std::abs(result.back().first + 1) < 10 * precision);
    CHECK(result.back().second < precision);
}

TEST_CASE("first function") {
    double precision = 0.001;
    Iterations result;
    Polynomial function = Polynomial({
        {3, 1},
        {2, -6},
        {1, 5},
        {0, 12},
    });
    SUBCASE("relaxation") {
        result = Math::Relaxation(function, precision, { -2, 0 });
    }
    SUBCASE("newton") {
        result = Math::Newton(function, precision, { -1.3, -0.7 });
    }
    SUBCASE("secant") {
        result = Math::Secant(function, precision, { -1.3, -0.7 });
    }
    REQUIRE(!result.empty());
    CHECK(std::abs(result.back().first + 1) < 10 * precision);
    CHECK(result.back().second < precision);
}

TEST_CASE("second function") {
    double precision = 0.001;
    Iterations result;
    Polynomial function = Polynomial({
        {3, 1},
        {2, 3},
        {1, -1},
        {0, -3},
    });
    SUBCASE("relaxation") {
        result = Math::Relaxation(function, precision, { 0, 2 });
    }
    SUBCASE("newton") {
        result = Math::Newton(function, precision, { 0.6, 1.4 });
    }
    SUBCASE("secant") {
        result = Math::Secant(function, precision, { 0.6, 1.4 });
    }
    REQUIRE(!result.empty());
    CHECK(std::abs(result.back().first - 1) < 10 * precision);
    CHECK(result.back().second < precision);
}

TEST_CASE("third function") {
    double precision = 0.001;
    Iterations result;
    Polynomial function = Polynomial({
        {3, 1},
        {2, 1},
        {1, -4},
        {0, -4},
    });
    SUBCASE("relaxation") {
        result = Math::Relaxation(function, precision, { 1, 3 });
    }
    SUBCASE("newton") {
        result = Math::Newton(function, precision, { 1.5, 2.5 });
    }
    SUBCASE("secant") {
        result = Math::Secant(function, precision, { 1.5, 2.5 });
    }
    REQUIRE(!result.empty());
    CHECK(std::abs(result.back().first - 2) < 10 * precision);
    CHECK(result.back().second < precision);
}