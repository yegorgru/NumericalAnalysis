#include "Math.h"

#include "doctest.h"

using namespace NumericalAnalysis;

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

TEST_CASE("test Newton") {
    double precision = 0.001;
    Iterations result;
    Polynomial function;
    Interval interval;
    double expectedValue = 0.0;
    SUBCASE("1") {
        function = Polynomial({
            {3, 1},
            {2, 2},
            {1, -1},
            {0, -2},
        });
        interval = { 0.6, 1.4 };
        expectedValue = 1;
    }
    SUBCASE("5") {
        function = Polynomial({
            {3, 1},
            {2, -6},
            {1, 5},
            {0, 12},
        });
        interval = { -1.5, -0.5 };
        expectedValue = -1;
    }
    SUBCASE("7") {
        function = Polynomial({
            {3, 1},
            {2, -4},
            {1, -4},
            {0, 16},
        });
        interval = { -2.5, -1.5 };
        expectedValue = -2;
    }
    SUBCASE("9") {
        function = Polynomial({
            {3, 1},
            {2, -7},
            {1, 7},
            {0, 15},
        });
        interval = { -1.5, -0.5 };
        expectedValue = -1;
    }
    result = Math::Newton(function, precision, interval);
    REQUIRE(!result.empty());
    CHECK(std::abs(result.back().first - expectedValue) < 10 * precision);
    CHECK(result.back().second < precision);
}

TEST_CASE("test secant") {
    double precision = 0.001;
    Iterations result;
    Polynomial function;
    Interval interval;
    double expectedValue = 0.0;
    SUBCASE("1") {
        function = Polynomial({
            {3, 1},
            {2, -2},
            {1, -1},
            {0, 2},
        });
        interval = { -1.3, -0.7 };
        expectedValue = -1;
    }
    SUBCASE("6") {
        function = Polynomial({
            {3, 1},
            {2, 3},
            {1, -1},
            {0, -3},
        });
        interval = { 0.7, 1.3 };
        expectedValue = 1;
    }
    SUBCASE("9") {
        function = Polynomial({
            {3, 1},
            {2, -5},
            {1, -4},
            {0, 20},
            });
        interval = { -2.5, -1.5 };
        expectedValue = -2;
    }
    result = Math::Secant(function, precision, interval);
    REQUIRE(!result.empty());
    CHECK(std::abs(result.back().first - expectedValue) < 10 * precision);
    CHECK(result.back().second < precision);
}

TEST_CASE("test relaxation") {
    double precision = 0.001;
    Iterations result;
    Polynomial function;
    Interval interval;
    double expectedValue = 0.0;
    SUBCASE("2, 3") {
        function = Polynomial({
            {3, 1},
            {2, -4},
            {1, 1},
            {0, 6},
        });
        interval = { -2.0, 0.0 };
        expectedValue = -1;
    }
    SUBCASE("4") {
        function = Polynomial({
            {3, 1},
            {2, 3},
            {1, -1},
            {0, -3},
        });
        interval = { 0.0, 2.0 };
        expectedValue = 1;
    }
    SUBCASE("6") {
        function = Polynomial({
            {3, 1},
            {2, 1},
            {1, -4},
            {0, -4},
            });
        interval = { 1.0, 3.0 };
        expectedValue = 2;
    }
    SUBCASE("6") {
        function = Polynomial({
            {3, 1},
            {2, -7},
            {1, 7},
            {0, 15},
            });
        interval = { -2.0, 0.0 };
        expectedValue = -1;
    }
    SUBCASE("10") {
        function = Polynomial({
            {3, 1},
            {2, -8},
            {1, 9},
            {0, 18},
            });
        interval = { -2.0, 0.0 };
        expectedValue = -1;
    }
    result = Math::Relaxation(function, precision, interval);
    REQUIRE(!result.empty());
    CHECK(std::abs(result.back().first - expectedValue) < 10 * precision);
    CHECK(result.back().second < precision);
}

TEST_CASE("test ModifiedNewton") {
    double precision = 0.001;
    Iterations result;
    Polynomial function;
    Interval interval;
    double expectedValue = 0.0;
    SUBCASE("1") {
        function = Polynomial({
            {3, 1},
            {2, 2},
            {1, -1},
            {0, -2},
            });
        interval = { 0.6, 1.4 };
        expectedValue = 1;
    }
    SUBCASE("5") {
        function = Polynomial({
            {3, 1},
            {2, -6},
            {1, 5},
            {0, 12},
            });
        interval = { -1.5, -0.5 };
        expectedValue = -1;
    }
    SUBCASE("7") {
        function = Polynomial({
            {3, 1},
            {2, -4},
            {1, -4},
            {0, 16},
            });
        interval = { -2.5, -1.5 };
        expectedValue = -2;
    }
    SUBCASE("9") {
        function = Polynomial({
            {3, 1},
            {2, -7},
            {1, 7},
            {0, 15},
            });
        interval = { -1.5, -0.5 };
        expectedValue = -1;
    }
    SUBCASE("var4") {
        function = Polynomial({
            {3, 1},
            {2, -6},
            {1, 5},
            {0, 12},
            });
        interval = { -1.5, -0.5 };
        expectedValue = -1;
    }
    result = Math::ModifiedNewton(function, precision, interval);
    REQUIRE(!result.empty());
    CHECK(std::abs(result.back().first - expectedValue) < 10 * precision);
    CHECK(result.back().second < precision);
}