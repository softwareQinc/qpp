#include <catch2/catch_all.hpp>

#include <qpp/qpp.hpp>

int fun(int x) { return x * x; }

TEST_CASE("Fun with Catch2 - I") {
    REQUIRE(fun(42) == 1764);
    REQUIRE(fun(43) == 1849);

    BENCHMARK("Fun") { return fun(20); };
}
