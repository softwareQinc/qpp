#include <catch2/catch_all.hpp>

#include <qpp/qpp.hpp>

int fun2(int x) { return x * x; }

TEST_CASE("Fun with Catch2 - II") {
    REQUIRE(fun2(42) == 1764);
    REQUIRE(fun2(43) == 1849);

    BENCHMARK("Fun2") { return fun2(20); };
}
