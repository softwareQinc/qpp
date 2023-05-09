#include <cmath>
#include <vector>
#include <list>
#include "gtest/gtest.h"
#include "qpp.h"

using namespace qpp;

// Unit testing "traits.hpp"

/******************************************************************************/
/// BEGIN template<typename T> struct is_complex
TEST(qpp_is_complex, AllTests) {
    EXPECT_TRUE(is_complex<std::complex<double>>::value);
    EXPECT_TRUE(is_complex<std::complex<int>>::value);

    EXPECT_FALSE(is_complex<double>::value);
    EXPECT_FALSE(is_complex<cmat>::value);
}
/******************************************************************************/
/// BEGIN template<typename T> struct is_iterable
TEST(qpp_is_iterable, AllTests) {
    EXPECT_TRUE(is_iterable<std::vector<int>>::value);
    EXPECT_TRUE(is_iterable<std::list<double>>::value);

    class X {};
    EXPECT_FALSE(is_iterable<X>::value);
    EXPECT_FALSE(is_iterable<States>::value);
}
/******************************************************************************/
/// BEGIN template<typename Derived> struct is_matrix_expression
TEST(qpp_is_matrix_expression, AllTests) {
    rmat A, B, C, D;
    int x{}, y{}, z{};

    EXPECT_TRUE(is_matrix_expression<decltype(3 * A)>::value);
    EXPECT_TRUE(is_matrix_expression<decltype(A + B)>::value);
    EXPECT_TRUE(is_matrix_expression<decltype(A + B * C)>::value);
    EXPECT_TRUE(is_matrix_expression<decltype(D * D * D)>::value);

    EXPECT_FALSE(is_matrix_expression<decltype(x + y * z)>::value);
}
/******************************************************************************/
