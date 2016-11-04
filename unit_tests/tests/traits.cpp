/*
 * Quantum++
 *
 * Copyright (c) 2013 - 2017 Vlad Gheorghiu (vgheorgh@gmail.com)
 *
 * This file is part of Quantum++.
 *
 * Quantum++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Quantum++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Quantum++.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <cmath>
#include <vector>
#include <list>
#include "gtest/gtest.h"
#include "qpp.h"

using namespace qpp;

// Unit testing "traits.h"

/******************************************************************************/
/// BEGIN template<typename T> struct qpp::is_complex
TEST(qpp_is_complex, AllTests)
{
    EXPECT_TRUE(qpp::is_complex<std::complex<double>>::value);
    EXPECT_TRUE(qpp::is_complex<std::complex<int>>::value);

    EXPECT_FALSE(qpp::is_complex<double>::value);
    EXPECT_FALSE(qpp::is_complex<cmat>::value);
}
/******************************************************************************/
/// BEGIN template<typename T> struct qpp::is_iterable
TEST(qpp_is_iterable, AllTests)
{
    EXPECT_TRUE(qpp::is_iterable<std::vector<int>>::value);
    EXPECT_TRUE(qpp::is_iterable<std::list<double>>::value);

    class X{};
    EXPECT_FALSE(qpp::is_iterable<X>::value);
    EXPECT_FALSE(qpp::is_iterable<qpp::States>::value);
}
/******************************************************************************/
/// BEGIN template<typename T> struct qpp::is_matrix_expression
TEST(qpp_is_matrix_expression, AllTests)
{
    dmat A, B, C, D;
    int x{}, y{}, z{};

    EXPECT_TRUE(qpp::is_matrix_expression<decltype(3 * A)>::value);
    EXPECT_TRUE(qpp::is_matrix_expression<decltype(A + B)>::value);
    EXPECT_TRUE(qpp::is_matrix_expression<decltype(A + B * C)>::value);
    EXPECT_TRUE(qpp::is_matrix_expression<decltype(D * D * D)>::value);

    EXPECT_FALSE(qpp::is_matrix_expression<decltype(x + y * z)>::value);
}
/******************************************************************************/
