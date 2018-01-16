/*
 * This file is part of Quantum++.
 *
 * MIT License
 *
 * Copyright (c) 2013 - 2018 Vlad Gheorghiu (vgheorgh@gmail.com)
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
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
TEST(qpp_is_complex, AllTests) {
    EXPECT_TRUE(qpp::is_complex<std::complex<double>>::value);
    EXPECT_TRUE(qpp::is_complex<std::complex<int>>::value);

    EXPECT_FALSE(qpp::is_complex<double>::value);
    EXPECT_FALSE(qpp::is_complex<cmat>::value);
}
/******************************************************************************/
/// BEGIN template<typename T> struct qpp::is_iterable
TEST(qpp_is_iterable, AllTests) {
    EXPECT_TRUE(qpp::is_iterable<std::vector<int>>::value);
    EXPECT_TRUE(qpp::is_iterable<std::list<double>>::value);

    class X {};
    EXPECT_FALSE(qpp::is_iterable<X>::value);
    EXPECT_FALSE(qpp::is_iterable<qpp::States>::value);
}
/******************************************************************************/
/// BEGIN template<typename T> struct qpp::is_matrix_expression
TEST(qpp_is_matrix_expression, AllTests) {
    dmat A, B, C, D;
    int x{}, y{}, z{};

    EXPECT_TRUE(qpp::is_matrix_expression<decltype(3 * A)>::value);
    EXPECT_TRUE(qpp::is_matrix_expression<decltype(A + B)>::value);
    EXPECT_TRUE(qpp::is_matrix_expression<decltype(A + B * C)>::value);
    EXPECT_TRUE(qpp::is_matrix_expression<decltype(D * D * D)>::value);

    EXPECT_FALSE(qpp::is_matrix_expression<decltype(x + y * z)>::value);
}
/******************************************************************************/
