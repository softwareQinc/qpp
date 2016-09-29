/*
 * Quantum++
 *
 * Copyright (c) 2013 - 2016 Vlad Gheorghiu (vgheorgh@gmail.com)
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

#include <gtest/gtest.h>
#include <qpp.h>
using namespace qpp;

// Write your unit tests here. Some examples are provided below.

///// BEGIN template<typename InputIterator>
///// typename std::iterator_traits<InputIterator>::value_type
///// qpp::sum(InputIterator first, InputIterator last)
TEST(qpp_sum_test, PositiveNumbers)
{
    std::vector<int> v{0, 1, 2, 3};
    EXPECT_EQ (6, qpp::sum(v.begin(), v.end()));
}

TEST(qpp_sum_test, NegativeNumbers)
{
    std::vector<int> v{0, -1, -2, -3};
    EXPECT_EQ (-6, qpp::sum(v.begin(), v.end()));
}

TEST(qpp_sum_test, MixedNumbers)
{
    std::vector<int> v{ -3, -2, -1, 0, 1, 2};
    EXPECT_EQ (-3, qpp::sum(v.begin(), v.end()));
}
///// END template<typename InputIterator>
///// typename std::iterator_traits<InputIterator>::value_type
///// qpp::sum(InputIterator first, InputIterator last)

///// BEGIN template<typename InputIterator>
///// typename std::iterator_traits<InputIterator>::value_type
///// qpp::prod(InputIterator first, InputIterator last)
TEST(qpp_prod_test, PositiveNumbers)
{
    std::vector<int> v{1, 2, 3, 4};
    EXPECT_EQ (24, qpp::prod(v.begin(), v.end()));
}
///// END template<typename InputIterator>
///// typename std::iterator_traits<InputIterator>::value_type
///// qpp::prod(InputIterator first, InputIterator last)

