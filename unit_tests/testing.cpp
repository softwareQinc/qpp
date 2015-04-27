/*
 * Quantum++
 *
 * Copyright (c) 2013 - 2015 Vlad Gheorghiu (vgheorgh@gmail.com)
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
// #include <MATLAB/matlab.h> // support for MATLAB
// #include <experimental/test.h> // support for testing features

using namespace qpp;

// TODO: write your unit tests here

// ********** qpp::sum() **********
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
// ********** END qpp::sum() **********

// ********** qpp::prod() **********
TEST(qpp_prod_test, PositiveNumbers)
{
    std::vector<int> v{1, 2, 3, 4};
    EXPECT_EQ (24, qpp::prod(v.begin(), v.end()));
}
// ********** END qpp::prod() **********

// ********** qpp::modpow() **********
TEST(qpp_modpow_test, PositiveNumbers)
{
    EXPECT_EQ (0, qpp::modpow(2, 3, 4));
    EXPECT_EQ (34, qpp::modpow(17, 176, 37));
    EXPECT_EQ (4042, qpp::modpow(178373, 9281623, 6217));
}
// ********** END qpp::modpow() **********

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
