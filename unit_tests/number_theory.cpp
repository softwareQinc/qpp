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

// ********** qpp::modpow() **********
TEST(qpp_modpow_test, PositiveNumbers)
{
    EXPECT_EQ (0, qpp::modpow(2, 3, 4));
    EXPECT_EQ (34, qpp::modpow(17, 176, 37));
    EXPECT_EQ (4042, qpp::modpow(178373, 9281623, 6217));
}
// ********** END qpp::modpow() **********

// ********** qpp::modinv() **********
TEST(qpp_modinv_test, NonNegativeNumbers)
{
    EXPECT_EQ (1, qpp::modinv(1, 1));
    EXPECT_EQ (62, qpp::modinv(2, 123));
    EXPECT_EQ (21, qpp::modinv(1231, 22));
}
// ********** END qpp::modinv() **********

// ********** qpp::egcd() **********
TEST(qpp_egcd_test, NonNegativeNumbers)
{
    EXPECT_EQ (std::make_tuple(1, 0, 10), qpp::egcd(10, 0));
    EXPECT_EQ (std::make_tuple(0, 1, 10), qpp::egcd(0, 10));
    EXPECT_EQ (std::make_tuple(0, 1, 1), qpp::egcd(1, 1));
    EXPECT_EQ (std::make_tuple(-103, 11, 4), qpp::egcd(120, 1124));
    EXPECT_EQ (std::make_tuple(54, -7, 1), qpp::egcd(17, 131));
}

TEST(qpp_egcd_test, MixedNumbers)
{
    EXPECT_EQ (std::make_tuple(-1, 0, 10), qpp::egcd(-10, 0));
    EXPECT_EQ (std::make_tuple(0, -1, 10), qpp::egcd(0, -10));
    EXPECT_EQ (std::make_tuple(0, -1, 1), qpp::egcd(1, -1));
    EXPECT_EQ (std::make_tuple(0, 1, 1), qpp::egcd(-1, 1));
    EXPECT_EQ (std::make_tuple(-103, -11, 4), qpp::egcd(120, -1124));
    EXPECT_EQ (std::make_tuple(-54, -7, 1), qpp::egcd(-17, 131));
}

TEST(qpp_egcd_test, NegativeNumbers)
{
    EXPECT_EQ (std::make_tuple(0, -1, 1), qpp::egcd(-1, -1));
    EXPECT_EQ (std::make_tuple(103, -11, 4), qpp::egcd(-120, -1124));
    EXPECT_EQ (std::make_tuple(-54, 7, 1), qpp::egcd(-17, -131));
}
// ********** END qpp::egcd() **********

