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

#include <limits>
#include <string>
#include <tuple>
#include <vector>
#include "gtest/gtest.h"
#include "qpp.h"

using namespace qpp;

// Unit testing "number_theory.h"

/******************************************************************************/
/// BEGIN std::vector<idx> qpp::compperm(const std::vector<idx>& perm,
///       const std::vector<idx>& sigma)
TEST(qpp_compperm, NonNegativeNumbers) {
    std::vector<idx> sigma1{0};
    std::vector<idx> tau1{0};
    std::vector<idx> result1{0};
    EXPECT_EQ(result1, qpp::compperm(sigma1, tau1));

    std::vector<idx> sigma2{0, 1, 2};
    std::vector<idx> tau2{0, 1, 2};
    std::vector<idx> result2{0, 1, 2};
    EXPECT_EQ(result2, qpp::compperm(sigma2, tau2));

    std::vector<idx> sigma3{3, 0, 1, 4, 5, 2, 6};
    std::vector<idx> tau3{1, 2, 5, 0, 3, 4, 6}; // inverse permutation
    std::vector<idx> result3{0, 1, 2, 3, 4, 5, 6};
    EXPECT_EQ(result3, qpp::compperm(sigma3, tau3));
}
/******************************************************************************/
/// BEGIN double qpp::contfrac2x(const std::vector<int>& cf, idx N = idx(-1))
TEST(qpp_contfrac2x, DefaultN) {
    EXPECT_NEAR(0, qpp::contfrac2x({0}), 1e-7);
    EXPECT_NEAR(42, qpp::contfrac2x({42}), 1e-7);
    EXPECT_NEAR(-42, qpp::contfrac2x({-42}), 1e-7);
    EXPECT_NEAR(qpp::pi, qpp::contfrac2x({3, 7, 15, 1, 292, 1}), 1e-5);
    EXPECT_NEAR(-qpp::pi, qpp::contfrac2x({-3, -7, -15, -1, -292, -1}), 1e-5);
    EXPECT_NEAR(0.1234, qpp::contfrac2x({0, 8, 9, 1, 1, 1}), 1e-5);
    EXPECT_NEAR(-0.4321, qpp::contfrac2x({0, -2, -3, -5, -2}), 1e-5);
    EXPECT_NEAR(std::sqrt(19), qpp::contfrac2x({4, 2, 1, 3, 1, 2, 8}), 1e-3);
}

TEST(qpp_contfrac2x, UserSpecifiedN) {
    EXPECT_NEAR(0, qpp::contfrac2x({0}, 1), 1e-7);
    EXPECT_NEAR(42, qpp::contfrac2x({42}, 1), 1e-7);
    EXPECT_NEAR(-42, qpp::contfrac2x({-42}, 1), 1e-7);
    EXPECT_NEAR(3.141, qpp::contfrac2x({3, 7, 15, 1, 292, 1}, 3), 1e-2);
    EXPECT_NEAR(-3.141, qpp::contfrac2x({-3, -7, -15, -1, -292, -1}, 3), 1e-2);
    EXPECT_NEAR(0.1234, qpp::contfrac2x({0, 8, 9, 1, 1, 1}, 4), 1e-3);
    EXPECT_NEAR(-0.4321, qpp::contfrac2x({0, -2, -3, -5, -2}, 4), 1e-3);
}
/******************************************************************************/
/// BEGIN std::tuple<bigint, bigint, bigint> qpp::egcd(bigint m, bigint n)
TEST(qpp_egcd, NonNegativeNumbers) {
    EXPECT_EQ(std::make_tuple(1, 0, 10), qpp::egcd(10, 0));
    EXPECT_EQ(std::make_tuple(0, 1, 10), qpp::egcd(0, 10));
    EXPECT_EQ(std::make_tuple(0, 1, 1), qpp::egcd(1, 1));
    EXPECT_EQ(std::make_tuple(-103, 11, 4), qpp::egcd(120, 1124));
    EXPECT_EQ(std::make_tuple(54, -7, 1), qpp::egcd(17, 131));
}

TEST(qpp_egcd, MixedNumbers) {
    EXPECT_EQ(std::make_tuple(-1, 0, 10), qpp::egcd(-10, 0));
    EXPECT_EQ(std::make_tuple(0, -1, 10), qpp::egcd(0, -10));
    EXPECT_EQ(std::make_tuple(0, -1, 1), qpp::egcd(1, -1));
    EXPECT_EQ(std::make_tuple(0, 1, 1), qpp::egcd(-1, 1));
    EXPECT_EQ(std::make_tuple(-103, -11, 4), qpp::egcd(120, -1124));
    EXPECT_EQ(std::make_tuple(-54, -7, 1), qpp::egcd(-17, 131));
}

TEST(qpp_egcd, NegativeNumbers) {
    EXPECT_EQ(std::make_tuple(0, -1, 1), qpp::egcd(-1, -1));
    EXPECT_EQ(std::make_tuple(103, -11, 4), qpp::egcd(-120, -1124));
    EXPECT_EQ(std::make_tuple(-54, 7, 1), qpp::egcd(-17, -131));
}
/******************************************************************************/
/// BEGIN std::vector<bigint> qpp::factors(bigint n)
TEST(qpp_factors, MixedNumbers) {
    bigint n1 = 2;
    std::vector<bigint> result1{2};
    EXPECT_EQ(result1, qpp::factors(n1));

    bigint n2 = -10;
    std::vector<bigint> result2{2, 5};
    EXPECT_EQ(result2, qpp::factors(n2));

    bigint n3 = 110500;
    std::vector<bigint> result3{2, 2, 5, 5, 5, 13, 17};
    EXPECT_EQ(result3, qpp::factors(n3));

    bigint n4 = -35750;
    std::vector<bigint> result4{2, 5, 5, 5, 11, 13};
    EXPECT_EQ(result4, qpp::factors(n4));

    // test some really large numbers
    bigint n5 = 1000000000001112;
    std::vector<bigint> result5{2, 2, 2, 3, 11, 29, 1217, 1637, 65563};
    EXPECT_EQ(result5, qpp::factors(n5));

    bigint n6 = 1000000000001118;
    std::vector<bigint> result6{2, 3, 897649, 185670197};
    EXPECT_EQ(result6, qpp::factors(n6));
}
/******************************************************************************/
/// BEGIN bigint qpp::gcd(bigint m, bigint n)
TEST(qpp_gcd, NonNegativeNumbers) {
    EXPECT_EQ(10, qpp::gcd(10, 0));
    EXPECT_EQ(10, qpp::gcd(0, 10));
    EXPECT_EQ(1, qpp::gcd(1, 1));
    EXPECT_EQ(4, qpp::gcd(120, 1124));
    EXPECT_EQ(1, qpp::gcd(17, 131));

    // test some really large numbers
    // assumes that qpp::bigint is 64 bits
    bigint maxbigint = std::numeric_limits<bigint>::max();
    EXPECT_EQ(maxbigint, qpp::gcd(maxbigint, maxbigint));
    EXPECT_EQ(1, qpp::gcd(maxbigint, maxbigint - 1));
    EXPECT_EQ(3, qpp::gcd(maxbigint - 1, maxbigint - 10));
}

TEST(qpp_gcd, MixedNumbers) {
    EXPECT_EQ(10, qpp::gcd(-10, 0));
    EXPECT_EQ(10, qpp::gcd(0, -10));
    EXPECT_EQ(1, qpp::gcd(1, -1));
    EXPECT_EQ(1, qpp::gcd(-1, 1));
    EXPECT_EQ(4, qpp::gcd(120, -1124));
    EXPECT_EQ(1, qpp::gcd(-17, 131));
}

TEST(qpp_gcd, NegativeNumbers) {
    EXPECT_EQ(1, qpp::gcd(-1, -1));
    EXPECT_EQ(4, qpp::gcd(-120, -1124));
    EXPECT_EQ(1, qpp::gcd(-17, -131));
    EXPECT_EQ(20, qpp::gcd(-170180, -13100));
}
/******************************************************************************/
/// BEGIN bigint qpp::gcd(const std::vector<bigint>& ns)
TEST(qpp_gcd_list, MixedNumbers) {
    std::vector<bigint> v1{1, 0, -1, 0, 1};
    EXPECT_EQ(1, qpp::gcd(v1));

    std::vector<bigint> v2{1, 2, 3, 4, 5};
    EXPECT_EQ(1, qpp::gcd(v2));

    std::vector<bigint> v3{125000, -12680, 25040, -4};
    EXPECT_EQ(4, qpp::gcd(v3));
}
/******************************************************************************/
/// BEGIN std::vector<idx> qpp::invperm(const std::vector<idx>& perm)
TEST(qpp_invperm, NonNegativeNumbers) {
    std::vector<idx> perm1{0};
    std::vector<idx> result1{0};
    EXPECT_EQ(result1, qpp::invperm(perm1));

    std::vector<idx> perm2{0, 1, 2, 3, 4, 5, 6};
    std::vector<idx> result2{0, 1, 2, 3, 4, 5, 6};
    EXPECT_EQ(result2, qpp::invperm(perm2));

    std::vector<idx> perm3{6, 5, 4, 3, 2, 1, 0};
    std::vector<idx> result3{6, 5, 4, 3, 2, 1, 0};
    EXPECT_EQ(result3, qpp::invperm(perm3));

    std::vector<idx> perm4{3, 0, 1, 4, 5, 2, 6};
    std::vector<idx> result4{1, 2, 5, 0, 3, 4, 6};
    EXPECT_EQ(result4, qpp::invperm(perm4));
}
/******************************************************************************/
/// BEGIN bool qpp::isprime(bigint n, idx k = 80)
TEST(qpp_isprime, MixedPrimeNumbers) {
    EXPECT_TRUE(qpp::isprime(2));
    EXPECT_TRUE(qpp::isprime(-2));
    EXPECT_TRUE(qpp::isprime(11));
    EXPECT_TRUE(qpp::isprime(-17));
    EXPECT_TRUE(qpp::isprime(19));
    EXPECT_TRUE(qpp::isprime(127));
    EXPECT_TRUE(qpp::isprime(541));
    EXPECT_TRUE(qpp::isprime(-104729));
    EXPECT_TRUE(qpp::isprime(10000000019));
}

TEST(qpp_isprime, MixedNonPrimeNumbers) {
    EXPECT_FALSE(qpp::isprime(4));
    EXPECT_FALSE(qpp::isprime(-4));
    EXPECT_FALSE(qpp::isprime(110));
    EXPECT_FALSE(qpp::isprime(-4891));
    EXPECT_FALSE(qpp::isprime(-13101));
    EXPECT_FALSE(qpp::isprime(110011));
    EXPECT_FALSE(qpp::isprime(10000000119));

    // test some Charmichael numbers
    EXPECT_FALSE(qpp::isprime(561));
    EXPECT_FALSE(qpp::isprime(6601));
    EXPECT_FALSE(qpp::isprime(8911));
    EXPECT_FALSE(qpp::isprime(41041));
}
/******************************************************************************/
/// BEGIN bigint qpp::lcm(bigint m, bigint n)
TEST(qpp_lcm, NonNegativeNumbers) {
    EXPECT_EQ(0, qpp::lcm(10, 0));
    EXPECT_EQ(0, qpp::lcm(0, 10));
    EXPECT_EQ(1, qpp::lcm(1, 1));
    EXPECT_EQ(33720, qpp::lcm(120, 1124));
    EXPECT_EQ(2227, qpp::lcm(17, 131));
}

TEST(qpp_lcm, MixedNumbers) {
    EXPECT_EQ(0, qpp::lcm(-10, 0));
    EXPECT_EQ(0, qpp::lcm(0, -10));
    EXPECT_EQ(1, qpp::lcm(1, -1));
    EXPECT_EQ(1, qpp::lcm(-1, 1));
    EXPECT_EQ(33720, qpp::lcm(120, -1124));
    EXPECT_EQ(2227, qpp::lcm(-17, 131));
}

TEST(qpp_lcm, NegativeNumbers) {
    EXPECT_EQ(1, qpp::lcm(-1, -1));
    EXPECT_EQ(33720, qpp::lcm(-120, -1124));
    EXPECT_EQ(2227, qpp::lcm(-17, -131));
    EXPECT_EQ(111467900, qpp::lcm(-170180, -13100));
}
/******************************************************************************/
/// BEGIN bigint qpp::lcm(const std::vector<bigint>& ns)
TEST(qpp_lcm_list, MixedNumbers) {
    std::vector<bigint> v1{10, -10, 10, -10};
    EXPECT_EQ(10, qpp::lcm(v1));

    std::vector<bigint> v2{1, 2, 3, 4, 5};
    EXPECT_EQ(60, qpp::lcm(v2));

    std::vector<bigint> v3{12500, -1268, 2504, -4};
    EXPECT_EQ(2480525000, qpp::lcm(v3));
}
/******************************************************************************/
/// BEGIN bigint qpp::modinv(bigint a, bigint p)
TEST(qpp_modinv, NonNegativeNumbers) {
    EXPECT_EQ(1, qpp::modinv(1, 1));
    EXPECT_EQ(62, qpp::modinv(2, 123));
    EXPECT_EQ(21, qpp::modinv(1231, 22));
}
/******************************************************************************/
/// BEGIN bigint qpp::modmul(bigint a, bigint n, bigint p)
TEST(qpp_modmul, NonNegativeNumbers) {
    EXPECT_EQ(0, qpp::modmul(0, 0, 1));
    EXPECT_EQ(0, qpp::modmul(0, 0, 2));
    EXPECT_EQ(0, qpp::modmul(1, 0, 2));
    EXPECT_EQ(0, qpp::modmul(0, 1, 2));
    EXPECT_EQ(2611, qpp::modmul(12127, 71623, 12345));

    // test some really large numbers
    // assumes that qpp::bigint is 64 bits
    bigint maxbigint = std::numeric_limits<bigint>::max();

    EXPECT_EQ(42, qpp::modmul(maxbigint - 1, maxbigint, 123));
    EXPECT_EQ(49, qpp::modmul(maxbigint, maxbigint, 123));
    EXPECT_EQ(2262, qpp::modmul(maxbigint - 189, maxbigint - 2345, 7891));
    EXPECT_EQ(0, qpp::modmul(maxbigint, maxbigint, maxbigint));
    EXPECT_EQ(14884,
              qpp::modmul(maxbigint - 1, maxbigint - 1, maxbigint - 123));
}

TEST(qpp_modmul, MixedNumbers) {
    EXPECT_EQ(0, qpp::modmul(-1, 0, 2));
    EXPECT_EQ(0, qpp::modmul(0, -1, 2));
    EXPECT_EQ(9734, qpp::modmul(12127, -71623, 12345));
    EXPECT_EQ(9734, qpp::modmul(-12127, 71623, 12345));

    // test some really large numbers
    // assumes that qpp::bigint is 64 bits
    bigint minbigint = std::numeric_limits<bigint>::min();
    bigint maxbigint = std::numeric_limits<bigint>::max();

    EXPECT_EQ(1114, qpp::modmul(minbigint, maxbigint, 2314));
    EXPECT_EQ(21, qpp::modmul(-maxbigint, maxbigint, 34));
    EXPECT_EQ(240, qpp::modmul(maxbigint, -1234567890, 314));
    EXPECT_EQ(219, qpp::modmul(maxbigint, (maxbigint - 1) / 2, 1314));
}

TEST(qpp_modmul, NegativeNumbers) {
    // test some really large numbers
    // assumes that qpp::bigint is 64 bits
    bigint minbigint = std::numeric_limits<bigint>::min();
    bigint maxbigint = std::numeric_limits<bigint>::max();
    EXPECT_EQ(13, qpp::modmul(-maxbigint, -maxbigint, 34));
    EXPECT_EQ(1, qpp::modmul(-maxbigint, -maxbigint, maxbigint - 1));
    EXPECT_EQ(74, qpp::modmul(-maxbigint, -1234567890, 314));
    EXPECT_EQ(1799, qpp::modmul(minbigint + 1234, -maxbigint + 2345, 7891));

    EXPECT_EQ(64, qpp::modmul(minbigint, minbigint, 123));
    EXPECT_EQ(56, qpp::modmul(minbigint + 1, minbigint, 123));
    EXPECT_EQ(1799, qpp::modmul(minbigint + 1234, -maxbigint + 2345, 7891));
}
/******************************************************************************/
/// BEGIN bigint qpp::modpow(bigint a, bigint n, bigint p)
TEST(qpp_modpow, PositiveNumbers) {
    bigint maxbigint = std::numeric_limits<bigint>::max();

    EXPECT_EQ(0, qpp::modpow(0, 100, 1));
    EXPECT_EQ(0, qpp::modpow(0, 100, 5));
    EXPECT_EQ(0, qpp::modpow(100, 0, 1));
    EXPECT_EQ(1, qpp::modpow(100, 0, 5));
    EXPECT_EQ(0, qpp::modpow(12, 23, 1));
    EXPECT_EQ(0, qpp::modpow(2, 3, 4));
    EXPECT_EQ(34, qpp::modpow(17, 176, 37));
    EXPECT_EQ(4042, qpp::modpow(178373, 9281623, 6217));

    // test some really large numbers
    // assumes that qpp::bigint is 64 bits
    EXPECT_EQ(24502114, qpp::modpow(10000000019, 10000000019, 26527121));
    EXPECT_EQ(1847779, qpp::modpow(9000000019, 10000000119, 2652711));
    EXPECT_EQ(1099, qpp::modpow(9000000019, 10000000119, 1980));
    EXPECT_EQ(1, qpp::modpow(6897998630, 10000000018, 10000000019));
    EXPECT_EQ(0, qpp::modpow(maxbigint, maxbigint, maxbigint));
    EXPECT_EQ(1, qpp::modpow(maxbigint - 1, maxbigint - 1, maxbigint));
    EXPECT_EQ(32, qpp::modpow(maxbigint - 1, maxbigint - 2, maxbigint - 3));
}

TEST(qpp_modpow_exception, ParameterOutOfRange) {
    // example of how to unit-test for an exception
    EXPECT_THROW(qpp::modpow(0, 0, 2), qpp::exception::OutOfRange);
}
/******************************************************************************/
/// BEGIN bigint qpp::randprime(bigint a, bigint b, idx N = 1000)
TEST(qpp_randprime, AllTests) {
    EXPECT_TRUE(qpp::isprime(qpp::randprime(0, 100)));
    EXPECT_TRUE(qpp::isprime(qpp::randprime(100, 1000)));
    EXPECT_TRUE(qpp::isprime(qpp::randprime(10000, 10100)));
}
/******************************************************************************/
/// BEGIN std::vector<int> qpp::x2contfrac(double x, idx n, idx cut = 1e5)
TEST(qpp_x2contfrac, AllTests) {
    EXPECT_EQ(std::vector<int>({0}), qpp::x2contfrac(0, 3));
    EXPECT_EQ(std::vector<int>({1}), qpp::x2contfrac(1, 3));
    EXPECT_EQ(std::vector<int>({3, 7, 15, 1, 292, 1, 1, 1, 2, 1}),
              qpp::x2contfrac(qpp::pi, 10));
    EXPECT_EQ(std::vector<int>({0, -8, -9, -1, -135665, -1}),
              qpp::x2contfrac(-0.123456789, 6,
                              1e6)); // due to large term in c.f.e.
    EXPECT_EQ(std::vector<int>({0, -8, -9, -1}),
              qpp::x2contfrac(-0.123456789,
                              6)); // cuts the expansion at the 4th term
    EXPECT_EQ(std::vector<int>({0, -1, -80}),
              qpp::x2contfrac(-0.987654321, 10));
    EXPECT_EQ(std::vector<int>({-1, -4, -3, -1, -3, -1, -13565, -1, -8}),
              qpp::x2contfrac(-1.23456789, 9, 1e7));
}
/******************************************************************************/
