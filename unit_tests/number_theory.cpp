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

///// BEGIN std::vector<idx> qpp::compperm(const std::vector<idx>& perm,
/////       const std::vector<idx>& sigma)
TEST(qpp_compperm_test, NonNegativeNumbers)
{
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
///// END std::vector<idx> qpp::compperm(const std::vector<idx>& perm,
/////     const std::vector<idx>& sigma)

///// BEGIN double qpp::contfrac2x(const std::vector<int>& cf)
TEST(qpp_contfrac2x_test, AllTests)
{

}
///// END double qpp::contfrac2x(const std::vector<int>& cf)

///// BEGIN double qpp::contfrac2x(const std::vector<int>& cf, idx n)
TEST(qpp_contfrac2x_n_test, AllTests)
{

}
///// END double qpp::contfrac2x(const std::vector<int>& cf, idx n)

///// BEGIN std::tuple<bigint, bigint, bigint> qpp::egcd(bigint m, bigint n)
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
///// END std::tuple<bigint, bigint, bigint> qpp::egcd(bigint m, bigint n)

///// BEGIN std::vector<bigint> qpp::factors(bigint n)
TEST(qpp_factors_test, MixedNumbers)
{
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
    std::vector<bigint> result4{2, 5, 5, 5,11, 13};
    EXPECT_EQ(result4, qpp::factors(n4));
}
///// END std::vector<bigint> qpp::factors(bigint n)

///// BEGIN bigint qpp::gcd(bigint m, bigint n)
TEST(qpp_gcd_test, NonNegativeNumbers)
{
    EXPECT_EQ (10, qpp::gcd(10, 0));
    EXPECT_EQ (10, qpp::gcd(0, 10));
    EXPECT_EQ (1, qpp::gcd(1, 1));
    EXPECT_EQ (4, qpp::gcd(120, 1124));
    EXPECT_EQ (1, qpp::gcd(17, 131));
}

TEST(qpp_gcd_test, MixedNumbers)
{
    EXPECT_EQ (10, qpp::gcd(-10, 0));
    EXPECT_EQ (10, qpp::gcd(0, -10));
    EXPECT_EQ (1, qpp::gcd(1, -1));
    EXPECT_EQ (1, qpp::gcd(-1, 1));
    EXPECT_EQ (4, qpp::gcd(120, -1124));
    EXPECT_EQ (1, qpp::gcd(-17, 131));
}

TEST(qpp_gcd_test, NegativeNumbers)
{
    EXPECT_EQ (1, qpp::gcd(-1, -1));
    EXPECT_EQ (4, qpp::gcd(-120, -1124));
    EXPECT_EQ (1, qpp::gcd(-17, -131));
    EXPECT_EQ (20, qpp::gcd(-170180, -13100));
}
///// END bigint qpp::gcd(bigint m, bigint n)

///// BEGIN bigint qpp::gcd(const std::vector<bigint>& ns)
TEST(qpp_gcd_list_test, MixedNumbers)
{
    std::vector<bigint> v1{1, 0, -1, 0, 1};
    EXPECT_EQ(1, qpp::gcd(v1));

    std::vector<bigint> v2{1, 2, 3, 4, 5};
    EXPECT_EQ(1, qpp::gcd(v2));

    std::vector<bigint> v3{125000, -12680, 25040, -4};
    EXPECT_EQ(4, qpp::gcd(v3));
}
///// END bigint qpp::gcd(const std::vector<bigint>& ns)

///// BEGIN std::vector<idx> qpp::invperm(const std::vector<idx>& perm)
TEST(qpp_invperm_test, NonNegativeNumbers)
{
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
///// END std::vector<idx> qpp::invperm(const std::vector<idx>& perm)

///// BEGIN bool qpp::isprime(bigint n)
TEST(qpp_isprime_test, MixedPrimeNumbers)
{
    EXPECT_EQ(true, qpp::isprime(2));
    EXPECT_EQ(true, qpp::isprime(-2));
    EXPECT_EQ(true, qpp::isprime(11));
    EXPECT_EQ(true, qpp::isprime(541));
    EXPECT_EQ(true, qpp::isprime(-104729));
}

TEST(qpp_isprime_test, MixedNonPrimeNumbers)
{
    EXPECT_EQ(false, qpp::isprime(4));
    EXPECT_EQ(false, qpp::isprime(-4));
    EXPECT_EQ(false, qpp::isprime(110));
    EXPECT_EQ(false, qpp::isprime(110011));
    EXPECT_EQ(false, qpp::isprime(-13101));
}
///// END bool qpp::isprime(bigint n)

///// BEGIN bigint qpp::lcm(bigint m, bigint n)
TEST(qpp_lcm_test, NonNegativeNumbers)
{
    EXPECT_EQ (0, qpp::lcm(10, 0));
    EXPECT_EQ (0, qpp::lcm(0, 10));
    EXPECT_EQ (1, qpp::lcm(1, 1));
    EXPECT_EQ (33720, qpp::lcm(120, 1124));
    EXPECT_EQ (2227, qpp::lcm(17, 131));
}

TEST(qpp_lcm_test, MixedNumbers)
{
    EXPECT_EQ (0, qpp::lcm(-10, 0));
    EXPECT_EQ (0, qpp::lcm(0, -10));
    EXPECT_EQ (1, qpp::lcm(1, -1));
    EXPECT_EQ (1, qpp::lcm(-1, 1));
    EXPECT_EQ (33720, qpp::lcm(120, -1124));
    EXPECT_EQ (2227, qpp::lcm(-17, 131));
}

TEST(qpp_lcm_test, NegativeNumbers)
{
    EXPECT_EQ (1, qpp::lcm(-1, -1));
    EXPECT_EQ (33720, qpp::lcm(-120, -1124));
    EXPECT_EQ (2227, qpp::lcm(-17, -131));
    EXPECT_EQ (111467900, qpp::lcm(-170180, -13100));
}
///// END bigint qpp::lcm(bigint m, bigint n)

///// BEGIN bigint qpp::lcm(const std::vector<bigint>& ns)
TEST(qpp_lcm_list_test, MixedNumbers)
{
    std::vector<bigint> v1{10, -10, 10, -10};
    EXPECT_EQ(10, qpp::lcm(v1));

    std::vector<bigint> v2{1, 2, 3, 4, 5};
    EXPECT_EQ(60, qpp::lcm(v2));

    std::vector<bigint> v3{12500, -1268, 2504, -4};
    EXPECT_EQ(2480525000, qpp::lcm(v3));
}
///// END bigint qpp::lcm(const std::vector<bigint>& ns)

///// BEGIN bigint qpp::modinv(bigint a, bigint p)
TEST(qpp_modinv_test, NonNegativeNumbers)
{
    EXPECT_EQ (1, qpp::modinv(1, 1));
    EXPECT_EQ (62, qpp::modinv(2, 123));
    EXPECT_EQ (21, qpp::modinv(1231, 22));
}
///// END bigint qpp::modinv(bigint a, bigint p)

///// BEGIN bigint qpp::modpow(bigint a, bigint n, bigint p)
TEST(qpp_modpow_test, PositiveNumbers)
{
    EXPECT_EQ (0, qpp::modpow(0, 100, 1));
    EXPECT_EQ (0, qpp::modpow(0, 100, 5));
    EXPECT_EQ (0, qpp::modpow(100, 0, 1));
    EXPECT_EQ (1, qpp::modpow(100, 0, 5));
    EXPECT_EQ (0, qpp::modpow(2, 3, 4));
    EXPECT_EQ (34, qpp::modpow(17, 176, 37));
    EXPECT_EQ (4042, qpp::modpow(178373, 9281623, 6217));
}

TEST(qpp_modpow_exception_test, ParameterOutOfRange)
{
    // example of how to unit-test for an exception
    try
    {
        qpp::modpow(0, 0, 2); // 0^0 is not defined
        FAIL() << "Exception: qpp::modpow() ";
    }
    catch(const qpp::Exception& err)
    {
        std::string err_what{"IN qpp::modpow(): Parameter out of range!"};
        EXPECT_EQ(err.what(), err_what);
    }
    catch(...) {
        FAIL() << "Other exception";
    }
}
///// END bigint qpp::modpow(bigint a, bigint n, bigint p)

///// BEGIN std::vector<int> qpp::x2contfrac(double x, idx n, idx cut = 1e5)
TEST(qpp_x2contfrac_test, AllTests)
{

}
///// END std::vector<int> qpp::x2contfrac(double x, idx n, idx cut = 1e5)