#include <limits>
#include <tuple>
#include <vector>

#include "gtest/gtest.h"

#include "qpp.h"

using namespace qpp;

// Unit testing "number_theory.hpp"

/******************************************************************************/
/// BEGIN std::vector<idx> compperm(const std::vector<idx>& perm,
///       const std::vector<idx>& sigma)
TEST(qpp_compperm, AllTests) {
    std::vector<idx> sigma1{0};
    std::vector<idx> tau1{0};
    std::vector<idx> result1{0};
    EXPECT_EQ(result1, compperm(sigma1, tau1));

    std::vector<idx> sigma2{0, 1, 2};
    std::vector<idx> tau2{0, 1, 2};
    std::vector<idx> result2{0, 1, 2};
    EXPECT_EQ(result2, compperm(sigma2, tau2));

    std::vector<idx> sigma3{3, 0, 1, 4, 5, 2, 6};
    std::vector<idx> tau3{1, 2, 5, 0, 3, 4, 6}; // inverse permutation
    std::vector<idx> result3{0, 1, 2, 3, 4, 5, 6};
    EXPECT_EQ(result3, compperm(sigma3, tau3));
}
/******************************************************************************/
/// BEGIN realT contfrac2x(const std::vector<bigint>& cf, idx N = idx(-1))
TEST(qpp_contfrac2x, DefaultN) {
    EXPECT_NEAR(0, contfrac2x({0}), 1e-5);
    EXPECT_NEAR(42, contfrac2x({42}), 1e-5);
    EXPECT_NEAR(-42, contfrac2x({-42}), 1e-5);
    EXPECT_NEAR(pi, contfrac2x({3, 7, 15, 1, 292, 1}), 1e-5);
    EXPECT_NEAR(-pi, contfrac2x({-3, -7, -15, -1, -292, -1}), 1e-5);
    EXPECT_NEAR(0.1234, contfrac2x({0, 8, 9, 1, 1, 1}), 1e-5);
    EXPECT_NEAR(-0.4321, contfrac2x({0, -2, -3, -5, -2}), 1e-5);
    EXPECT_NEAR(std::sqrt(19), contfrac2x({4, 2, 1, 3, 1, 2, 8}), 1e-3);
}
/******************************************************************************/
TEST(qpp_contfrac2x, UserSpecifiedN) {
    EXPECT_NEAR(0, contfrac2x({0}, 1), 1e-5);
    EXPECT_NEAR(42, contfrac2x({42}, 1), 1e-5);
    EXPECT_NEAR(-42, contfrac2x({-42}, 1), 1e-5);
    EXPECT_NEAR(3.141, contfrac2x({3, 7, 15, 1, 292, 1}, 3), 1e-2);
    EXPECT_NEAR(-3.141, contfrac2x({-3, -7, -15, -1, -292, -1}, 3), 1e-2);
    EXPECT_NEAR(0.1234, contfrac2x({0, 8, 9, 1, 1, 1}, 4), 1e-3);
    EXPECT_NEAR(-0.4321, contfrac2x({0, -2, -3, -5, -2}, 4), 1e-3);
}
/******************************************************************************/
/// BEGIN std::vector<std::pair<bigint, bigint>>
///       convergents(const std::vector<bigint>& cf)
TEST(qpp_convergents, ContinuedFraction) {}
/******************************************************************************/
/// BEGIN std::vector<std::pair<bigint, bigint>>
///       convergents(realT x, idx N)
TEST(qpp_convergents, RealNumber) {}
/******************************************************************************/
/// BEGIN std::tuple<bigint, bigint, bigint> egcd(bigint m, bigint n)
TEST(qpp_egcd, NonNegativeNumbers) {
    EXPECT_EQ(std::make_tuple(1, 0, 10), egcd(10, 0));
    EXPECT_EQ(std::make_tuple(0, 1, 10), egcd(0, 10));
    EXPECT_EQ(std::make_tuple(0, 1, 1), egcd(1, 1));
    EXPECT_EQ(std::make_tuple(-103, 11, 4), egcd(120, 1124));
    EXPECT_EQ(std::make_tuple(54, -7, 1), egcd(17, 131));
}
/******************************************************************************/
TEST(qpp_egcd, MixedNumbers) {
    EXPECT_EQ(std::make_tuple(-1, 0, 10), egcd(-10, 0));
    EXPECT_EQ(std::make_tuple(0, -1, 10), egcd(0, -10));
    EXPECT_EQ(std::make_tuple(0, -1, 1), egcd(1, -1));
    EXPECT_EQ(std::make_tuple(0, 1, 1), egcd(-1, 1));
    EXPECT_EQ(std::make_tuple(-103, -11, 4), egcd(120, -1124));
    EXPECT_EQ(std::make_tuple(-54, -7, 1), egcd(-17, 131));
}
/******************************************************************************/
TEST(qpp_egcd, NegativeNumbers) {
    EXPECT_EQ(std::make_tuple(0, -1, 1), egcd(-1, -1));
    EXPECT_EQ(std::make_tuple(103, -11, 4), egcd(-120, -1124));
    EXPECT_EQ(std::make_tuple(-54, 7, 1), egcd(-17, -131));
}
/******************************************************************************/
/// BEGIN std::vector<bigint> factors(bigint n)
TEST(qpp_factors, AllTests) {
    bigint n1 = 2;
    std::vector<bigint> result1{2};
    EXPECT_EQ(result1, factors(n1));

    bigint n2 = -10;
    std::vector<bigint> result2{2, 5};
    EXPECT_EQ(result2, factors(n2));

    bigint n3 = 110500;
    std::vector<bigint> result3{2, 2, 5, 5, 5, 13, 17};
    EXPECT_EQ(result3, factors(n3));

    bigint n4 = -35750;
    std::vector<bigint> result4{2, 5, 5, 5, 11, 13};
    EXPECT_EQ(result4, factors(n4));

    // test some really large numbers
    bigint n5 = 1000000000001112;
    std::vector<bigint> result5{2, 2, 2, 3, 11, 29, 1217, 1637, 65563};
    EXPECT_EQ(result5, factors(n5));

    bigint n6 = 1000000000001118;
    std::vector<bigint> result6{2, 3, 897649, 185670197};
    EXPECT_EQ(result6, factors(n6));
}
/******************************************************************************/
/// BEGIN bigint gcd(bigint m, bigint n)
TEST(qpp_gcd, NonNegativeNumbers) {
    EXPECT_EQ(10, gcd(10, 0));
    EXPECT_EQ(10, gcd(0, 10));
    EXPECT_EQ(1, gcd(1, 1));
    EXPECT_EQ(4, gcd(120, 1124));
    EXPECT_EQ(1, gcd(17, 131));

    // test some really large numbers
    // assumes that bigint is 64 bits
    bigint maxbigint = std::numeric_limits<bigint>::max();
    EXPECT_EQ(maxbigint, gcd(maxbigint, maxbigint));
    EXPECT_EQ(1, gcd(maxbigint, maxbigint - 1));
    EXPECT_EQ(3, gcd(maxbigint - 1, maxbigint - 10));
}
/******************************************************************************/
TEST(qpp_gcd, MixedNumbers) {
    EXPECT_EQ(10, gcd(-10, 0));
    EXPECT_EQ(10, gcd(0, -10));
    EXPECT_EQ(1, gcd(1, -1));
    EXPECT_EQ(1, gcd(-1, 1));
    EXPECT_EQ(4, gcd(120, -1124));
    EXPECT_EQ(1, gcd(-17, 131));
}
/******************************************************************************/
TEST(qpp_gcd, NegativeNumbers) {
    EXPECT_EQ(1, gcd(-1, -1));
    EXPECT_EQ(4, gcd(-120, -1124));
    EXPECT_EQ(1, gcd(-17, -131));
    EXPECT_EQ(20, gcd(-170180, -13100));
}
/******************************************************************************/
/// BEGIN bigint gcd(const std::vector<bigint>& ns)
TEST(qpp_gcd, ListMixedNumbers) {
    std::vector<bigint> v1{1, 0, -1, 0, 1};
    EXPECT_EQ(1, gcd(v1));

    std::vector<bigint> v2{1, 2, 3, 4, 5};
    EXPECT_EQ(1, gcd(v2));

    std::vector<bigint> v3{125000, -12680, 25040, -4};
    EXPECT_EQ(4, gcd(v3));
}
/******************************************************************************/
/// BEGIN std::vector<idx> invperm(const std::vector<idx>& perm)
TEST(qpp_invperm, AllTests) {
    std::vector<idx> perm1{0};
    std::vector<idx> result1{0};
    EXPECT_EQ(result1, invperm(perm1));

    std::vector<idx> perm2{0, 1, 2, 3, 4, 5, 6};
    std::vector<idx> result2{0, 1, 2, 3, 4, 5, 6};
    EXPECT_EQ(result2, invperm(perm2));

    std::vector<idx> perm3{6, 5, 4, 3, 2, 1, 0};
    std::vector<idx> result3{6, 5, 4, 3, 2, 1, 0};
    EXPECT_EQ(result3, invperm(perm3));

    std::vector<idx> perm4{3, 0, 1, 4, 5, 2, 6};
    std::vector<idx> result4{1, 2, 5, 0, 3, 4, 6};
    EXPECT_EQ(result4, invperm(perm4));
}
/******************************************************************************/
/// BEGIN bool isprime(bigint n, idx k = 80)
TEST(qpp_isprime, MixedPrimeNumbers) {
    EXPECT_TRUE(isprime(2));
    EXPECT_TRUE(isprime(-2));
    EXPECT_TRUE(isprime(11));
    EXPECT_TRUE(isprime(-17));
    EXPECT_TRUE(isprime(19));
    EXPECT_TRUE(isprime(127));
    EXPECT_TRUE(isprime(541));
    EXPECT_TRUE(isprime(-104729));
    EXPECT_TRUE(isprime(10000000019));
}
/******************************************************************************/
TEST(qpp_isprime, MixedNonPrimeNumbers) {
    EXPECT_FALSE(isprime(4));
    EXPECT_FALSE(isprime(-4));
    EXPECT_FALSE(isprime(110));
    EXPECT_FALSE(isprime(-4891));
    EXPECT_FALSE(isprime(-13101));
    EXPECT_FALSE(isprime(110011));
    EXPECT_FALSE(isprime(10000000119));

    // test some Charmichael numbers
    EXPECT_FALSE(isprime(561));
    EXPECT_FALSE(isprime(6601));
    EXPECT_FALSE(isprime(8911));
    EXPECT_FALSE(isprime(41041));
}
/******************************************************************************/
/// BEGIN bigint lcm(bigint m, bigint n)
TEST(qpp_lcm, NonNegativeNumbers) {
    EXPECT_EQ(0, lcm(10, 0));
    EXPECT_EQ(0, lcm(0, 10));
    EXPECT_EQ(1, lcm(1, 1));
    EXPECT_EQ(33720, lcm(120, 1124));
    EXPECT_EQ(2227, lcm(17, 131));
}
/******************************************************************************/
TEST(qpp_lcm, MixedNumbers) {
    EXPECT_EQ(0, lcm(-10, 0));
    EXPECT_EQ(0, lcm(0, -10));
    EXPECT_EQ(1, lcm(1, -1));
    EXPECT_EQ(1, lcm(-1, 1));
    EXPECT_EQ(33720, lcm(120, -1124));
    EXPECT_EQ(2227, lcm(-17, 131));
}
/******************************************************************************/
TEST(qpp_lcm, NegativeNumbers) {
    EXPECT_EQ(1, lcm(-1, -1));
    EXPECT_EQ(33720, lcm(-120, -1124));
    EXPECT_EQ(2227, lcm(-17, -131));
    EXPECT_EQ(111467900, lcm(-170180, -13100));
}
/******************************************************************************/
/// BEGIN bigint lcm(const std::vector<bigint>& ns)
TEST(qpp_lcm, ListMixedNumbers) {
    std::vector<bigint> v1{10, -10, 10, -10};
    EXPECT_EQ(10, lcm(v1));

    std::vector<bigint> v2{1, 2, 3, 4, 5};
    EXPECT_EQ(60, lcm(v2));

    std::vector<bigint> v3{12500, -1268, 2504, -4};
    EXPECT_EQ(2480525000, lcm(v3));
}
/******************************************************************************/
/// BEGIN bigint modinv(bigint a, bigint p)
TEST(qpp_modinv, AllTests) {
    EXPECT_EQ(1, modinv(1, 1));
    EXPECT_EQ(62, modinv(2, 123));
    EXPECT_EQ(21, modinv(1231, 22));
}
/******************************************************************************/
/// BEGIN bigint modmul(bigint a, bigint n, bigint p)
TEST(qpp_modmul, NonNegativeNumbers) {
    EXPECT_EQ(0, modmul(0, 0, 1));
    EXPECT_EQ(0, modmul(0, 0, 2));
    EXPECT_EQ(0, modmul(1, 0, 2));
    EXPECT_EQ(0, modmul(0, 1, 2));
    EXPECT_EQ(2611, modmul(12127, 71623, 12345));

    // test some really large numbers
    // assumes that bigint is 64 bits
    bigint maxbigint = std::numeric_limits<bigint>::max();

    EXPECT_EQ(42, modmul(maxbigint - 1, maxbigint, 123));
    EXPECT_EQ(49, modmul(maxbigint, maxbigint, 123));
    EXPECT_EQ(2262, modmul(maxbigint - 189, maxbigint - 2345, 7891));
    EXPECT_EQ(0, modmul(maxbigint, maxbigint, maxbigint));
    EXPECT_EQ(14884, modmul(maxbigint - 1, maxbigint - 1, maxbigint - 123));
}
/******************************************************************************/
TEST(qpp_modmul, MixedNumbers) {
    EXPECT_EQ(0, modmul(-1, 0, 2));
    EXPECT_EQ(0, modmul(0, -1, 2));
    EXPECT_EQ(9734, modmul(12127, -71623, 12345));
    EXPECT_EQ(9734, modmul(-12127, 71623, 12345));

    // test some really large numbers
    // assumes that bigint is 64 bits
    bigint minbigint = std::numeric_limits<bigint>::min();
    bigint maxbigint = std::numeric_limits<bigint>::max();

    EXPECT_EQ(1114, modmul(minbigint, maxbigint, 2314));
    EXPECT_EQ(21, modmul(-maxbigint, maxbigint, 34));
    EXPECT_EQ(240, modmul(maxbigint, -1234567890, 314));
    EXPECT_EQ(219, modmul(maxbigint, (maxbigint - 1) / 2, 1314));
}
/******************************************************************************/
TEST(qpp_modmul, NegativeNumbers) {
    // test some really large numbers
    // assumes that bigint is 64 bits
    bigint minbigint = std::numeric_limits<bigint>::min();
    bigint maxbigint = std::numeric_limits<bigint>::max();
    EXPECT_EQ(13, modmul(-maxbigint, -maxbigint, 34));
    EXPECT_EQ(1, modmul(-maxbigint, -maxbigint, maxbigint - 1));
    EXPECT_EQ(74, modmul(-maxbigint, -1234567890, 314));
    EXPECT_EQ(1799, modmul(minbigint + 1234, -maxbigint + 2345, 7891));

    EXPECT_EQ(64, modmul(minbigint, minbigint, 123));
    EXPECT_EQ(56, modmul(minbigint + 1, minbigint, 123));
    EXPECT_EQ(1799, modmul(minbigint + 1234, -maxbigint + 2345, 7891));
}
/******************************************************************************/
/// BEGIN bigint modpow(bigint a, bigint n, bigint p)
TEST(qpp_modpow, PositiveNumbers) {
    bigint maxbigint = std::numeric_limits<bigint>::max();

    EXPECT_EQ(0, modpow(0, 100, 1));
    EXPECT_EQ(0, modpow(0, 100, 5));
    EXPECT_EQ(0, modpow(100, 0, 1));
    EXPECT_EQ(1, modpow(100, 0, 5));
    EXPECT_EQ(0, modpow(12, 23, 1));
    EXPECT_EQ(0, modpow(2, 3, 4));
    EXPECT_EQ(34, modpow(17, 176, 37));
    EXPECT_EQ(4042, modpow(178373, 9281623, 6217));

    // test some really large numbers
    // assumes that bigint is 64 bits
    EXPECT_EQ(24502114, modpow(10000000019, 10000000019, 26527121));
    EXPECT_EQ(1847779, modpow(9000000019, 10000000119, 2652711));
    EXPECT_EQ(1099, modpow(9000000019, 10000000119, 1980));
    EXPECT_EQ(1, modpow(6897998630, 10000000018, 10000000019));
    EXPECT_EQ(0, modpow(maxbigint, maxbigint, maxbigint));
    EXPECT_EQ(1, modpow(maxbigint - 1, maxbigint - 1, maxbigint));
    EXPECT_EQ(32, modpow(maxbigint - 1, maxbigint - 2, maxbigint - 3));
}
/******************************************************************************/
TEST(qpp_modpow, ExceptionParameterOutOfRange) {
    // example of how to unit-test for an exception
    EXPECT_THROW(modpow(0, 0, 2), exception::OutOfRange);
}
/******************************************************************************/
/// BEGIN bigint randprime(bigint a, bigint b, idx N = 1000)
TEST(qpp_randprime, AllTests) {
    EXPECT_TRUE(isprime(randprime(0, 100)));
    EXPECT_TRUE(isprime(randprime(100, 1000)));
    EXPECT_TRUE(isprime(randprime(10000, 10100)));
}
/******************************************************************************/
/// BEGIN std::vector<bigint> x2contfrac(realT x, idx n, idx cut = 1e5)
TEST(qpp_x2contfrac, AllTests) {
    EXPECT_EQ(std::vector<bigint>({0}), x2contfrac(0, 3));
    EXPECT_EQ(std::vector<bigint>({1}), x2contfrac(1, 3));
    EXPECT_EQ(std::vector<bigint>({3, 7, 15, 1, 292, 1, 1, 1, 2, 1}),
              x2contfrac(pi, 10));

    EXPECT_EQ(std::vector<bigint>({0, -8, -9, -1, -135665, -1}),
              x2contfrac(-0.123456789, 6,
                         1000000)); // due to large term in c.f.e.

    EXPECT_EQ(
        std::vector<bigint>({0, -8, -9, -1}),
        x2contfrac(-0.123456789, 6)); // cuts the expansion at the 4th term
    EXPECT_EQ(std::vector<bigint>({0, -1, -80}), x2contfrac(-0.987654321, 10));
    EXPECT_EQ(std::vector<bigint>({-1, -4, -3, -1, -3, -1, -13565, -1, -8}),
              x2contfrac(-1.23456789, 9, 10000000));
}
/******************************************************************************/
