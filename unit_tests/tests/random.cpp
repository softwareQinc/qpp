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

#include <algorithm>
#include "gtest/gtest.h"
#include "qpp.h"

using namespace qpp;

// Unit testing "random.h"

/******************************************************************************/
/// BEGIN inline bigint qpp::rand(bigint a, bigint b)
TEST(qpp_rand_integer, AllTests)
{
    // 1 element equal boundaries
    bigint a = 42, b = a;
    EXPECT_EQ(a, qpp::rand(a, b));

    // 1000 elements between 0 and 1
    idx N = 1000;
    a = 0, b = 1;
    for (idx i = 0; i < N; ++i)
    {
        auto n = qpp::rand(a, b);
        EXPECT_GE(n, a);
        EXPECT_LE(n, b);
    }

    // 10000 elements between -10 and 10
    N = 10000;
    a = -10, b = 10;
    double average = 0;
    for (idx i = 0; i < N; ++i)
    {
        auto n = qpp::rand(a, b);

        EXPECT_GE(n, a);
        EXPECT_LE(n, b);
        average += n;
    }
    average /= N;
    EXPECT_NEAR(0, average, 0.2); // very likely
}
/******************************************************************************/
/// BEGIN inline double qpp::rand(double a, double b)
TEST(qpp_rand_double, AllTests)
{
    // 1000 elements between 0 and 1
    idx N = 1000;
    double a = 0, b = 1;
    for (idx i = 0; i < N; ++i)
    {
        auto n = qpp::rand(a, b);
        EXPECT_GE(n, a);
        EXPECT_LE(n, b);
    }

    // 10000 elements between -10 and 10
    N = 10000;
    a = -10, b = 10;
    double average = 0;
    for (idx i = 0; i < N; ++i)
    {
        auto n = qpp::rand(a, b);
        EXPECT_GE(n, a);
        EXPECT_LE(n, b);
        average += n;
    }
    average /= N;
    EXPECT_NEAR(0, average, 0.2); // very likely
}
/******************************************************************************/
/// BEGIN template<> inline cmat qpp::rand(
///       idx rows, idx cols, double a, double b)
TEST(qpp_rand_cmat, AllTests)
{
    // 1 x 1 matrix with elements between -1 and 1
    idx NA = 1, NB = 1;
    double a = -1, b = 1;
    cmat A = qpp::rand<cmat>(NA, NB, a, b);
    for (idx i = 0; i < NA * NB; ++i)
    {
        cplx elem = A.data()[i];
        EXPECT_GE(elem.real(), a);
        EXPECT_GE(elem.imag(), a);
        EXPECT_LT(elem.real(), b);
        EXPECT_LT(elem.imag(), b);
    }

    // 1 x 10 matrix with elements between -10 and 10
    NA = 1, NB = 10;
    a = -10, b = 10;
    A = qpp::rand<cmat>(NA, NB, a, b);
    for (idx i = 0; i < NA * NB; ++i)
    {
        cplx elem = A.data()[i];
        EXPECT_GE(elem.real(), a);
        EXPECT_GE(elem.imag(), a);
        EXPECT_LT(elem.real(), b);
        EXPECT_LT(elem.imag(), b);
    }

    // 25 x 15 matrix with elements between -10 and 10
    NA = 25, NB = 15;
    a = -10, b = 10;
    A = qpp::rand<cmat>(NA, NB, a, b);
    cplx average = 0;
    for (idx i = 0; i < NA * NB; ++i)
    {
        cplx elem = A.data()[i];
        EXPECT_GE(elem.real(), a);
        EXPECT_GE(elem.imag(), a);
        EXPECT_LT(elem.real(), b);
        EXPECT_LT(elem.imag(), b);
        average += elem;
    }
    average /= (NA * NB);
    EXPECT_NEAR(0, average.real(), 1); // very likely
    EXPECT_NEAR(0, average.imag(), 1); // very likely
}
/******************************************************************************/
/// BEGIN template<> inline dmat qpp::rand(
///       idx rows, idx cols, double a, double b)
TEST(qpp_rand_dmat, AllTests)
{
    // 1 x 1 matrix with elements between -1 and 1
    idx NA = 1, NB = 1;
    double a = -1, b = 1;
    dmat A = qpp::rand<dmat>(NA, NB, a, b);
    for (idx i = 0; i < NA * NB; ++i)
    {
        double elem = A.data()[i];
        EXPECT_GE(elem, a);
        EXPECT_LT(elem, b);
    }

    // 1 x 10 matrix with elements between -10 and 10
    NA = 1, NB = 10;
    a = -10, b = 10;
    A = qpp::rand<dmat>(NA, NB, a, b);
    for (idx i = 0; i < NA * NB; ++i)
    {
        double elem = A.data()[i];
        EXPECT_GE(elem, a);
        EXPECT_LT(elem, b);
    }

    // 25 x 15 matrix with elements between -10 and 10
    NA = 25, NB = 15;
    a = -10, b = 10;
    A = qpp::rand<dmat>(NA, NB, a, b);
    double average = 0;
    for (idx i = 0; i < NA * NB; ++i)
    {
        double elem = A.data()[i];
        EXPECT_GE(elem, a);
        EXPECT_LT(elem, b);
        average += elem;
    }
    average /= (NA * NB);
    EXPECT_NEAR(0, average, 1); // very likely
}
/******************************************************************************/
/// BEGIN inline cmat qpp::randH(idx D = 2)
TEST(qpp_randH, AllTests)
{
    // D = 1 degenerate case
    idx D = 1;
    cmat A = randH(D);
    EXPECT_NEAR(0, norm(A - adjoint(A)), 1e-7);

    // D = 2
    D = 2;
    A = randH(D);
    EXPECT_NEAR(0, norm(A - adjoint(A)), 1e-7);

    // D = 10
    D = 10;
    A = randH(D);
    EXPECT_NEAR(0, norm(A - adjoint(A)), 1e-7);
}
/******************************************************************************/
/// BEGIN inline idx qpp::randidx(
///       idx a = std::numeric_limits<idx>::min(),
///       idx b = std::numeric_limits<idx>::max())
TEST(qpp_randidx, AllTests)
{
    // 1 element equal boundaries
    idx a = 42, b = a;
    EXPECT_EQ(a, qpp::randidx(a, b));

    // 1000 elements between 0 and 1
    idx N = 1000;
    a = 0, b = 1;
    for (idx i = 0; i < N; ++i)
    {
        auto n = qpp::randidx(a, b);
        EXPECT_GE(n, a);
        EXPECT_LE(n, b);
    }

    // 10000 elements between -10 and 10
    N = 10000;
    a = 0, b = 10;
    double average = 0;
    for (idx i = 0; i < N; ++i)
    {
        auto n = qpp::randidx(a, b);
        EXPECT_GE(n, a);
        EXPECT_LE(n, b);
        average += n;
    }
    average /= N;
    EXPECT_NEAR(5, average, 0.2); // very likely
}
/******************************************************************************/
/// BEGIN inline ket qpp::randket(idx D = 2)
TEST(qpp_randket, AllTests)
{
    // D = 1 degenerate case
    idx D = 1;
    ket expected = ket::Ones(D);
    EXPECT_NEAR(0, norm(expected) - norm(qpp::randket(D)), 1e-7);


}
/******************************************************************************/
/// BEGIN inline std::vector<cmat> qpp::randkraus(idx N, idx D = 2)
TEST(qpp_randkraus, AllTests)
{

}
/******************************************************************************/
/// BEGIN inline double qpp::randn(double mean = 0, double sigma = 1)
TEST(qpp_randn_double, AllTests)
{

}
/******************************************************************************/
/// BEGIN template<> inline cmat qpp::randn(
///       idx rows, idx cols, double mean, double sigma)
TEST(qpp_randn_cmat, AllTests)
{

}
/******************************************************************************/
/// BEGIN template<> inline dmat qpp::randn(
///       idx rows, idx cols, double mean, double sigma)
TEST(qpp_randn_dmat, AllTests)
{

}
/******************************************************************************/
/// BEGIN inline std::vector<idx> qpp::randperm(idx N)
TEST(qpp_randperm, AllTests)
{
    idx N = 1;
    auto result = qpp::randperm(N);
    std::sort(std::begin(result), std::end(result));
    std::vector<idx> expected(N);
    std::iota(std::begin(expected), std::end(expected), 0u);
    EXPECT_TRUE(result == expected);

    N = 2;
    result = qpp::randperm(N);
    std::sort(std::begin(result), std::end(result));
    expected.resize(N);
    std::iota(std::begin(expected), std::end(expected), 0u);
    EXPECT_TRUE(result == expected);

    N = 20;
    result = qpp::randperm(N);
    expected.resize(N);
    std::iota(std::begin(expected), std::end(expected), 0u);
    EXPECT_FALSE(result == expected); // very very likely
    std::sort(std::begin(result), std::end(result));
    EXPECT_TRUE(result == expected);
}
/******************************************************************************/
/// BEGIN inline std::vector<double> qpp::randprob(idx N)
TEST(qpp_randprob, AllTests)
{
    idx N = 1;
    auto result = qpp::randprob(N);
    EXPECT_EQ(1, result.size());
    EXPECT_NEAR(1, sum(result), 1e-7);

    N = 2;
    result = qpp::randprob(N);
    for (idx i = 0; i < N; ++i)
        EXPECT_GE(result[i], 0);
    EXPECT_EQ(2, result.size());
    EXPECT_NEAR(1, sum(result), 1e-7);

    N = 10;
    result = qpp::randprob(N);
    for (idx i = 0; i < N; ++i)
        EXPECT_GE(result[i], 0);
    EXPECT_EQ(10, result.size());
    EXPECT_NEAR(1, sum(result), 1e-7);
}
/******************************************************************************/
/// BEGIN inline cmat qpp::randrho(idx D = 2)
TEST(qpp_randrho, AllTests)
{

}
/******************************************************************************/
/// BEGIN inline cmat qpp::randU(idx D = 2)
TEST(qpp_randU, AllTests)
{

}
/******************************************************************************/
/// BEGIN inline cmat qpp::randV(idx Din, idx Dout)
TEST(qpp_randV, AllTests)
{

}
/******************************************************************************/
