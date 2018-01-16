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

#include <algorithm>
#include "gtest/gtest.h"
#include "qpp.h"

using namespace qpp;

// Unit testing "random.h"

/******************************************************************************/
/// BEGIN inline bigint qpp::rand(bigint a, bigint b)
TEST(qpp_rand_integer, AllTests) {
    // 1 element equal boundaries
    bigint a = 42, b = a;
    EXPECT_EQ(a, qpp::rand(a, b));

    // 1000 elements between 0 and 1
    idx N = 1000;
    a = 0, b = 1;
    for (idx i = 0; i < N; ++i) {
        bigint n = qpp::rand(a, b);
        EXPECT_GE(n, a);
        EXPECT_LE(n, b);
    }

    // 10000 elements between -10 and 10
    N = 10000;
    a = -10, b = 10;
    double average = 0;
    for (idx i = 0; i < N; ++i) {
        bigint n = qpp::rand(a, b);

        EXPECT_GE(n, a);
        EXPECT_LE(n, b);
        average += n;
    }
    average /= N;
    EXPECT_NEAR(0, average, 2e-1); // very likely
}
/******************************************************************************/
/// BEGIN inline double qpp::rand(double a, double b)
TEST(qpp_rand_double, AllTests) {
    // 1000 elements between 0 and 1
    idx N = 1000;
    double a = 0, b = 1;
    for (idx i = 0; i < N; ++i) {
        double n = qpp::rand(a, b);
        EXPECT_GE(n, a);
        EXPECT_LE(n, b);
    }

    // 10000 elements between -10 and 10
    N = 10000;
    a = -10, b = 10;
    double average = 0;
    for (idx i = 0; i < N; ++i) {
        double n = qpp::rand(a, b);
        EXPECT_GE(n, a);
        EXPECT_LE(n, b);
        average += n;
    }
    average /= N;
    EXPECT_NEAR(0, average, 2e-1); // very likely
}
/******************************************************************************/
/// BEGIN template<> inline cmat qpp::rand(
///       idx rows, idx cols, double a, double b)
TEST(qpp_rand_cmat, AllTests) {
    // 1 x 1 matrix with elements between -1 and 1
    idx NA = 1, NB = 1;
    double a = -1, b = 1;
    cmat A = qpp::rand<cmat>(NA, NB, a, b);
    for (idx i = 0; i < NA * NB; ++i) {
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
    for (idx i = 0; i < NA * NB; ++i) {
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
    for (idx i = 0; i < NA * NB; ++i) {
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
TEST(qpp_rand_dmat, AllTests) {
    // 1 x 1 matrix with elements between -1 and 1
    idx NA = 1, NB = 1;
    double a = -1, b = 1;
    dmat A = qpp::rand<dmat>(NA, NB, a, b);
    for (idx i = 0; i < NA * NB; ++i) {
        double elem = A.data()[i];
        EXPECT_GE(elem, a);
        EXPECT_LT(elem, b);
    }

    // 1 x 10 matrix with elements between -10 and 10
    NA = 1, NB = 10;
    a = -10, b = 10;
    A = qpp::rand<dmat>(NA, NB, a, b);
    for (idx i = 0; i < NA * NB; ++i) {
        double elem = A.data()[i];
        EXPECT_GE(elem, a);
        EXPECT_LT(elem, b);
    }

    // 25 x 15 matrix with elements between -10 and 10
    NA = 25, NB = 15;
    a = -10, b = 10;
    A = qpp::rand<dmat>(NA, NB, a, b);
    double average = 0;
    for (idx i = 0; i < NA * NB; ++i) {
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
TEST(qpp_randH, AllTests) {
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
TEST(qpp_randidx, AllTests) {
    // 1 element equal boundaries
    idx a = 42, b = a;
    EXPECT_EQ(a, qpp::randidx(a, b));

    // 1000 elements between 0 and 1
    idx N = 1000;
    a = 0, b = 1;
    for (idx i = 0; i < N; ++i) {
        idx n = qpp::randidx(a, b);
        EXPECT_GE(n, a);
        EXPECT_LE(n, b);
    }

    // 10000 elements between -10 and 10
    N = 10000;
    a = 0, b = 10;
    double average = 0;
    for (idx i = 0; i < N; ++i) {
        idx n = qpp::randidx(a, b);
        EXPECT_GE(n, a);
        EXPECT_LE(n, b);
        average += n;
    }
    average /= N;
    EXPECT_NEAR(5, average, 1e-1); // very likely
}
/******************************************************************************/
/// BEGIN inline ket qpp::randket(idx D = 2)
TEST(qpp_randket, AllTests) {
    // D = 1 degenerate case
    idx D = 1;
    ket expected = ket::Ones(D);
    EXPECT_NEAR(0, norm(prj(expected) - prj(qpp::randket(D))), 1e-7);

    // D = 2, 100000 randomly generated kets, should average to zero
    D = 2;
    idx N = 100000; // number of runs
    ket avg_state = ket::Zero(D);
    for (idx i = 0; i < N; ++i) {
        avg_state += qpp::randket(D);
    }
    expected = ket::Zero(D);
    EXPECT_NEAR(0, norm(expected - avg_state / N), 2e-2);

    // D = 5, 100000 randomly generated kets, should average to zero
    D = 5;
    N = 100000; // number of runs
    avg_state = ket::Zero(D);
    for (idx i = 0; i < N; ++i) {
        avg_state += qpp::randket(D);
    }
    expected = ket::Zero(D);
    EXPECT_NEAR(0, norm(expected - avg_state / N), 2e-2);
}
/******************************************************************************/
/// BEGIN inline std::vector<cmat> qpp::randkraus(idx N, idx D = 2)
TEST(qpp_randkraus, AllTests) {
    // D = 1, N = 1 degenerate case
    idx D = 1, N = 1;
    std::vector<cmat> Ks = qpp::randkraus(N, D);
    cmat closure = cmat::Zero(D, D);
    for (idx i = 0; i < N; ++i) {
        closure += adjoint(Ks[i]) * Ks[i];
    }
    EXPECT_NEAR(0, norm(closure - gt.Id(D)), 1e-7);

    // D = 1, N = 10 degenerate case
    D = 1, N = 10;
    Ks = qpp::randkraus(N, D);
    closure = cmat::Zero(D, D);
    for (idx i = 0; i < N; ++i) {
        closure += adjoint(Ks[i]) * Ks[i];
    }
    EXPECT_NEAR(0, norm(closure - gt.Id(D)), 1e-7);

    // D = 2, N = 1
    D = 2, N = 1;
    Ks = qpp::randkraus(N);
    closure = cmat::Zero(D, D);
    for (idx i = 0; i < N; ++i) {
        closure += adjoint(Ks[i]) * Ks[i];
    }
    EXPECT_NEAR(0, norm(closure - gt.Id(D)), 1e-7);

    // D = 2, N = 10
    D = 2, N = 10;
    Ks = qpp::randkraus(N);
    closure = cmat::Zero(D, D);
    for (idx i = 0; i < N; ++i) {
        closure += adjoint(Ks[i]) * Ks[i];
    }
    EXPECT_NEAR(0, norm(closure - gt.Id(D)), 1e-7);

    // D = 5, N = 5
    D = 5, N = 10;
    Ks = qpp::randkraus(N, D);
    closure = cmat::Zero(D, D);
    for (idx i = 0; i < N; ++i) {
        closure += adjoint(Ks[i]) * Ks[i];
    }
    EXPECT_NEAR(0, norm(closure - gt.Id(D)), 1e-7);
}
/******************************************************************************/
/// BEGIN inline double qpp::randn(double mean = 0, double sigma = 1)
TEST(qpp_randn_double, AllTests) {
    // 10000 elements with mean = 0 and sigma = 1
    double mean = 0, sigma = 1;
    idx N = 10000;
    double average = 0;
    for (idx i = 0; i < N; ++i) {
        double n = qpp::randn(mean, sigma);
        average += n;
    }
    average /= N;
    EXPECT_NEAR(mean, average, 2.58 * sigma); // very likely (99%)

    // 10000 elements with mean = 1 and sigma = 2
    mean = 1, sigma = 2;
    N = 10000;
    average = 0;
    for (idx i = 0; i < N; ++i) {
        double n = qpp::randn(mean, sigma);
        average += n;
    }
    average /= N;
    EXPECT_NEAR(mean, average, 2.58 * sigma); // very likely (99%)

    // 10000 elements with mean = -1 and sigma = 1e-7
    mean = 1, sigma = 1e-7;
    N = 10000;
    average = 0;
    for (idx i = 0; i < N; ++i) {
        double n = qpp::randn(mean, sigma);
        average += n;
    }
    average /= N;
    EXPECT_NEAR(mean, average, 2.58 * sigma); // very likely (99%)
}
/******************************************************************************/
/// BEGIN template<> inline cmat qpp::randn(
///       idx rows, idx cols, double mean, double sigma)
TEST(qpp_randn_cmat, AllTests) {
    // 1 x 1 matrix with elements of mean = 0 and sigma = 0.1
    idx NA = 1, NB = 1;
    double mean = 0, sigma = 0.1;
    cmat A = qpp::randn<cmat>(NA, NB, mean, sigma);
    cplx average = 0;
    for (idx i = 0; i < NA * NB; ++i) {
        cplx elem = A.data()[i];
        average += elem;
    }
    average /= NA * NB;
    EXPECT_NEAR(mean, average.real(), 2.58 * sigma); // very likely (99%)
    EXPECT_NEAR(mean, average.imag(), 2.58 * sigma); // very likely (99%)

    // 1 x 10 matrix with elements of mean = -1 and sigma = 0.1
    NA = 1, NB = 10;
    mean = -1, sigma = 0.1;
    A = qpp::randn<cmat>(NA, NB, mean, sigma);
    average = 0;
    for (idx i = 0; i < NA * NB; ++i) {
        cplx elem = A.data()[i];
        average += elem;
    }
    average /= NA * NB;
    EXPECT_NEAR(mean, average.real(), 2.58 * sigma); // very likely (99%)
    EXPECT_NEAR(mean, average.imag(), 2.58 * sigma); // very likely (99%)

    // 25 x 15 matrix with elements of mean = 10 and sigma = 1e-7
    NA = 25, NB = 15;
    mean = 10, sigma = 1e-7;
    A = qpp::randn<cmat>(NA, NB, mean, sigma);
    average = 0;
    for (idx i = 0; i < NA * NB; ++i) {
        cplx elem = A.data()[i];
        average += elem;
    }
    average /= NA * NB;
    EXPECT_NEAR(mean, average.real(), 2.58 * sigma); // very likely (99%)
    EXPECT_NEAR(mean, average.imag(), 2.58 * sigma); // very likely (99%)
}
/******************************************************************************/
/// BEGIN template<> inline dmat qpp::randn(
///       idx rows, idx cols, double mean, double sigma)
TEST(qpp_randn_dmat, AllTests) {
    // 1 x 1 matrix with elements of mean = 0 and sigma = 0.1
    idx NA = 1, NB = 1;
    double mean = 0, sigma = 0.1;
    dmat A = qpp::randn<dmat>(NA, NB, mean, sigma);
    double average = 0;
    for (idx i = 0; i < NA * NB; ++i) {
        double elem = A.data()[i];
        average += elem;
    }
    average /= NA * NB;
    EXPECT_NEAR(mean, average, 2.58 * sigma); // very likely (99%)

    // 1 x 10 matrix with elements of mean = -1 and sigma = 0.1
    NA = 1, NB = 10;
    mean = -1, sigma = 0.1;
    A = qpp::randn<dmat>(NA, NB, mean, sigma);
    average = 0;
    for (idx i = 0; i < NA * NB; ++i) {
        double elem = A.data()[i];
        average += elem;
    }
    average /= NA * NB;
    EXPECT_NEAR(mean, average, 2.58 * sigma); // very likely (99%)

    // 25 x 15 matrix with elements of mean = 10 and sigma = 1e-7
    NA = 25, NB = 15;
    mean = 10, sigma = 1e-7;
    A = qpp::randn<dmat>(NA, NB, mean, sigma);
    average = 0;
    for (idx i = 0; i < NA * NB; ++i) {
        double elem = A.data()[i];
        average += elem;
    }
    average /= NA * NB;
    EXPECT_NEAR(mean, average, 2.58 * sigma); // very likely (99%)
}
/******************************************************************************/
/// BEGIN inline std::vector<idx> qpp::randperm(idx N)
TEST(qpp_randperm, AllTests) {
    idx N = 1;
    std::vector<idx> result = qpp::randperm(N);
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
TEST(qpp_randprob, AllTests) {
    idx N = 1;
    std::vector<double> result = qpp::randprob(N);
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
TEST(qpp_randrho, AllTests) {
    // D = 1 degenerate case
    idx D = 1;
    cmat rho = qpp::randrho(D);
    EXPECT_NEAR(0, norm(rho - adjoint(rho)), 1e-7);
    dyn_col_vect<double> lambdas = hevals(rho);
    for (idx i = 0; i < lambdas.size(); ++i) {
        EXPECT_GE(lambdas[i], 0);
        EXPECT_LE(lambdas[i], 1);
    }
    EXPECT_NEAR(1, trace(rho).real(), 1e-7);
    EXPECT_NEAR(0, trace(rho).imag(), 1e-7);

    // D = 2
    D = 2;
    rho = qpp::randrho(D);
    EXPECT_NEAR(0, norm(rho - adjoint(rho)), 1e-7);
    lambdas = hevals(rho);
    for (idx i = 0; i < lambdas.size(); ++i) {
        EXPECT_GE(lambdas[i], 0);
        EXPECT_LE(lambdas[i], 1);
    }
    EXPECT_NEAR(1, trace(rho).real(), 1e-7);
    EXPECT_NEAR(0, trace(rho).imag(), 1e-7);

    // D = 5
    D = 5;
    rho = qpp::randrho(D);
    EXPECT_NEAR(0, norm(rho - adjoint(rho)), 1e-7);
    lambdas = hevals(rho);
    for (idx i = 0; i < lambdas.size(); ++i) {
        EXPECT_GE(lambdas[i], 0);
        EXPECT_LE(lambdas[i], 1);
    }
    EXPECT_NEAR(1, trace(rho).real(), 1e-7);
    EXPECT_NEAR(0, trace(rho).imag(), 1e-7);
}
/******************************************************************************/
/// BEGIN inline cmat qpp::randU(idx D = 2)
TEST(qpp_randU, AllTests) {
    // D = 1 degenerate case
    idx D = 1;
    cmat U = qpp::randU(D);
    EXPECT_NEAR(0, norm(U * adjoint(U) - gt.Id(D)), 1e-7);

    // D = 2
    D = 2;
    U = qpp::randU();
    EXPECT_NEAR(0, norm(U * adjoint(U) - gt.Id(D)), 1e-7);

    idx N = 10000; // number of runs, average over random unitaries
    ket avg_state = ket::Zero(D);
    for (idx i = 0; i < N; ++i) {
        avg_state += qpp::randU() * st.z0;
    }
    ket expected = ket::Zero(D);
    EXPECT_NEAR(0, norm(expected - avg_state / N), 2e-2);

    // D = 10
    D = 10;
    U = qpp::randU(D);
    EXPECT_NEAR(0, norm(U * adjoint(U) - gt.Id(D)), 1e-7);
}
/******************************************************************************/
/// BEGIN inline cmat qpp::randV(idx Din, idx Dout)
TEST(qpp_randV, AllTests) {
    // Din = 1, Dout = 1 degenerate case
    idx Din = 1, Dout = 1;
    cmat V = qpp::randV(Din, Dout);
    EXPECT_NEAR(0, norm(adjoint(V) * V - gt.Id(Din)), 1e-7);

    // Din = 1, Dout = 10 degenerate case
    Din = 1, Dout = 10;
    V = qpp::randV(Din, Dout);
    EXPECT_NEAR(0, norm(adjoint(V) * V - gt.Id(Din)), 1e-7);

    // Din = 2, Dout = 2
    Din = 2, Dout = 2;
    V = qpp::randV(Din, Dout);
    EXPECT_NEAR(0, norm(adjoint(V) * V - gt.Id(Din)), 1e-7);

    // Din = 3, Dout = 5
    Din = 3, Dout = 5;
    V = qpp::randV(Din, Dout);
    EXPECT_NEAR(0, norm(adjoint(V) * V - gt.Id(Din)), 1e-7);
}
/******************************************************************************/
