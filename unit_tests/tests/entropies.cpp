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
#include "gtest/gtest.h"
#include "qpp.h"

using namespace qpp;

// Unit testing "entropies.h"

/******************************************************************************/
/// BEGIN template<typename Derived> double qpp::entropy(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_entropy_matrix, AllTests) {
    // 1 x 1 case
    cmat A(1, 1);
    A << 1.;
    EXPECT_NEAR(0, qpp::entropy(A), 1e-7);

    // 2 x 2 random matrix with fixed only 1 non-zero eigenvalue
    idx D = 2;
    cmat evals = cmat::Zero(D, 1);
    evals << 1., 0.;
    A = evals.asDiagonal();
    cmat U = randU(D);
    EXPECT_NEAR(0, qpp::entropy(A), 1e-7);

    // 2 x 2 random matrix with fixed eigenvalues
    D = 2;
    evals = cmat::Zero(D, 1);
    evals << 0.6, 0.4;
    A = evals.asDiagonal();
    U = randU(D);
    EXPECT_NEAR(0.970950594455, qpp::entropy(A), 1e-7);

    // 2 x 2 random matrix with fixed equal eigenvalues
    D = 2;
    evals = cmat::Zero(D, 1);
    evals << 0.5, 0.5;
    A = evals.asDiagonal();
    U = randU(D);
    EXPECT_NEAR(1, qpp::entropy(A), 1e-7);

    // 3 x 3 random matrix with fixed equal eigenvalues
    D = 3;
    evals = cmat::Zero(D, 1);
    evals << 1 / 3., 1 / 3., 1 / 3.;
    A = evals.asDiagonal();
    U = randU(D);
    EXPECT_NEAR(std::log2(3), qpp::entropy(A), 1e-7);
}
/******************************************************************************/
/// BEGIN inline double qpp::entropy(const std::vector<double>& prob)
TEST(qpp_entropy_vector, AllTests) {
    // 1 value
    std::vector<double> v = {1};
    EXPECT_NEAR(0, qpp::entropy(v), 1e-7);

    // 2 values, only 1 non-zero
    v = {1, 0};
    EXPECT_NEAR(0, qpp::entropy(v), 1e-7);

    // 2 fixed values
    v = {0.6, 0.4};
    EXPECT_NEAR(0.970950594455, qpp::entropy(v), 1e-7);

    // 2 equal values
    v = {0.5, 0.5};
    EXPECT_NEAR(1, qpp::entropy(v), 1e-7);

    // 3 equal values
    v = {1 / 3., 1 / 3., 1 / 3.};
    idx D = v.size();
    EXPECT_NEAR(std::log2(D), qpp::entropy(v), 1e-7);
}
/******************************************************************************/
/// BEGIN template<typename Derived> double qpp::qmutualinfo(
///       const Eigen::MatrixBase<Derived>& A,
///       const std::vector<idx>& subsysA,
///       const std::vector<idx>& subsysB,
///       const std::vector<idx>& dims)
TEST(qpp_qmutualinfo, AllTests) {
    // 1 x 1 degenerate product state
    idx dA = 1, dB = 1;
    cmat rhoA = randrho(dA), rhoB = randrho(dB);
    cmat rho = kron(rhoA, rhoB);
    double result = qpp::qmutualinfo(rho, {0}, {1}, {dA, dB});
    double expected = 0;
    EXPECT_NEAR(result, expected, 1e-7);

    // 2 x 3 product state
    dA = 2, dB = 3;
    rhoA = randrho(dA), rhoB = randrho(dB);
    rho = kron(rhoA, rhoB);
    result = qpp::qmutualinfo(rho, {0}, {1}, {dA, dB});
    expected = 0;
    EXPECT_NEAR(result, expected, 1e-7);

    // 3 x 3 maximally entangled state
    dA = dB = 3;
    rho = prj(st.mes(dA));
    result = qpp::qmutualinfo(rho, {0}, {1}, {dA, dB});
    expected = 2 * std::log2(dA);
    EXPECT_NEAR(result, expected, 1e-7);

    // 3 x 3 maximally entangled state tensored with a 2 x 2 product state
    rho = prj(st.mes(3)); // MES
    rho = kron(rho, kron(randrho(2), randrho(2)));
    rho = syspermute(rho, {0, 2, 1, 3}, {3, 3, 2, 2});
    result = qpp::qmutualinfo(rho, {0}, {2}, {3, 2, 3, 2});
    expected = 2 * std::log2(3);
    EXPECT_NEAR(result, expected, 1e-7);

    // random 3 x 4 state
    dA = 3, dB = 4;
    rho = randrho(dA * dB);
    rhoA = ptrace2(rho, {dA, dB});
    rhoB = ptrace1(rho, {dA, dB});
    result = qpp::qmutualinfo(rho, {0}, {1}, {dA, dB});
    expected = entropy(rhoA) + entropy(rhoB) - entropy(rho);
    EXPECT_NEAR(result, expected, 1e-7);
}
/******************************************************************************/
/// BEGIN template<typename Derived> double qpp::qmutualinfo(
///       const Eigen::MatrixBase<Derived>& A,
///       const std::vector<idx>& subsysA,
///       const std::vector<idx>& subsysB,
///       idx d = 2)
TEST(qpp_qmutualinfo_qubits, AllTests) {
    // 2 x 2 product state
    idx d = 2;
    cmat rhoA = randrho(d), rhoB = randrho(d);
    cmat rho = kron(rhoA, rhoB);
    double result = qpp::qmutualinfo(rho, {0}, {1});
    double expected = 0;
    EXPECT_NEAR(result, expected, 1e-7);

    // 3 x 3 maximally entangled state
    d = 3;
    rho = prj(st.mes(d));
    result = qpp::qmutualinfo(rho, {0}, {1}, d);
    expected = 2 * std::log2(d);
    EXPECT_NEAR(result, expected, 1e-7);

    // 3 x 3 maximally entangled state tensored with a 3 x 3 product state
    rho = prj(st.mes(3)); // MES
    rho = kron(rho, kron(randrho(3), randrho(3)));
    rho = syspermute(rho, {0, 2, 1, 3}, {3, 3, 3, 3});
    result = qpp::qmutualinfo(rho, {0}, {2}, 3);
    expected = 2 * std::log2(3);
    EXPECT_NEAR(result, expected, 1e-7);

    // random 3 x 3 state
    d = 3;
    rho = randrho(d * d);
    rhoA = ptrace2(rho, d);
    rhoB = ptrace1(rho, d);
    result = qpp::qmutualinfo(rho, {0}, {1}, d);
    expected = entropy(rhoA) + entropy(rhoB) - entropy(rho);
    EXPECT_NEAR(result, expected, 1e-7);
}
/******************************************************************************/
/// BEGIN template<typename Derived> double qpp::renyi(
///       const Eigen::MatrixBase<Derived>& A, double alpha)
TEST(qpp_renyi_matrix, AllTests) {
    // 1 x 1 case
    cmat A(1, 1);
    A << 1.;
    EXPECT_NEAR(0, qpp::renyi(A, 0), 1e-7);
    EXPECT_NEAR(0, qpp::renyi(A, 1 / 2.), 1e-7);
    EXPECT_NEAR(0, qpp::renyi(A, 1), 1e-7);
    EXPECT_NEAR(0, qpp::renyi(A, 2), 1e-7);
    EXPECT_NEAR(0, qpp::renyi(A, qpp::infty), 1e-7);

    // 2 x 2 random matrix with fixed eigenvalues
    idx D = 2;
    cmat evals = cmat::Zero(D, 1);
    evals << 0.6, 0.4;
    A = evals.asDiagonal();
    cmat U = randU(D);
    EXPECT_NEAR(1, qpp::renyi(A, 0), 1e-7);
    EXPECT_NEAR(0.985351706365, qpp::renyi(A, 1 / 2.), 1e-7);
    EXPECT_NEAR(0.970950594455, qpp::renyi(A, 1), 1e-7);
    EXPECT_NEAR(0.943416471634, qpp::renyi(A, 2), 1e-7);
    EXPECT_NEAR(0.918250633859, qpp::renyi(A, 3), 1e-7);
    EXPECT_NEAR(-std::log2(0.6), qpp::renyi(A, qpp::infty), 1e-7);

    // 2 x 2 random matrix with fixed equal eigenvalues
    D = 2;
    evals = cmat::Zero(D, 1);
    evals << 0.5, 0.5;
    A = evals.asDiagonal();
    U = randU(D);
    EXPECT_NEAR(1, qpp::renyi(A, 0), 1e-7);
    EXPECT_NEAR(1, qpp::renyi(A, 1 / 2.), 1e-7);
    EXPECT_NEAR(1, qpp::renyi(A, 1), 1e-7);
    EXPECT_NEAR(1, qpp::renyi(A, 2), 1e-7);
    EXPECT_NEAR(1, qpp::renyi(A, 3), 1e-7);
    EXPECT_NEAR(1, qpp::renyi(A, qpp::infty), 1e-7);

    // 3 x 3 random matrix with fixed only 1 non-zero eigenvalue
    D = 3;
    evals = cmat::Zero(D, 1);
    evals << 1., 0., 0.;
    A = evals.asDiagonal();
    U = randU(D);
    EXPECT_NEAR(std::log2(D), qpp::renyi(A, 0), 1e-7);
    EXPECT_NEAR(0, qpp::renyi(A, 1 / 2.), 1e-7);
    EXPECT_NEAR(0, qpp::renyi(A, 1), 1e-7);
    EXPECT_NEAR(0, qpp::renyi(A, 2), 1e-7);
    EXPECT_NEAR(0, qpp::renyi(A, qpp::infty), 1e-7);

    // 3 x 3 random matrix with fixed equal eigenvalues
    D = 3;
    evals = cmat::Zero(D, 1);
    evals << 1 / 3., 1 / 3., 1 / 3.;
    A = evals.asDiagonal();
    U = randU(D);
    EXPECT_NEAR(std::log2(D), qpp::renyi(A, 0), 1e-7);
    EXPECT_NEAR(std::log2(D), qpp::renyi(A, 1 / 2.), 1e-7);
    EXPECT_NEAR(std::log2(D), qpp::renyi(A, 1), 1e-7);
    EXPECT_NEAR(std::log2(D), qpp::renyi(A, 2), 1e-7);
    EXPECT_NEAR(std::log2(D), qpp::renyi(A, 3), 1e-7);
    EXPECT_NEAR(std::log2(D), qpp::renyi(A, qpp::infty), 1e-7);
}
/******************************************************************************/
/// BEGIN qpp::inline double renyi(const std::vector<double>& prob,
///       double alpha)
TEST(qpp_renyi_vector, AllTests) {
    // 1 value
    std::vector<double> v = {1};
    EXPECT_NEAR(0, qpp::renyi(v, 0), 1e-7);
    EXPECT_NEAR(0, qpp::renyi(v, 1 / 2.), 1e-7);
    EXPECT_NEAR(0, qpp::renyi(v, 1), 1e-7);
    EXPECT_NEAR(0, qpp::renyi(v, 2), 1e-7);
    EXPECT_NEAR(0, qpp::renyi(v, qpp::infty), 1e-7);

    // 2 fixed values
    v = {0.6, 0.4};
    EXPECT_NEAR(1, qpp::renyi(v, 0), 1e-7);
    EXPECT_NEAR(0.985351706365, qpp::renyi(v, 1 / 2.), 1e-7);
    EXPECT_NEAR(0.970950594455, qpp::renyi(v, 1), 1e-7);
    EXPECT_NEAR(0.943416471634, qpp::renyi(v, 2), 1e-7);
    EXPECT_NEAR(0.918250633859, qpp::renyi(v, 3), 1e-7);
    EXPECT_NEAR(-std::log2(0.6), qpp::renyi(v, qpp::infty), 1e-7);

    // 2 equal values
    v = {0.5, 0.5};
    EXPECT_NEAR(1, qpp::renyi(v, 0), 1e-7);
    EXPECT_NEAR(1, qpp::renyi(v, 1 / 2.), 1e-7);
    EXPECT_NEAR(1, qpp::renyi(v, 1), 1e-7);
    EXPECT_NEAR(1, qpp::renyi(v, 2), 1e-7);
    EXPECT_NEAR(1, qpp::renyi(v, 3), 1e-7);
    EXPECT_NEAR(1, qpp::renyi(v, qpp::infty), 1e-7);

    // 3 values, only 1 non-zero
    v = {1, 0, 0};
    idx D = 3;
    EXPECT_NEAR(std::log2(D), qpp::renyi(v, 0), 1e-7);
    EXPECT_NEAR(0, qpp::renyi(v, 1 / 2.), 1e-7);
    EXPECT_NEAR(0, qpp::renyi(v, 1), 1e-7);
    EXPECT_NEAR(0, qpp::renyi(v, 2), 1e-7);
    EXPECT_NEAR(0, qpp::renyi(v, qpp::infty), 1e-7);

    // 3 equal values
    v = {1 / 3., 1 / 3., 1 / 3.};
    EXPECT_NEAR(std::log2(D), qpp::renyi(v, 0), 1e-7);
    EXPECT_NEAR(std::log2(D), qpp::renyi(v, 1 / 2.), 1e-7);
    EXPECT_NEAR(std::log2(D), qpp::renyi(v, 1), 1e-7);
    EXPECT_NEAR(std::log2(D), qpp::renyi(v, 2), 1e-7);
    EXPECT_NEAR(std::log2(D), qpp::renyi(v, 3), 1e-7);
    EXPECT_NEAR(std::log2(D), qpp::renyi(v, qpp::infty), 1e-7);
}
/******************************************************************************/
/// BEGIN template<typename Derived> double qpp::tsallis(
///       const Eigen::MatrixBase<Derived>& A, double q)
TEST(qpp_tsallis_matrix, AllTests) {
    // 1 x 1 case
    cmat A(1, 1);
    A << 1.;
    EXPECT_NEAR(0, qpp::tsallis(A, 0), 1e-7);
    EXPECT_NEAR(0, qpp::tsallis(A, 1 / 2.), 1e-7);
    EXPECT_NEAR(0, qpp::tsallis(A, 1), 1e-7);
    EXPECT_NEAR(0, qpp::tsallis(A, 2), 1e-7);
    EXPECT_NEAR(0, qpp::tsallis(A, qpp::infty), 1e-7);

    // 2 x 2 random matrix with fixed eigenvalues
    idx D = 2;
    cmat evals = cmat::Zero(D, 1);
    evals << 0.6, 0.4;
    A = evals.asDiagonal();
    cmat U = randU(D);
    EXPECT_NEAR(1, qpp::tsallis(A, 0), 1e-7);
    EXPECT_NEAR(0.81410440255, qpp::tsallis(A, 1 / 2.), 1e-7);
    EXPECT_NEAR(qpp::entropy(A) * std::log(2), qpp::tsallis(A, 1), 1e-7);
    EXPECT_NEAR(0.48, qpp::tsallis(A, 2), 1e-7);
    EXPECT_NEAR(0.36, qpp::tsallis(A, 3), 1e-7);
    EXPECT_NEAR(0, qpp::tsallis(A, qpp::infty), 1e-7);

    // 2 x 2 random matrix with fixed equal eigenvalues
    D = 2;
    evals = cmat::Zero(D, 1);
    evals << 0.5, 0.5;
    A = evals.asDiagonal();
    U = randU(D);
    EXPECT_NEAR(1, qpp::tsallis(A, 0), 1e-7);
    EXPECT_NEAR(0.828427124746, qpp::tsallis(A, 1 / 2.), 1e-7);
    EXPECT_NEAR(qpp::entropy(A) * std::log(2), qpp::tsallis(A, 1), 1e-7);
    EXPECT_NEAR(0.5, qpp::tsallis(A, 2), 1e-7);
    EXPECT_NEAR(0.375, qpp::tsallis(A, 3), 1e-7);
    EXPECT_NEAR(0, qpp::tsallis(A, qpp::infty), 1e-7);

    // 3 x 3 random matrix with fixed only 1 non-zero eigenvalue
    D = 3;
    evals = cmat::Zero(D, 1);
    evals << 1., 0., 0.;
    A = evals.asDiagonal();
    U = randU(D);
    EXPECT_NEAR(0, qpp::tsallis(A, 0), 1e-7);
    EXPECT_NEAR(0, qpp::tsallis(A, 1 / 2.), 1e-7);
    EXPECT_NEAR(0, qpp::tsallis(A, 1), 1e-7);
    EXPECT_NEAR(0, qpp::tsallis(A, 2), 1e-7);
    EXPECT_NEAR(0, qpp::tsallis(A, qpp::infty), 1e-7);

    // 3 x 3 random matrix with fixed equal eigenvalues
    D = 3;
    evals = cmat::Zero(D, 1);
    evals << 1 / 3., 1 / 3., 1 / 3.;
    A = evals.asDiagonal();
    U = randU(D);
    EXPECT_NEAR(2, qpp::tsallis(A, 0), 1e-7);
    EXPECT_NEAR(1.46410161514, qpp::tsallis(A, 1 / 2.), 1e-7);
    EXPECT_NEAR(qpp::entropy(A) * std::log(2), qpp::tsallis(A, 1), 1e-7);
    EXPECT_NEAR(2 / 3., qpp::tsallis(A, 2), 1e-7);
    EXPECT_NEAR(4 / 9., qpp::tsallis(A, 3), 1e-7);
    EXPECT_NEAR(0, qpp::tsallis(A, qpp::infty), 1e-7);
}
/******************************************************************************/
/// BEGIN inline double qpp::tsallis(const std::vector<double>& prob, double q)
TEST(qpp_tsallis_vector, AllTests) {
    // 1 value
    std::vector<double> v = {1};
    EXPECT_NEAR(0, qpp::tsallis(v, 0), 1e-7);
    EXPECT_NEAR(0, qpp::tsallis(v, 1 / 2.), 1e-7);
    EXPECT_NEAR(0, qpp::tsallis(v, 1), 1e-7);
    EXPECT_NEAR(0, qpp::tsallis(v, 2), 1e-7);
    EXPECT_NEAR(0, qpp::tsallis(v, qpp::infty), 1e-7);

    // 2 fixed values
    v = {0.6, 0.4};
    EXPECT_NEAR(1, qpp::tsallis(v, 0), 1e-7);
    EXPECT_NEAR(0.81410440255, qpp::tsallis(v, 1 / 2.), 1e-7);
    EXPECT_NEAR(qpp::entropy(v) * std::log(2), qpp::tsallis(v, 1), 1e-7);
    EXPECT_NEAR(0.48, qpp::tsallis(v, 2), 1e-7);
    EXPECT_NEAR(0.36, qpp::tsallis(v, 3), 1e-7);
    EXPECT_NEAR(0, qpp::tsallis(v, qpp::infty), 1e-7);

    // 2 equal values
    v = {0.5, 0.5};
    EXPECT_NEAR(1, qpp::tsallis(v, 0), 1e-7);
    EXPECT_NEAR(0.828427124746, qpp::tsallis(v, 1 / 2.), 1e-7);
    EXPECT_NEAR(qpp::entropy(v) * std::log(2), qpp::tsallis(v, 1), 1e-7);
    EXPECT_NEAR(0.5, qpp::tsallis(v, 2), 1e-7);
    EXPECT_NEAR(0.375, qpp::tsallis(v, 3), 1e-7);
    EXPECT_NEAR(0, qpp::tsallis(v, qpp::infty), 1e-7);

    // 3 values, only 1 non-zero
    v = {1, 0, 0};
    EXPECT_NEAR(0, qpp::tsallis(v, 0), 1e-7);
    EXPECT_NEAR(0, qpp::tsallis(v, 1 / 2.), 1e-7);
    EXPECT_NEAR(0, qpp::tsallis(v, 1), 1e-7);
    EXPECT_NEAR(0, qpp::tsallis(v, 2), 1e-7);
    EXPECT_NEAR(0, qpp::tsallis(v, qpp::infty), 1e-7);

    // 3 equal values
    v = {1 / 3., 1 / 3., 1 / 3.};
    EXPECT_NEAR(2, qpp::tsallis(v, 0), 1e-7);
    EXPECT_NEAR(1.46410161514, qpp::tsallis(v, 1 / 2.), 1e-7);
    EXPECT_NEAR(qpp::entropy(v) * std::log(2), qpp::tsallis(v, 1), 1e-7);
    EXPECT_NEAR(2 / 3., qpp::tsallis(v, 2), 1e-7);
    EXPECT_NEAR(4 / 9., qpp::tsallis(v, 3), 1e-7);
    EXPECT_NEAR(0, qpp::tsallis(v, qpp::infty), 1e-7);
}
/******************************************************************************/
