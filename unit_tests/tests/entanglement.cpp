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
#include <cmath>
#include "gtest/gtest.h"
#include "qpp.h"

using namespace qpp;

// Unit testing "entanglement.h"

/******************************************************************************/
/// BEGIN template<typename Derived> double qpp::concurrence(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_concurrence, AllTests) {
    // random qubit product state
    cmat rho1 = prj(kron(randU(), randU()) * mket({0, 0}));
    EXPECT_NEAR(0, qpp::concurrence(rho1), 1e-7);

    // random maximally entangled 2-qubit state
    cmat rho2 = prj(kron(randU(), randU()) * st.b00);
    EXPECT_NEAR(1, qpp::concurrence(rho2), 1e-7);

    // random 2-qubit state with Schmidt coefficients 0.8 and 0.6
    cmat rho3 =
        prj(kron(randU(), randU()) * (0.8 * mket({0, 0}) + 0.6 * mket({1, 1})));
    EXPECT_NEAR(2 * 0.8 * 0.6, qpp::concurrence(rho3), 1e-7);
}
/******************************************************************************/
/// BEGIN template<typename Derived> double qpp::entanglement(
///       const Eigen::MatrixBase<Derived>& A, const std::vector<idx>& dims)
TEST(qpp_entanglement, AllTests) {
    // random qutrit product state
    ket psi1 = kron(randU(3), randU(3)) * mket({0, 0}, 3);
    EXPECT_NEAR(0, qpp::entanglement(psi1, {3, 3}), 1e-7);

    // random 2-qutrit state with Schmidt coefficients 0.36, 0.09 and 0.01
    ket psi2 =
        kron(randU(3), randU(3)) *
        (0.6 * mket({0, 0}, 3) + 0.3 * mket({1, 1}, 3) + 0.1 * mket({2, 2}, 3));
    EXPECT_NEAR(-0.36 * std::log2(0.36) - 0.09 * std::log2(0.09) -
                    0.01 * std::log2(0.01),
                qpp::entanglement(psi2, {3, 3}), 1e-7);

    // random maximally entangled 2-qutrit state
    ket psi3 = kron(randU(3), randU(3)) * st.mes(3);
    EXPECT_NEAR(std::log2(3), qpp::entanglement(psi3, {3, 3}), 1e-7);
}
/******************************************************************************/
/// BEGIN template<typename Derived> double qpp::entanglement(
///       const Eigen::MatrixBase<Derived>& A, idx d = 2)
TEST(qpp_entanglement_qubits, AllTests) {
    // random qubit product state
    ket psi1 = kron(randU(), randU()) * mket({0, 0});
    EXPECT_NEAR(0, qpp::entanglement(psi1), 1e-7);

    // random maximally entangled 2-qubit state
    ket psi2 = kron(randU(), randU()) * st.b00;
    EXPECT_NEAR(1, qpp::entanglement(psi2), 1e-7);

    // random 2-qubit state with Schmidt coefficients 0.8 and 0.6
    ket psi3 =
        kron(randU(), randU()) * (0.8 * mket({0, 0}) + 0.6 * mket({1, 1}));
    EXPECT_NEAR(-0.64 * std::log2(0.64) - 0.36 * std::log2(0.36),
                qpp::entanglement(psi3), 1e-7);
}
/******************************************************************************/
/// BEGIN template<typename Derived> double qpp::gconcurrence(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_gconcurrence, AllTests) {
    // random qubit product state
    ket psi1 = kron(randU(), randU()) * mket({0, 0});
    EXPECT_NEAR(0, qpp::gconcurrence(psi1), 1e-7);

    // random maximally entangled 2-qubit state
    ket psi2 = kron(randU(), randU()) * st.b00;
    EXPECT_NEAR(1, qpp::gconcurrence(psi2), 1e-7);

    // random 2-qubit state with Schmidt coefficients 0.8 and 0.6
    ket psi3 =
        kron(randU(), randU()) * (0.8 * mket({0, 0}) + 0.6 * mket({1, 1}));
    EXPECT_NEAR(2 * 0.8 * 0.6, qpp::gconcurrence(psi3), 1e-7);

    // random maximally entangled 2-qutrit state
    ket psi4 = kron(randU(3), randU(3)) * st.mes(3);
    EXPECT_NEAR(1, qpp::gconcurrence(psi4), 1e-7);

    // random 2-qutrit state with Schmidt coefficients 0.36, 0.09 and 0.01
    ket psi5 =
        kron(randU(3), randU(3)) *
        (0.6 * mket({0, 0}, 3) + 0.3 * mket({1, 1}, 3) + 0.1 * mket({2, 2}, 3));
    EXPECT_NEAR(3 * std::pow(0.36 * 0.09 * 0.01, 1. / 3),
                qpp::gconcurrence(psi5), 1e-7);

    // for qubits, expect the gconcurrence to be the same as concurrence
    ket psi6 = randket(4);
    EXPECT_NEAR(gconcurrence(psi6), concurrence(prj(psi6)), 1e-7);
}
/******************************************************************************/
/// BEGIN template<typename Derived> double qpp::lognegativity(
///       const Eigen::MatrixBase<Derived>& A, const std::vector<idx>& dims)
TEST(qpp_lognegativity, AllTests) {
    // zero on product states (2 qutrits)
    cmat rho = randrho(3);
    cmat sigma = randrho(3);
    EXPECT_NEAR(0, qpp::lognegativity(kron(rho, sigma), {3, 3}), 1e-7);

    // additivity (2 ququads)
    rho = randrho(4);
    sigma = randrho(4);
    EXPECT_NEAR(
        qpp::lognegativity(syspermute(kron(rho, sigma), {0, 2, 1, 3}), {4, 4}),
        qpp::lognegativity(rho, {2, 2}) + qpp::lognegativity(sigma, {2, 2}),
        1e-7);

    // random maximally entangled 2-qutrit state
    idx d = 3;
    rho = prj(kron(randU(d), randU(d)) * st.mes(d));
    EXPECT_NEAR(std::log2(d), qpp::lognegativity(rho, {d, d}), 1e-7);

    // random maximally entangled 2-ququad state
    d = 4;
    rho = prj(kron(randU(d), randU(d)) * st.mes(d));
    EXPECT_NEAR(std::log2(d), qpp::lognegativity(rho, {d, d}), 1e-7);

    // random maximally entangled state (d = 7)
    d = 7;
    rho = prj(kron(randU(d), randU(d)) * st.mes(d));
    EXPECT_NEAR(std::log2(d), qpp::lognegativity(rho, {d, d}), 1e-7);
}
/******************************************************************************/
/// BEGIN template<typename Derived> double qpp::lognegativity(
///       const Eigen::MatrixBase<Derived>& A, idx d = 2)
TEST(qpp_lognegativity_qubits, AllTests) {
    // zero on product states (2 qubits)
    cmat rho = randrho();
    cmat sigma = randrho();
    EXPECT_NEAR(0, qpp::lognegativity(kron(rho, sigma)), 1e-7);

    // random 2-qubit state with Schmidt coefficients 0.8 and 0.6
    ket psi =
        kron(randU(), randU()) * (0.8 * mket({0, 0}) + 0.6 * mket({1, 1}));
    rho = prj(psi);
    EXPECT_NEAR(0.9708536, qpp::lognegativity(rho), 1e-7);

    // random maximally entangled 2-qubit state
    idx d = 2;
    rho = prj(kron(randU(d), randU(d)) * st.mes(d));
    EXPECT_NEAR(std::log2(d), qpp::lognegativity(rho), 1e-7);
}
/******************************************************************************/
/// BEGIN template<typename Derived> double qpp::negativity(
///       const Eigen::MatrixBase<Derived>& A, const std::vector<idx>& dims)
TEST(qpp_negativity, AllTests) {
    // Must be (d - 1)/2 on MES

    // zero on product states (2 qutrits)
    cmat rho = randrho(3);
    cmat sigma = randrho(3);
    EXPECT_NEAR(0, qpp::negativity(kron(rho, sigma), {3, 3}), 1e-7);

    // convexity (10 ququads)
    idx N = 10;
    idx d = 4;
    std::vector<cmat> rhos(N);
    std::vector<double> probs = randprob(N);
    rho = cmat::Zero(d * d, d * d);
    double sum_neg = 0;
    for (idx i = 0; i < N; ++i) {
        cmat rho_i = randrho(d * d);
        sum_neg += probs[i] * qpp::negativity(rho_i, {d, d});
        rho += probs[i] * rho_i;
    }
    EXPECT_LE(qpp::negativity(rho, {d, d}), sum_neg);

    // random maximally entangled 2-qutrit state
    d = 3;
    rho = prj(kron(randU(d), randU(d)) * st.mes(d));
    EXPECT_NEAR((d - 1) / 2., qpp::negativity(rho, {d, d}), 1e-7);

    // random maximally entangled 2-ququad state
    d = 4;
    rho = prj(kron(randU(d), randU(d)) * st.mes(d));
    EXPECT_NEAR((d - 1) / 2., qpp::negativity(rho, {d, d}), 1e-7);

    // random maximally entangled state (d = 7)
    d = 7;
    rho = prj(kron(randU(d), randU(d)) * st.mes(d));
    EXPECT_NEAR((d - 1) / 2., qpp::negativity(rho, {d, d}), 1e-7);
}
/******************************************************************************/
/// BEGIN template<typename Derived> double qpp::negativity(
///       const Eigen::MatrixBase<Derived>& A, idx d = 2)
TEST(qpp_negativity_qubits, AllTests) {
    // Must be (d - 1)/2 on MES

    // zero on product states (2 qutrits)
    cmat rho = randrho();
    cmat sigma = randrho();
    EXPECT_NEAR(0, qpp::negativity(kron(rho, sigma)), 1e-7);

    // random 2-qubit state with Schmidt coefficients 0.8 and 0.6
    ket psi =
        kron(randU(), randU()) * (0.8 * mket({0, 0}) + 0.6 * mket({1, 1}));
    rho = prj(psi);
    EXPECT_NEAR(0.48, qpp::negativity(rho), 1e-7);

    // random maximally entangled 2-qubit state
    psi = kron(randU(), randU()) * st.mes();
    rho = prj(psi);
    EXPECT_NEAR(0.5, qpp::negativity(rho), 1e-7);
}
/******************************************************************************/
/// BEGIN template<typename Derived> cmat qpp::schmidtA/B(
///       const Eigen::MatrixBase<Derived>& A, const std::vector<idx>& dims)
TEST(qpp_schmidtA_schmidtB, AllTests) {
    // random degenerate 1 x 1 product state
    idx dA = 1, dB = 1, D = dA * dB, minD = std::min(dA, dB);
    cmat UA = randU(dA);
    cmat UB = randU(dB);
    ket psi = kron(UA, UB) * mket({0, 0}, {dA, dB});
    cmat basisA = qpp::schmidtA(psi, {dA, dB});
    cmat basisB = qpp::schmidtB(psi, {dA, dB});
    // unitarity
    EXPECT_NEAR(0, norm(adjoint(basisA) * basisA - gt.Id(dA)), 1e-7);
    EXPECT_NEAR(0, norm(adjoint(basisB) * basisB - gt.Id(dB)), 1e-7);
    // get the Schmidt coefficients and test the result
    dyn_col_vect<double> scf = schmidtcoeffs(psi, {dA, dB});
    ket expected = ket::Zero(D);
    for (idx i = 0; i < minD; ++i) {
        expected += scf(i) * kron(basisA.col(i), basisB.col(i));
    }
    EXPECT_NEAR(0, norm(expected - psi), 1e-7);

    // random degenerate 3 x 1 product state
    dA = 3, dB = 1, D = dA * dB, minD = std::min(dA, dB);
    UA = randU(dA);
    UB = randU(dB);
    psi = kron(UA, UB) * mket({0, 0}, {dA, dB});
    basisA = qpp::schmidtA(psi, {dA, dB});
    basisB = qpp::schmidtB(psi, {dA, dB});
    // unitarity
    EXPECT_NEAR(0, norm(adjoint(basisA) * basisA - gt.Id(dA)), 1e-7);
    EXPECT_NEAR(0, norm(adjoint(basisB) * basisB - gt.Id(dB)), 1e-7);
    // get the Schmidt coefficients and test the result
    scf = schmidtcoeffs(psi, {dA, dB});
    expected = ket::Zero(D);
    for (idx i = 0; i < minD; ++i) {
        expected += scf(i) * kron(basisA.col(i), basisB.col(i));
    }
    EXPECT_NEAR(0, norm(expected - psi), 1e-7);

    // random 2 x 4 product state
    dA = 2, dB = 4, D = dA * dB, minD = std::min(dA, dB);
    UA = randU(dA);
    UB = randU(dB);
    psi = kron(UA, UB) * mket({0, 0}, {dA, dB});
    basisA = qpp::schmidtA(psi, {dA, dB});
    basisB = qpp::schmidtB(psi, {dA, dB});
    // unitarity
    EXPECT_NEAR(0, norm(adjoint(basisA) * basisA - gt.Id(dA)), 1e-7);
    EXPECT_NEAR(0, norm(adjoint(basisB) * basisB - gt.Id(dB)), 1e-7);
    // get the Schmidt coefficients and test the result
    scf = schmidtcoeffs(psi, {dA, dB});
    expected = ket::Zero(D);
    for (idx i = 0; i < minD; ++i) {
        expected += scf(i) * kron(basisA.col(i), basisB.col(i));
    }
    EXPECT_NEAR(0, norm(expected - psi), 1e-7);

    // random 3 x 4 state
    dA = 3, dB = 4, D = dA * dB, minD = std::min(dA, dB);
    psi = randket(D);
    basisA = qpp::schmidtA(psi, {dA, dB});
    basisB = qpp::schmidtB(psi, {dA, dB});
    // unitarity
    EXPECT_NEAR(0, norm(adjoint(basisA) * basisA - gt.Id(dA)), 1e-7);
    EXPECT_NEAR(0, norm(adjoint(basisB) * basisB - gt.Id(dB)), 1e-7);
    // get the Schmidt coefficients and test the result
    scf = schmidtcoeffs(psi, {dA, dB});
    expected = ket::Zero(D);
    for (idx i = 0; i < minD; ++i) {
        expected += scf(i) * kron(basisA.col(i), basisB.col(i));
    }
    EXPECT_NEAR(0, norm(expected - psi), 1e-7);
}
/******************************************************************************/
/// BEGIN template<typename Derived> cmat qpp::schmidtA/B(
///       const Eigen::MatrixBase<Derived>& A, idx d = 2)
TEST(qpp_schmidtA_schmidtB_qubits, AllTests) {
    // random 2 x 2 product state
    idx d = 2, D = d * d;
    cmat UA = randU(d);
    cmat UB = randU(d);
    ket psi = kron(UA, UB) * st.zero(2, d);
    cmat basisA = qpp::schmidtA(psi, d);
    cmat basisB = qpp::schmidtB(psi, d);
    // unitarity
    EXPECT_NEAR(0, norm(adjoint(basisA) * basisA - gt.Id(d)), 1e-7);
    EXPECT_NEAR(0, norm(adjoint(basisB) * basisB - gt.Id(d)), 1e-7);
    // get the Schmidt coefficients and test the result
    dyn_col_vect<double> scf = schmidtcoeffs(psi, d);
    ket expected = ket::Zero(D);
    for (idx i = 0; i < d; ++i) {
        expected += scf(i) * kron(basisA.col(i), basisB.col(i));
    }
    EXPECT_NEAR(0, norm(expected - psi), 1e-7);

    // random 4 x 4 product state
    d = 4, D = d * d;
    UA = randU(d);
    UB = randU(d);
    psi = kron(UA, UB) * st.zero(2, d);
    basisA = qpp::schmidtA(psi, d);
    basisB = qpp::schmidtB(psi, d);
    // unitarity
    EXPECT_NEAR(0, norm(adjoint(basisA) * basisA - gt.Id(d)), 1e-7);
    EXPECT_NEAR(0, norm(adjoint(basisB) * basisB - gt.Id(d)), 1e-7);
    // get the Schmidt coefficients and test the result
    scf = schmidtcoeffs(psi, d);
    expected = ket::Zero(D);
    for (idx i = 0; i < d; ++i) {
        expected += scf(i) * kron(basisA.col(i), basisB.col(i));
    }
    EXPECT_NEAR(0, norm(expected - psi), 1e-7);

    // random 5 x 5 state
    d = 5, D = d * d;
    psi = randket(D);
    basisA = qpp::schmidtA(psi, d);
    basisB = qpp::schmidtB(psi, d);
    // unitarity
    EXPECT_NEAR(0, norm(adjoint(basisA) * basisA - gt.Id(d)), 1e-7);
    EXPECT_NEAR(0, norm(adjoint(basisB) * basisB - gt.Id(d)), 1e-7);
    // get the Schmidt coefficients and test the result
    scf = schmidtcoeffs(psi, d);
    expected = ket::Zero(D);
    for (idx i = 0; i < d; ++i) {
        expected += scf(i) * kron(basisA.col(i), basisB.col(i));
    }
    EXPECT_NEAR(0, norm(expected - psi), 1e-7);
}
/******************************************************************************/
/// BEGIN template<typename Derived> dyn_col_vect<double> qpp::schmidtcoeffs(
///       const Eigen::MatrixBase<Derived>& A, const std::vector<idx>& dims)
TEST(qpp_schmidtcoeffs, AllTests) {
    // random degenerate 1 x 1 product state
    idx dA = 1, dB = 1, D = dA * dB, minD = std::min(dA, dB);
    cmat UA = randU(dA);
    cmat UB = randU(dB);
    ket psi = kron(UA, UB) * mket({0, 0}, {dA, dB});
    dyn_col_vect<double> result = qpp::schmidtcoeffs(psi, {dA, dB});
    dyn_col_vect<double> expected(minD);
    expected << 1;
    EXPECT_NEAR(0, norm(result - expected), 1e-7);

    // random degenerate 3 x 1 product state
    dA = 3, dB = 1, D = dA * dB, minD = std::min(dA, dB);
    UA = randU(dA);
    UB = randU(dB);
    psi = kron(UA, UB) * mket({0, 0}, {dA, dB});
    result = qpp::schmidtcoeffs(psi, {dA, dB});
    expected = dyn_col_vect<double>(minD);
    expected << 1;
    EXPECT_NEAR(0, norm(result - expected), 1e-7);

    // random 3 x 4 product state
    dA = 3, dB = 4, D = dA * dB, minD = std::min(dA, dB);
    UA = randU(dA);
    UB = randU(dB);
    psi = kron(UA, UB) * mket({0, 0}, {dA, dB});
    result = qpp::schmidtcoeffs(psi, {dA, dB});
    expected = dyn_col_vect<double>::Zero(minD);
    expected(0) = 1;
    EXPECT_NEAR(0, norm(result - expected), 1e-7);

    // random 2 x 3 state with fixed Schmidt coefficients
    dA = 2, dB = 3, D = dA * dB, minD = std::min(dA, dB);
    double c0 = 0.8, c1 = 0.6;
    UA = randU(dA);
    UB = randU(dB);
    psi = kron(UA, UB) *
          (c0 * mket({0, 0}, {dA, dB}) + c1 * mket({1, 1}, {dA, dB}));
    result = qpp::schmidtcoeffs(psi, {dA, dB});
    expected = dyn_col_vect<double>::Zero(minD);
    expected << c0, c1;
    EXPECT_NEAR(0, norm(result - expected), 1e-7);

    // random 5 x 3 state with fixed Schmidt coefficients
    dA = 5, dB = 3, D = dA * dB, minD = std::min(dA, dB);
    c0 = 0.8, c1 = 0.5;
    double c2 = std::sqrt(1 - c0 * c0 - c1 * c1);
    UA = randU(dA);
    UB = randU(dB);
    psi = kron(UA, UB) *
          (c0 * mket({0, 0}, {dA, dB}) + c1 * mket({1, 1}, {dA, dB}) +
           c2 * mket({2, 2}, {dA, dB}));
    result = qpp::schmidtcoeffs(psi, {dA, dB});
    expected = dyn_col_vect<double>::Zero(minD);
    expected << c0, c1, c2;
    EXPECT_NEAR(0, norm(result - expected), 1e-7);
}
/******************************************************************************/
/// BEGIN template<typename Derived> dyn_col_vect<double> qpp::schmidtcoeffs(
///       const Eigen::MatrixBase<Derived>& A, idx d = 2)
TEST(qpp_schmidtcoeffs_qubits, AllTests) {
    // random 2 x 2 product state
    idx d = 2;
    cmat UA = randU(d);
    cmat UB = randU(d);
    ket psi = kron(UA, UB) * mket({0, 0});
    dyn_col_vect<double> result = qpp::schmidtcoeffs(psi);
    dyn_col_vect<double> expected = dyn_col_vect<double>::Zero(d);
    expected(0) = 1;
    EXPECT_NEAR(0, norm(result - expected), 1e-7);

    // random 3 x 3 product state
    d = 3;
    UA = randU(d);
    UB = randU(d);
    psi = kron(UA, UB) * mket({1, 1}, d);
    result = qpp::schmidtcoeffs(psi, d);
    expected = dyn_col_vect<double>::Zero(d);
    expected(0) = 1;
    EXPECT_NEAR(0, norm(result - expected), 1e-7);

    // random 2 x 2 state with fixed Schmidt coefficients
    d = 2;
    double c0 = 0.8, c1 = 0.6;
    UA = randU(d);
    UB = randU(d);
    psi = kron(UA, UB) * (c0 * st.zero(2) + c1 * st.one(2));
    result = qpp::schmidtcoeffs(psi);
    expected = dyn_col_vect<double>::Zero(d);
    expected << c0, c1;
    EXPECT_NEAR(0, norm(result - expected), 1e-7);

    // random 3 x 3 state with fixed Schmidt coefficients
    d = 3;
    c0 = 0.8, c1 = 0.5;
    double c2 = std::sqrt(1 - c0 * c0 - c1 * c1);
    UA = randU(d);
    UB = randU(d);
    psi = kron(UA, UB) *
          (c0 * mket({0, 0}, d) + c1 * mket({1, 1}, d) + c2 * mket({2, 2}, d));
    result = qpp::schmidtcoeffs(psi, d);
    expected = dyn_col_vect<double>::Zero(d);
    expected << c0, c1, c2;
    EXPECT_NEAR(0, norm(result - expected), 1e-7);
}
/******************************************************************************/
/// BEGIN template<typename Derived> std::vector<double> qpp::schmidtprobs(
///       const Eigen::MatrixBase<Derived>& A, const std::vector<idx>& dims)
TEST(qpp_schmidtprobs, AllTests) {
    // random degenerate 1 x 1 product state
    idx dA = 1, dB = 1, D = dA * dB, minD = std::min(dA, dB);
    cmat UA = randU(dA);
    cmat UB = randU(dB);
    ket psi = kron(UA, UB) * mket({0, 0}, {dA, dB});
    std::vector<double> result_vect = qpp::schmidtprobs(psi, {dA, dB});
    dyn_col_vect<double> result = Eigen::Map<dyn_col_vect<double>>(
        result_vect.data(), result_vect.size());
    dyn_col_vect<double> expected(minD);
    expected << 1;
    EXPECT_NEAR(0, norm(result - expected), 1e-7);

    // random degenerate 3 x 1 product state
    dA = 3, dB = 1, D = dA * dB, minD = std::min(dA, dB);
    UA = randU(dA);
    UB = randU(dB);
    psi = kron(UA, UB) * mket({0, 0}, {dA, dB});
    result_vect = qpp::schmidtprobs(psi, {dA, dB});
    result = Eigen::Map<dyn_col_vect<double>>(result_vect.data(),
                                              result_vect.size());
    expected = dyn_col_vect<double>(minD);
    expected << 1;
    EXPECT_NEAR(0, norm(result - expected), 1e-7);

    // random 3 x 4 product state
    dA = 3, dB = 4, D = dA * dB, minD = std::min(dA, dB);
    UA = randU(dA);
    UB = randU(dB);
    psi = kron(UA, UB) * mket({0, 0}, {dA, dB});
    result_vect = qpp::schmidtprobs(psi, {dA, dB});
    result = Eigen::Map<dyn_col_vect<double>>(result_vect.data(),
                                              result_vect.size());
    expected = dyn_col_vect<double>::Zero(minD);
    expected(0) = 1;
    EXPECT_NEAR(0, norm(result - expected), 1e-7);

    // random 2 x 3 state with fixed Schmidt coefficients
    dA = 2, dB = 3, D = dA * dB, minD = std::min(dA, dB);
    double c0 = 0.8, c1 = 0.6;
    UA = randU(dA);
    UB = randU(dB);
    psi = kron(UA, UB) *
          (c0 * mket({0, 0}, {dA, dB}) + c1 * mket({1, 1}, {dA, dB}));
    result_vect = qpp::schmidtprobs(psi, {dA, dB});
    result = Eigen::Map<dyn_col_vect<double>>(result_vect.data(),
                                              result_vect.size());
    expected = dyn_col_vect<double>::Zero(minD);
    expected << c0 * c0, c1 * c1;
    EXPECT_NEAR(0, norm(result - expected), 1e-7);

    // random 5 x 3 state with fixed Schmidt coefficients
    dA = 5, dB = 3, D = dA * dB, minD = std::min(dA, dB);
    c0 = 0.8, c1 = 0.5;
    double c2 = std::sqrt(1 - c0 * c0 - c1 * c1);
    UA = randU(dA);
    UB = randU(dB);
    psi = kron(UA, UB) *
          (c0 * mket({0, 0}, {dA, dB}) + c1 * mket({1, 1}, {dA, dB}) +
           c2 * mket({2, 2}, {dA, dB}));
    result_vect = qpp::schmidtprobs(psi, {dA, dB});
    result = Eigen::Map<dyn_col_vect<double>>(result_vect.data(),
                                              result_vect.size());
    expected = dyn_col_vect<double>::Zero(minD);
    expected << c0 * c0, c1 * c1, c2 * c2;
    EXPECT_NEAR(0, norm(result - expected), 1e-7);
}
/******************************************************************************/
/// BEGIN template<typename Derived> std::vector<double> qpp::schmidtprobs(
///       const Eigen::MatrixBase<Derived>& A, idx d = 2)
TEST(qpp_schmidtprobs_qubits, AllTests) {
    // random 2 x 2 product state
    idx d = 2;
    cmat UA = randU(d);
    cmat UB = randU(d);
    ket psi = kron(UA, UB) * mket({0, 0});
    std::vector<double> result_vect = qpp::schmidtprobs(psi);
    dyn_col_vect<double> result = Eigen::Map<dyn_col_vect<double>>(
        result_vect.data(), result_vect.size());
    dyn_col_vect<double> expected = dyn_col_vect<double>::Zero(d);
    expected(0) = 1;
    EXPECT_NEAR(0, norm(result - expected), 1e-7);

    // random 3 x 3 product state
    d = 3;
    UA = randU(d);
    UB = randU(d);
    psi = kron(UA, UB) * mket({1, 1}, d);
    result_vect = qpp::schmidtprobs(psi, d);
    result = Eigen::Map<dyn_col_vect<double>>(result_vect.data(),
                                              result_vect.size());
    expected = dyn_col_vect<double>::Zero(d);
    expected(0) = 1;
    EXPECT_NEAR(0, norm(result - expected), 1e-7);

    // random 2 x 2 state with fixed Schmidt coefficients
    d = 2;
    double c0 = 0.8, c1 = 0.6;
    UA = randU(d);
    UB = randU(d);
    psi = kron(UA, UB) * (c0 * st.zero(2) + c1 * st.one(2));
    result_vect = qpp::schmidtprobs(psi);
    result = Eigen::Map<dyn_col_vect<double>>(result_vect.data(),
                                              result_vect.size());
    expected = dyn_col_vect<double>::Zero(d);
    expected << c0 * c0, c1 * c1;
    EXPECT_NEAR(0, norm(result - expected), 1e-7);

    // random 3 x 3 state with fixed Schmidt coefficients
    d = 3;
    c0 = 0.8, c1 = 0.5;
    double c2 = std::sqrt(1 - c0 * c0 - c1 * c1);
    UA = randU(d);
    UB = randU(d);
    psi = kron(UA, UB) *
          (c0 * mket({0, 0}, d) + c1 * mket({1, 1}, d) + c2 * mket({2, 2}, d));
    result_vect = qpp::schmidtprobs(psi, d);
    result = Eigen::Map<dyn_col_vect<double>>(result_vect.data(),
                                              result_vect.size());
    expected = dyn_col_vect<double>::Zero(d);
    expected << c0 * c0, c1 * c1, c2 * c2;
    EXPECT_NEAR(0, norm(result - expected), 1e-7);
}
/******************************************************************************/
