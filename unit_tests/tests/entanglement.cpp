#include <algorithm>
#include <cmath>

#include "gtest/gtest.h"

#include "qpp.h"

using namespace qpp;

// Unit testing "entanglement.hpp"

/******************************************************************************/
/// BEGIN template <typename Derived> realT concurrence(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_concurrence, AllTests) {
    // random qubit product state
    cmat rho1 = prj(kron(randU(), randU()) * mket({0, 0}));
    EXPECT_NEAR(0, concurrence(rho1), 1e-3);

    // random maximally entangled 2-qubit state
    cmat rho2 = prj(kron(randU(), randU()) * st.b00);
    EXPECT_NEAR(1, concurrence(rho2), 1e-3);

    // random 2-qubit state with Schmidt coefficients 0.8 and 0.6
    cmat rho3 =
        prj(kron(randU(), randU()) * (0.8 * mket({0, 0}) + 0.6 * mket({1, 1})));
    EXPECT_NEAR(2 * 0.8 * 0.6, concurrence(rho3), 1e-3);
}
/******************************************************************************/
/// BEGIN template <typename Derived> realT entanglement(
///       const Eigen::MatrixBase<Derived>& A, const std::vector<idx>& dims)
TEST(qpp_entanglement, Qudits) {
    // random qutrit product state
    ket psi1 = kron(randU(3), randU(3)) * mket({0, 0}, 3);
    EXPECT_NEAR(0, entanglement(psi1, {3, 3}), 1e-5);

    // random 2-qutrit state with Schmidt coefficients 0.36, 0.09 and 0.01
    ket psi2 =
        kron(randU(3), randU(3)) *
        (0.6 * mket({0, 0}, 3) + 0.3 * mket({1, 1}, 3) + 0.1 * mket({2, 2}, 3));
    EXPECT_NEAR(-0.36 * std::log2(0.36) - 0.09 * std::log2(0.09) -
                    0.01 * std::log2(0.01),
                entanglement(psi2, {3, 3}), 1e-5);

    // random maximally entangled 2-qutrit state
    ket psi3 = kron(randU(3), randU(3)) * st.mes(3);
    EXPECT_NEAR(std::log2(3), entanglement(psi3, {3, 3}), 1e-5);
}
/******************************************************************************/
/// BEGIN template <typename Derived> realT entanglement(
///       const Eigen::MatrixBase<Derived>& A, idx d = 2)
TEST(qpp_entanglement, Qubits) {
    // random qubit product state
    ket psi1 = kron(randU(), randU()) * mket({0, 0});
    EXPECT_NEAR(0, entanglement(psi1), 1e-5);

    // random maximally entangled 2-qubit state
    ket psi2 = kron(randU(), randU()) * st.b00;
    EXPECT_NEAR(1, entanglement(psi2), 1e-5);

    // random 2-qubit state with Schmidt coefficients 0.8 and 0.6
    ket psi3 =
        kron(randU(), randU()) * (0.8 * mket({0, 0}) + 0.6 * mket({1, 1}));
    EXPECT_NEAR(-0.64 * std::log2(0.64) - 0.36 * std::log2(0.36),
                entanglement(psi3), 1e-5);
}
/******************************************************************************/
/// BEGIN template <typename Derived> realT gconcurrence(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_gconcurrence, AllTests) {
    // random qubit product state
    ket psi1 = kron(randU(), randU()) * mket({0, 0});
    EXPECT_NEAR(0, gconcurrence(psi1), 1e-3);

    // random maximally entangled 2-qubit state
    ket psi2 = kron(randU(), randU()) * st.b00;
    EXPECT_NEAR(1, gconcurrence(psi2), 1e-3);

    // random 2-qubit state with Schmidt coefficients 0.8 and 0.6
    ket psi3 =
        kron(randU(), randU()) * (0.8 * mket({0, 0}) + 0.6 * mket({1, 1}));
    EXPECT_NEAR(2 * 0.8 * 0.6, gconcurrence(psi3), 1e-3);

    // random maximally entangled 2-qutrit state
    ket psi4 = kron(randU(3), randU(3)) * st.mes(3);
    EXPECT_NEAR(1, gconcurrence(psi4), 1e-3);

    // random 2-qutrit state with Schmidt coefficients 0.36, 0.09 and 0.01
    ket psi5 =
        kron(randU(3), randU(3)) *
        (0.6 * mket({0, 0}, 3) + 0.3 * mket({1, 1}, 3) + 0.1 * mket({2, 2}, 3));
    EXPECT_NEAR(3 * std::pow(0.36 * 0.09 * 0.01, 1. / 3), gconcurrence(psi5),
                1e-3);

    // for qubits, expect the gconcurrence to be the same as concurrence
    ket psi6 = randket(4);
    EXPECT_NEAR(gconcurrence(psi6), concurrence(prj(psi6)), 1e-3);
}
/******************************************************************************/
/// BEGIN template <typename Derived> realT lognegativity(
///       const Eigen::MatrixBase<Derived>& A, const std::vector<idx>& dims)
TEST(qpp_lognegativity, Qudits) {
    // zero on product states (2 qutrits)
    cmat rho = randrho(3);
    cmat sigma = randrho(3);
    EXPECT_NEAR(0, lognegativity(kron(rho, sigma), {3, 3}), 1e-5);

    // additivity (2 ququads)
    rho = randrho(4);
    sigma = randrho(4);
    EXPECT_NEAR(
        lognegativity(syspermute(kron(rho, sigma), {0, 2, 1, 3}), {4, 4}),
        lognegativity(rho, {2, 2}) + lognegativity(sigma, {2, 2}), 1e-5);

    // random maximally entangled 2-qutrit state
    idx d = 3;
    rho = prj(kron(randU(d), randU(d)) * st.mes(d));
    EXPECT_NEAR(std::log2(d), lognegativity(rho, {d, d}), 1e-5);

    // random maximally entangled 2-ququad state
    d = 4;
    rho = prj(kron(randU(d), randU(d)) * st.mes(d));
    EXPECT_NEAR(std::log2(d), lognegativity(rho, {d, d}), 1e-5);

    // random maximally entangled state (d = 7)
    d = 7;
    rho = prj(kron(randU(d), randU(d)) * st.mes(d));
    EXPECT_NEAR(std::log2(d), lognegativity(rho, {d, d}), 1e-5);
}
/******************************************************************************/
/// BEGIN template <typename Derived> realT lognegativity(
///       const Eigen::MatrixBase<Derived>& A, idx d = 2)
TEST(qpp_lognegativity, Qubits) {
    // zero on product states (2 qubits)
    cmat rho = randrho();
    cmat sigma = randrho();
    EXPECT_NEAR(0, lognegativity(kron(rho, sigma)), 1e-5);

    // random 2-qubit state with Schmidt coefficients 0.8 and 0.6
    ket psi =
        kron(randU(), randU()) * (0.8 * mket({0, 0}) + 0.6 * mket({1, 1}));
    rho = prj(psi);
    EXPECT_NEAR(0.9708536, lognegativity(rho), 1e-5);

    // random maximally entangled 2-qubit state
    idx d = 2;
    rho = prj(kron(randU(d), randU(d)) * st.mes(d));
    EXPECT_NEAR(std::log2(d), lognegativity(rho), 1e-5);
}
/******************************************************************************/
/// BEGIN template <typename Derived> realT negativity(
///       const Eigen::MatrixBase<Derived>& A, const std::vector<idx>& dims)
TEST(qpp_negativity, Qudits) {
    // Must be (d - 1)/2 on MES

    // zero on product states (2 qutrits)
    cmat rho = randrho(3);
    cmat sigma = randrho(3);
    EXPECT_NEAR(0, negativity(kron(rho, sigma), {3, 3}), 1e-5);

    // convexity (10 ququads)
    idx N = 10;
    idx d = 4;
    std::vector<cmat> rhos(N);
    std::vector<realT> probs = randprob(N);
    rho = cmat::Zero(d * d, d * d);
    realT sum_neg = 0;
    for (idx i = 0; i < N; ++i) {
        cmat rho_i = randrho(d * d);
        sum_neg += probs[i] * negativity(rho_i, {d, d});
        rho += probs[i] * rho_i;
    }
    EXPECT_LE(negativity(rho, {d, d}), sum_neg);

    // random maximally entangled 2-qutrit state
    d = 3;
    rho = prj(kron(randU(d), randU(d)) * st.mes(d));
    EXPECT_NEAR((d - 1) / 2., negativity(rho, {d, d}), 1e-5);

    // random maximally entangled 2-ququad state
    d = 4;
    rho = prj(kron(randU(d), randU(d)) * st.mes(d));
    EXPECT_NEAR((d - 1) / 2., negativity(rho, {d, d}), 1e-5);

    // random maximally entangled state (d = 7)
    d = 7;
    rho = prj(kron(randU(d), randU(d)) * st.mes(d));
    EXPECT_NEAR((d - 1) / 2., negativity(rho, {d, d}), 1e-5);
}
/******************************************************************************/
/// BEGIN template <typename Derived> realT negativity(
///       const Eigen::MatrixBase<Derived>& A, idx d = 2)
TEST(qpp_negativity, Qubits) {
    // Must be (d - 1)/2 on MES

    // zero on product states (2 qutrits)
    cmat rho = randrho();
    cmat sigma = randrho();
    EXPECT_NEAR(0, negativity(kron(rho, sigma)), 1e-5);

    // random 2-qubit state with Schmidt coefficients 0.8 and 0.6
    ket psi =
        kron(randU(), randU()) * (0.8 * mket({0, 0}) + 0.6 * mket({1, 1}));
    rho = prj(psi);
    EXPECT_NEAR(0.48, negativity(rho), 1e-5);

    // random maximally entangled 2-qubit state
    psi = kron(randU(), randU()) * st.mes();
    rho = prj(psi);
    EXPECT_NEAR(0.5, negativity(rho), 1e-5);
}
/******************************************************************************/
/// BEGIN template <typename Derived>
///       std::tuple<cmat, cmat, dyn_col_vect<realT>, dyn_col_vect<realT>>
///       schmidt(const Eigen::MatrixBase<Derived>& A,
///       const std::vector<idx>& dims)
TEST(qpp_schmidt, Qudits) {
    // random 3 x 4 state
    idx dA = 3, dB = 4, D = dA * dB, minD = std::min(dA, dB);
    auto const psi = randket(D);

    auto const t = schmidt(psi, {dA, dB});
    auto const& basisA = std::get<0>(t);
    auto const& basisB = std::get<1>(t);
    auto const& coeffs = std::get<2>(t);
    auto const& probs = std::get<3>(t);

    auto const basisA_ref = schmidtA(psi, {dA, dB});
    EXPECT_NEAR(0, norm(basisA - basisA_ref), 1e-5);

    auto const basisB_ref = schmidtB(psi, {dA, dB});
    EXPECT_NEAR(0, norm(basisB - basisB_ref), 1e-5);

    auto const coeffs_ref = schmidtcoeffs(psi, {dA, dB});
    EXPECT_NEAR(0, norm(coeffs - coeffs_ref), 1e-5);

    auto const probs_ref_vect = schmidtprobs(psi, {dA, dB});
    auto const probs_ref =
        dyn_col_vect<realT>::Map(probs_ref_vect.data(), probs_ref_vect.size());
    EXPECT_NEAR(0, norm(probs - probs_ref), 1e-5);
}
/******************************************************************************/
/// BEGIN template <typename Derived>
///       std::tuple<cmat, cmat, dyn_col_vect<realT>, dyn_col_vect<realT>>
///       schmidt(const Eigen::MatrixBase<Derived>& A, idx d = 2)
TEST(qpp_schmidt, Qubits) {
    // random 5 x 5 state
    idx d = 5, D = d * d;
    auto const psi = randket(D);

    auto const t = schmidt(psi, d);
    auto const& basisA = std::get<0>(t);
    auto const& basisB = std::get<1>(t);
    auto const& coeffs = std::get<2>(t);
    auto const& probs = std::get<3>(t);

    auto const basisA_ref = schmidtA(psi, d);
    EXPECT_NEAR(0, norm(basisA - basisA_ref), 1e-5);

    auto const basisB_ref = schmidtB(psi, d);
    EXPECT_NEAR(0, norm(basisB - basisB_ref), 1e-5);

    auto const coeffs_ref = schmidtcoeffs(psi, d);
    EXPECT_NEAR(0, norm(coeffs - coeffs_ref), 1e-5);

    auto const probs_ref_vect = schmidtprobs(psi, d);
    auto const probs_ref =
        dyn_col_vect<realT>::Map(probs_ref_vect.data(), probs_ref_vect.size());
    EXPECT_NEAR(0, norm(probs - probs_ref), 1e-5);
}
/******************************************************************************/
/// BEGIN template <typename Derived> cmat schmidtA/B(
///       const Eigen::MatrixBase<Derived>& A, const std::vector<idx>& dims)
TEST(qpp_schmidtA_schmidtB, Qudits) {
    // random degenerate 1 x 1 product state
    idx dA = 1, dB = 1, D = dA * dB, minD = std::min(dA, dB);
    cmat UA = randU(dA);
    cmat UB = randU(dB);
    ket psi = kron(UA, UB) * mket({0, 0}, {dA, dB});
    cmat basisA = schmidtA(psi, {dA, dB});
    cmat basisB = schmidtB(psi, {dA, dB});
    // unitarity
    EXPECT_NEAR(0, norm(adjoint(basisA) * basisA - gt.Id(dA)), 1e-5);
    EXPECT_NEAR(0, norm(adjoint(basisB) * basisB - gt.Id(dB)), 1e-5);
    // get the Schmidt coefficients and test the result
    dyn_col_vect<realT> scf = schmidtcoeffs(psi, {dA, dB});
    ket expected = ket::Zero(D);
    for (idx i = 0; i < minD; ++i) {
        expected += scf(i) * kron(basisA.col(i), basisB.col(i));
    }
    EXPECT_NEAR(0, norm(expected - psi), 1e-5);

    // random degenerate 3 x 1 product state
    dA = 3, dB = 1, D = dA * dB, minD = std::min(dA, dB);
    UA = randU(dA);
    UB = randU(dB);
    psi = kron(UA, UB) * mket({0, 0}, {dA, dB});
    basisA = schmidtA(psi, {dA, dB});
    basisB = schmidtB(psi, {dA, dB});
    // unitarity
    EXPECT_NEAR(0, norm(adjoint(basisA) * basisA - gt.Id(dA)), 1e-5);
    EXPECT_NEAR(0, norm(adjoint(basisB) * basisB - gt.Id(dB)), 1e-5);
    // get the Schmidt coefficients and test the result
    scf = schmidtcoeffs(psi, {dA, dB});
    expected = ket::Zero(D);
    for (idx i = 0; i < minD; ++i) {
        expected += scf(i) * kron(basisA.col(i), basisB.col(i));
    }
    EXPECT_NEAR(0, norm(expected - psi), 1e-5);

    // random 2 x 4 product state
    dA = 2, dB = 4, D = dA * dB, minD = std::min(dA, dB);
    UA = randU(dA);
    UB = randU(dB);
    psi = kron(UA, UB) * mket({0, 0}, {dA, dB});
    basisA = schmidtA(psi, {dA, dB});
    basisB = schmidtB(psi, {dA, dB});
    // unitarity
    EXPECT_NEAR(0, norm(adjoint(basisA) * basisA - gt.Id(dA)), 1e-5);
    EXPECT_NEAR(0, norm(adjoint(basisB) * basisB - gt.Id(dB)), 1e-5);
    // get the Schmidt coefficients and test the result
    scf = schmidtcoeffs(psi, {dA, dB});
    expected = ket::Zero(D);
    for (idx i = 0; i < minD; ++i) {
        expected += scf(i) * kron(basisA.col(i), basisB.col(i));
    }
    EXPECT_NEAR(0, norm(expected - psi), 1e-5);

    // random 3 x 4 state
    dA = 3, dB = 4, D = dA * dB, minD = std::min(dA, dB);
    psi = randket(D);
    basisA = schmidtA(psi, {dA, dB});
    basisB = schmidtB(psi, {dA, dB});
    // unitarity
    EXPECT_NEAR(0, norm(adjoint(basisA) * basisA - gt.Id(dA)), 1e-5);
    EXPECT_NEAR(0, norm(adjoint(basisB) * basisB - gt.Id(dB)), 1e-5);
    // get the Schmidt coefficients and test the result
    scf = schmidtcoeffs(psi, {dA, dB});
    expected = ket::Zero(D);
    for (idx i = 0; i < minD; ++i) {
        expected += scf(i) * kron(basisA.col(i), basisB.col(i));
    }
    EXPECT_NEAR(0, norm(expected - psi), 1e-5);
}
/******************************************************************************/
/// BEGIN template <typename Derived> cmat schmidtA/B(
///       const Eigen::MatrixBase<Derived>& A, idx d = 2)
TEST(qpp_schmidtA_schmidtB, Qubits) {
    // random 2 x 2 product state
    idx d = 2, D = d * d;
    cmat UA = randU(d);
    cmat UB = randU(d);
    ket psi = kron(UA, UB) * st.zero(2, d);
    cmat basisA = schmidtA(psi, d);
    cmat basisB = schmidtB(psi, d);
    // unitarity
    EXPECT_NEAR(0, norm(adjoint(basisA) * basisA - gt.Id(d)), 1e-5);
    EXPECT_NEAR(0, norm(adjoint(basisB) * basisB - gt.Id(d)), 1e-5);
    // get the Schmidt coefficients and test the result
    dyn_col_vect<realT> scf = schmidtcoeffs(psi, d);
    ket expected = ket::Zero(D);
    for (idx i = 0; i < d; ++i) {
        expected += scf(i) * kron(basisA.col(i), basisB.col(i));
    }
    EXPECT_NEAR(0, norm(expected - psi), 1e-5);

    // random 4 x 4 product state
    d = 4, D = d * d;
    UA = randU(d);
    UB = randU(d);
    psi = kron(UA, UB) * st.zero(2, d);
    basisA = schmidtA(psi, d);
    basisB = schmidtB(psi, d);
    // unitarity
    EXPECT_NEAR(0, norm(adjoint(basisA) * basisA - gt.Id(d)), 1e-5);
    EXPECT_NEAR(0, norm(adjoint(basisB) * basisB - gt.Id(d)), 1e-5);
    // get the Schmidt coefficients and test the result
    scf = schmidtcoeffs(psi, d);
    expected = ket::Zero(D);
    for (idx i = 0; i < d; ++i) {
        expected += scf(i) * kron(basisA.col(i), basisB.col(i));
    }
    EXPECT_NEAR(0, norm(expected - psi), 1e-5);

    // random 5 x 5 state
    d = 5, D = d * d;
    psi = randket(D);
    basisA = schmidtA(psi, d);
    basisB = schmidtB(psi, d);
    // unitarity
    EXPECT_NEAR(0, norm(adjoint(basisA) * basisA - gt.Id(d)), 1e-5);
    EXPECT_NEAR(0, norm(adjoint(basisB) * basisB - gt.Id(d)), 1e-5);
    // get the Schmidt coefficients and test the result
    scf = schmidtcoeffs(psi, d);
    expected = ket::Zero(D);
    for (idx i = 0; i < d; ++i) {
        expected += scf(i) * kron(basisA.col(i), basisB.col(i));
    }
    EXPECT_NEAR(0, norm(expected - psi), 1e-5);
}
/******************************************************************************/
/// BEGIN template <typename Derived> dyn_col_vect<realT> schmidtcoeffs(
///       const Eigen::MatrixBase<Derived>& A, const std::vector<idx>& dims)
TEST(qpp_schmidtcoeffs, Qudits) {
    // random degenerate 1 x 1 product state
    idx dA = 1, dB = 1, minD = std::min(dA, dB);
    cmat UA = randU(dA);
    cmat UB = randU(dB);
    ket psi = kron(UA, UB) * mket({0, 0}, {dA, dB});
    dyn_col_vect<realT> result = schmidtcoeffs(psi, {dA, dB});
    dyn_col_vect<realT> expected(minD);
    expected << 1;
    EXPECT_NEAR(0, norm(result - expected), 1e-5);

    // random degenerate 3 x 1 product state
    dA = 3, dB = 1, minD = std::min(dA, dB);
    UA = randU(dA);
    UB = randU(dB);
    psi = kron(UA, UB) * mket({0, 0}, {dA, dB});
    result = schmidtcoeffs(psi, {dA, dB});
    expected = dyn_col_vect<realT>(minD);
    expected << 1;
    EXPECT_NEAR(0, norm(result - expected), 1e-5);

    // random 3 x 4 product state
    dA = 3, dB = 4, minD = std::min(dA, dB);
    UA = randU(dA);
    UB = randU(dB);
    psi = kron(UA, UB) * mket({0, 0}, {dA, dB});
    result = schmidtcoeffs(psi, {dA, dB});
    expected = dyn_col_vect<realT>::Zero(minD);
    expected(0) = 1;
    EXPECT_NEAR(0, norm(result - expected), 1e-5);

    // random 2 x 3 state with fixed Schmidt coefficients
    dA = 2, dB = 3, minD = std::min(dA, dB);
    realT c0 = 0.8, c1 = 0.6;
    UA = randU(dA);
    UB = randU(dB);
    psi = kron(UA, UB) *
          (c0 * mket({0, 0}, {dA, dB}) + c1 * mket({1, 1}, {dA, dB}));
    result = schmidtcoeffs(psi, {dA, dB});
    expected = dyn_col_vect<realT>::Zero(minD);
    expected << c0, c1;
    EXPECT_NEAR(0, norm(result - expected), 1e-5);

    // random 5 x 3 state with fixed Schmidt coefficients
    dA = 5, dB = 3, minD = std::min(dA, dB);
    c0 = 0.8, c1 = 0.5;
    realT c2 = std::sqrt(1 - c0 * c0 - c1 * c1);
    UA = randU(dA);
    UB = randU(dB);
    psi = kron(UA, UB) *
          (c0 * mket({0, 0}, {dA, dB}) + c1 * mket({1, 1}, {dA, dB}) +
           c2 * mket({2, 2}, {dA, dB}));
    result = schmidtcoeffs(psi, {dA, dB});
    expected = dyn_col_vect<realT>::Zero(minD);
    expected << c0, c1, c2;
    EXPECT_NEAR(0, norm(result - expected), 1e-5);
}
/******************************************************************************/
/// BEGIN template <typename Derived> dyn_col_vect<realT> schmidtcoeffs(
///       const Eigen::MatrixBase<Derived>& A, idx d = 2)
TEST(qpp_schmidtcoeffs, Qubits) {
    // random 2 x 2 product state
    idx d = 2;
    cmat UA = randU(d);
    cmat UB = randU(d);
    ket psi = kron(UA, UB) * mket({0, 0});
    dyn_col_vect<realT> result = schmidtcoeffs(psi);
    dyn_col_vect<realT> expected = dyn_col_vect<realT>::Zero(d);
    expected(0) = 1;
    EXPECT_NEAR(0, norm(result - expected), 1e-5);

    // random 3 x 3 product state
    d = 3;
    UA = randU(d);
    UB = randU(d);
    psi = kron(UA, UB) * mket({1, 1}, d);
    result = schmidtcoeffs(psi, d);
    expected = dyn_col_vect<realT>::Zero(d);
    expected(0) = 1;
    EXPECT_NEAR(0, norm(result - expected), 1e-5);

    // random 2 x 2 state with fixed Schmidt coefficients
    d = 2;
    realT c0 = 0.8, c1 = 0.6;
    UA = randU(d);
    UB = randU(d);
    psi = kron(UA, UB) * (c0 * st.zero(2) + c1 * st.one(2));
    result = schmidtcoeffs(psi);
    expected = dyn_col_vect<realT>::Zero(d);
    expected << c0, c1;
    EXPECT_NEAR(0, norm(result - expected), 1e-5);

    // random 3 x 3 state with fixed Schmidt coefficients
    d = 3;
    c0 = 0.8, c1 = 0.5;
    realT c2 = std::sqrt(1 - c0 * c0 - c1 * c1);
    UA = randU(d);
    UB = randU(d);
    psi = kron(UA, UB) *
          (c0 * mket({0, 0}, d) + c1 * mket({1, 1}, d) + c2 * mket({2, 2}, d));
    result = schmidtcoeffs(psi, d);
    expected = dyn_col_vect<realT>::Zero(d);
    expected << c0, c1, c2;
    EXPECT_NEAR(0, norm(result - expected), 1e-5);
}
/******************************************************************************/
/// BEGIN template <typename Derived> std::vector<realT> schmidtprobs(
///       const Eigen::MatrixBase<Derived>& A, const std::vector<idx>& dims)
TEST(qpp_schmidtprobs, Qudits) {
    // random degenerate 1 x 1 product state
    idx dA = 1, dB = 1, minD = std::min(dA, dB);
    cmat UA = randU(dA);
    cmat UB = randU(dB);
    ket psi = kron(UA, UB) * mket({0, 0}, {dA, dB});
    std::vector<realT> result_vect = schmidtprobs(psi, {dA, dB});
    dyn_col_vect<realT> result =
        Eigen::Map<dyn_col_vect<realT>>(result_vect.data(), result_vect.size());
    dyn_col_vect<realT> expected(minD);
    expected << 1;
    EXPECT_NEAR(0, norm(result - expected), 1e-5);

    // random degenerate 3 x 1 product state
    dA = 3, dB = 1, minD = std::min(dA, dB);
    UA = randU(dA);
    UB = randU(dB);
    psi = kron(UA, UB) * mket({0, 0}, {dA, dB});
    result_vect = schmidtprobs(psi, {dA, dB});
    result =
        Eigen::Map<dyn_col_vect<realT>>(result_vect.data(), result_vect.size());
    expected = dyn_col_vect<realT>(minD);
    expected << 1;
    EXPECT_NEAR(0, norm(result - expected), 1e-5);

    // random 3 x 4 product state
    dA = 3, dB = 4, minD = std::min(dA, dB);
    UA = randU(dA);
    UB = randU(dB);
    psi = kron(UA, UB) * mket({0, 0}, {dA, dB});
    result_vect = schmidtprobs(psi, {dA, dB});
    result =
        Eigen::Map<dyn_col_vect<realT>>(result_vect.data(), result_vect.size());
    expected = dyn_col_vect<realT>::Zero(minD);
    expected(0) = 1;
    EXPECT_NEAR(0, norm(result - expected), 1e-5);

    // random 2 x 3 state with fixed Schmidt coefficients
    dA = 2, dB = 3, minD = std::min(dA, dB);
    realT c0 = 0.8, c1 = 0.6;
    UA = randU(dA);
    UB = randU(dB);
    psi = kron(UA, UB) *
          (c0 * mket({0, 0}, {dA, dB}) + c1 * mket({1, 1}, {dA, dB}));
    result_vect = schmidtprobs(psi, {dA, dB});
    result =
        Eigen::Map<dyn_col_vect<realT>>(result_vect.data(), result_vect.size());
    expected = dyn_col_vect<realT>::Zero(minD);
    expected << c0 * c0, c1 * c1;
    EXPECT_NEAR(0, norm(result - expected), 1e-5);

    // random 5 x 3 state with fixed Schmidt coefficients
    dA = 5, dB = 3, minD = std::min(dA, dB);
    c0 = 0.8, c1 = 0.5;
    realT c2 = std::sqrt(1 - c0 * c0 - c1 * c1);
    UA = randU(dA);
    UB = randU(dB);
    psi = kron(UA, UB) *
          (c0 * mket({0, 0}, {dA, dB}) + c1 * mket({1, 1}, {dA, dB}) +
           c2 * mket({2, 2}, {dA, dB}));
    result_vect = schmidtprobs(psi, {dA, dB});
    result =
        Eigen::Map<dyn_col_vect<realT>>(result_vect.data(), result_vect.size());
    expected = dyn_col_vect<realT>::Zero(minD);
    expected << c0 * c0, c1 * c1, c2 * c2;
    EXPECT_NEAR(0, norm(result - expected), 1e-5);
}
/******************************************************************************/
/// BEGIN template <typename Derived> std::vector<realT> schmidtprobs(
///       const Eigen::MatrixBase<Derived>& A, idx d = 2)
TEST(qpp_schmidtprobs, Qubits) {
    // random 2 x 2 product state
    idx d = 2;
    cmat UA = randU(d);
    cmat UB = randU(d);
    ket psi = kron(UA, UB) * mket({0, 0});
    std::vector<realT> result_vect = schmidtprobs(psi);
    dyn_col_vect<realT> result =
        Eigen::Map<dyn_col_vect<realT>>(result_vect.data(), result_vect.size());
    dyn_col_vect<realT> expected = dyn_col_vect<realT>::Zero(d);
    expected(0) = 1;
    EXPECT_NEAR(0, norm(result - expected), 1e-5);

    // random 3 x 3 product state
    d = 3;
    UA = randU(d);
    UB = randU(d);
    psi = kron(UA, UB) * mket({1, 1}, d);
    result_vect = schmidtprobs(psi, d);
    result =
        Eigen::Map<dyn_col_vect<realT>>(result_vect.data(), result_vect.size());
    expected = dyn_col_vect<realT>::Zero(d);
    expected(0) = 1;
    EXPECT_NEAR(0, norm(result - expected), 1e-5);

    // random 2 x 2 state with fixed Schmidt coefficients
    d = 2;
    realT c0 = 0.8, c1 = 0.6;
    UA = randU(d);
    UB = randU(d);
    psi = kron(UA, UB) * (c0 * st.zero(2) + c1 * st.one(2));
    result_vect = schmidtprobs(psi);
    result =
        Eigen::Map<dyn_col_vect<realT>>(result_vect.data(), result_vect.size());
    expected = dyn_col_vect<realT>::Zero(d);
    expected << c0 * c0, c1 * c1;
    EXPECT_NEAR(0, norm(result - expected), 1e-5);

    // random 3 x 3 state with fixed Schmidt coefficients
    d = 3;
    c0 = 0.8, c1 = 0.5;
    realT c2 = std::sqrt(1 - c0 * c0 - c1 * c1);
    UA = randU(d);
    UB = randU(d);
    psi = kron(UA, UB) *
          (c0 * mket({0, 0}, d) + c1 * mket({1, 1}, d) + c2 * mket({2, 2}, d));
    result_vect = schmidtprobs(psi, d);
    result =
        Eigen::Map<dyn_col_vect<realT>>(result_vect.data(), result_vect.size());
    expected = dyn_col_vect<realT>::Zero(d);
    expected << c0 * c0, c1 * c1, c2 * c2;
    EXPECT_NEAR(0, norm(result - expected), 1e-5);
}
/******************************************************************************/
