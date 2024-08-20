#include <cmath>

#include "gtest/gtest.h"

#include "qpp/qpp.hpp"

using namespace qpp;

// Unit testing "entropies.hpp"

/// BEGIN template <typename Derived> realT entropy(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_entropy, Matrix) {
    // 1 x 1 case
    cmat A(1, 1);
    A << 1.;
    EXPECT_NEAR(0, entropy(A), 1e-5);

    // 2 x 2 random matrix with fixed only 1 non-zero eigenvalue
    idx D = 2;
    cmat evals = cmat::Zero(D, 1);
    evals << 1., 0.;
    A = evals.asDiagonal();
    cmat U = randU(D);
    EXPECT_NEAR(0, entropy(A), 1e-5);

    // 2 x 2 random matrix with fixed eigenvalues
    D = 2;
    evals = cmat::Zero(D, 1);
    evals << 0.6, 0.4;
    A = evals.asDiagonal();
    U = randU(D);
    EXPECT_NEAR(0.970950594455, entropy(A), 1e-5);

    // 2 x 2 random matrix with fixed equal eigenvalues
    D = 2;
    evals = cmat::Zero(D, 1);
    evals << 0.5, 0.5;
    A = evals.asDiagonal();
    U = randU(D);
    EXPECT_NEAR(1, entropy(A), 1e-5);

    // 3 x 3 random matrix with fixed equal eigenvalues
    D = 3;
    evals = cmat::Zero(D, 1);
    evals << 1 / 3., 1 / 3., 1 / 3.;
    A = evals.asDiagonal();
    U = randU(D);
    EXPECT_NEAR(std::log2(3), entropy(A), 1e-5);
}

/// BEGIN inline realT entropy(const std::vector<realT>& prob)
TEST(qpp_entropy, Vector) {
    // 1 value
    std::vector<realT> v = {1};
    EXPECT_NEAR(0, entropy(v), 1e-5);

    // 2 values, only 1 non-zero
    v = {1, 0};
    EXPECT_NEAR(0, entropy(v), 1e-5);

    // 2 fixed values
    v = {0.6, 0.4};
    EXPECT_NEAR(0.970950594455, entropy(v), 1e-5);

    // 2 equal values
    v = {0.5, 0.5};
    EXPECT_NEAR(1, entropy(v), 1e-5);

    // 3 equal values
    v = {1 / 3., 1 / 3., 1 / 3.};
    idx D = v.size();
    EXPECT_NEAR(std::log2(D), entropy(v), 1e-5);
}

/// BEGIN template <typename Derived> realT qmutualinfo(
///       const Eigen::MatrixBase<Derived>& A, const std::vector<idx>& subsysA,
///       const std::vector<idx>& subsysB, const std::vector<idx>& dims)
TEST(qpp_qmutualinfo, Qudits) {
    // 1 x 1 degenerate product state
    idx dA = 1, dB = 1;
    cmat rhoA = randrho(dA), rhoB = randrho(dB);
    cmat rho = kron(rhoA, rhoB);
    realT result = qmutualinfo(rho, {0}, {1}, {dA, dB});
    realT expected = 0;
    EXPECT_NEAR(result, expected, 1e-5);

    // 2 x 3 product state
    dA = 2, dB = 3;
    rhoA = randrho(dA), rhoB = randrho(dB);
    rho = kron(rhoA, rhoB);
    result = qmutualinfo(rho, {0}, {1}, {dA, dB});
    expected = 0;
    EXPECT_NEAR(result, expected, 1e-5);

    // 3 x 3 maximally entangled state
    dA = dB = 3;
    rho = prj(st.mes(dA));
    result = qmutualinfo(rho, {0}, {1}, {dA, dB});
    expected = 2 * std::log2(dA);
    EXPECT_NEAR(result, expected, 1e-5);

    // 3 x 3 maximally entangled state tensored with a 2 x 2 product state
    rho = prj(st.mes(3)); // MES
    rho = kron(rho, kron(randrho(2), randrho(2)));
    rho = syspermute(rho, {0, 2, 1, 3}, {3, 3, 2, 2});
    result = qmutualinfo(rho, {0}, {2}, {3, 2, 3, 2});
    expected = 2 * std::log2(3);
    EXPECT_NEAR(result, expected, 1e-5);

    // random 3 x 4 state
    dA = 3, dB = 4;
    rho = randrho(dA * dB);
    rhoA = ptrace2(rho, {dA, dB});
    rhoB = ptrace1(rho, {dA, dB});
    result = qmutualinfo(rho, {0}, {1}, {dA, dB});
    expected = entropy(rhoA) + entropy(rhoB) - entropy(rho);
    EXPECT_NEAR(result, expected, 1e-5);
}

/// BEGIN template <typename Derived> realT qmutualinfo(
///       const Eigen::MatrixBase<Derived>& A, const std::vector<idx>& subsysA,
///       const std::vector<idx>& subsysB, idx d = 2)
TEST(qpp_qmutualinfo, Qubits) {
    // 2 x 2 product state
    idx d = 2;
    cmat rhoA = randrho(d), rhoB = randrho(d);
    cmat rho = kron(rhoA, rhoB);
    realT result = qmutualinfo(rho, {0}, {1});
    realT expected = 0;
    EXPECT_NEAR(result, expected, 1e-5);

    // 3 x 3 maximally entangled state
    d = 3;
    rho = prj(st.mes(d));
    result = qmutualinfo(rho, {0}, {1}, d);
    expected = 2 * std::log2(d);
    EXPECT_NEAR(result, expected, 1e-5);

    // 3 x 3 maximally entangled state tensored with a 3 x 3 product state
    rho = prj(st.mes(3)); // MES
    rho = kron(rho, kron(randrho(3), randrho(3)));
    rho = syspermute(rho, {0, 2, 1, 3}, {3, 3, 3, 3});
    result = qmutualinfo(rho, {0}, {2}, 3);
    expected = 2 * std::log2(3);
    EXPECT_NEAR(result, expected, 1e-5);

    // random 3 x 3 state
    d = 3;
    rho = randrho(d * d);
    rhoA = ptrace2(rho, d);
    rhoB = ptrace1(rho, d);
    result = qmutualinfo(rho, {0}, {1}, d);
    expected = entropy(rhoA) + entropy(rhoB) - entropy(rho);
    EXPECT_NEAR(result, expected, 1e-5);
}

/// BEGIN template <typename Derived> realT renyi(
///       const Eigen::MatrixBase<Derived>& A, realT alpha)
TEST(qpp_renyi, Matrix) {
    // 1 x 1 case
    cmat A(1, 1);
    A << 1.;
    EXPECT_NEAR(0, renyi(A, 0), 1e-5);
    EXPECT_NEAR(0, renyi(A, 1 / 2.), 1e-5);
    EXPECT_NEAR(0, renyi(A, 1), 1e-5);
    EXPECT_NEAR(0, renyi(A, 2), 1e-5);
    EXPECT_NEAR(0, renyi(A, infty), 1e-5);

    // 2 x 2 random matrix with fixed eigenvalues
    idx D = 2;
    cmat evals = cmat::Zero(D, 1);
    evals << 0.6, 0.4;
    A = evals.asDiagonal();
    cmat U = randU(D);
    EXPECT_NEAR(1, renyi(A, 0), 1e-5);
    EXPECT_NEAR(0.985351706365, renyi(A, 1 / 2.), 1e-5);
    EXPECT_NEAR(0.970950594455, renyi(A, 1), 1e-5);
    EXPECT_NEAR(0.943416471634, renyi(A, 2), 1e-5);
    EXPECT_NEAR(0.918250633859, renyi(A, 3), 1e-5);
    EXPECT_NEAR(-std::log2(0.6), renyi(A, infty), 1e-5);

    // 2 x 2 random matrix with fixed equal eigenvalues
    D = 2;
    evals = cmat::Zero(D, 1);
    evals << 0.5, 0.5;
    A = evals.asDiagonal();
    U = randU(D);
    EXPECT_NEAR(1, renyi(A, 0), 1e-5);
    EXPECT_NEAR(1, renyi(A, 1 / 2.), 1e-5);
    EXPECT_NEAR(1, renyi(A, 1), 1e-5);
    EXPECT_NEAR(1, renyi(A, 2), 1e-5);
    EXPECT_NEAR(1, renyi(A, 3), 1e-5);
    EXPECT_NEAR(1, renyi(A, infty), 1e-5);

    // 3 x 3 random matrix with fixed only 1 non-zero eigenvalue
    D = 3;
    evals = cmat::Zero(D, 1);
    evals << 1., 0., 0.;
    A = evals.asDiagonal();
    U = randU(D);
    EXPECT_NEAR(std::log2(D), renyi(A, 0), 1e-5);
    EXPECT_NEAR(0, renyi(A, 1 / 2.), 1e-5);
    EXPECT_NEAR(0, renyi(A, 1), 1e-5);
    EXPECT_NEAR(0, renyi(A, 2), 1e-5);
    EXPECT_NEAR(0, renyi(A, infty), 1e-5);

    // 3 x 3 random matrix with fixed equal eigenvalues
    D = 3;
    evals = cmat::Zero(D, 1);
    evals << 1 / 3., 1 / 3., 1 / 3.;
    A = evals.asDiagonal();
    U = randU(D);
    EXPECT_NEAR(std::log2(D), renyi(A, 0), 1e-5);
    EXPECT_NEAR(std::log2(D), renyi(A, 1 / 2.), 1e-5);
    EXPECT_NEAR(std::log2(D), renyi(A, 1), 1e-5);
    EXPECT_NEAR(std::log2(D), renyi(A, 2), 1e-5);
    EXPECT_NEAR(std::log2(D), renyi(A, 3), 1e-5);
    EXPECT_NEAR(std::log2(D), renyi(A, infty), 1e-5);
}

/// BEGIN inline realT renyi(const std::vector<realT>& prob, realT alpha)
TEST(qpp_renyi, Vector) {
    // 1 value
    std::vector<realT> v = {1};
    EXPECT_NEAR(0, renyi(v, 0), 1e-5);
    EXPECT_NEAR(0, renyi(v, 1 / 2.), 1e-5);
    EXPECT_NEAR(0, renyi(v, 1), 1e-5);
    EXPECT_NEAR(0, renyi(v, 2), 1e-5);
    EXPECT_NEAR(0, renyi(v, infty), 1e-5);

    // 2 fixed values
    v = {0.6, 0.4};
    EXPECT_NEAR(1, renyi(v, 0), 1e-5);
    EXPECT_NEAR(0.985351706365, renyi(v, 1 / 2.), 1e-5);
    EXPECT_NEAR(0.970950594455, renyi(v, 1), 1e-5);
    EXPECT_NEAR(0.943416471634, renyi(v, 2), 1e-5);
    EXPECT_NEAR(0.918250633859, renyi(v, 3), 1e-5);
    EXPECT_NEAR(-std::log2(0.6), renyi(v, infty), 1e-5);

    // 2 equal values
    v = {0.5, 0.5};
    EXPECT_NEAR(1, renyi(v, 0), 1e-5);
    EXPECT_NEAR(1, renyi(v, 1 / 2.), 1e-5);
    EXPECT_NEAR(1, renyi(v, 1), 1e-5);
    EXPECT_NEAR(1, renyi(v, 2), 1e-5);
    EXPECT_NEAR(1, renyi(v, 3), 1e-5);
    EXPECT_NEAR(1, renyi(v, infty), 1e-5);

    // 3 values, only 1 non-zero
    v = {1, 0, 0};
    idx D = 3;
    EXPECT_NEAR(std::log2(D), renyi(v, 0), 1e-5);
    EXPECT_NEAR(0, renyi(v, 1 / 2.), 1e-5);
    EXPECT_NEAR(0, renyi(v, 1), 1e-5);
    EXPECT_NEAR(0, renyi(v, 2), 1e-5);
    EXPECT_NEAR(0, renyi(v, infty), 1e-5);

    // 3 equal values
    v = {1 / 3., 1 / 3., 1 / 3.};
    EXPECT_NEAR(std::log2(D), renyi(v, 0), 1e-5);
    EXPECT_NEAR(std::log2(D), renyi(v, 1 / 2.), 1e-5);
    EXPECT_NEAR(std::log2(D), renyi(v, 1), 1e-5);
    EXPECT_NEAR(std::log2(D), renyi(v, 2), 1e-5);
    EXPECT_NEAR(std::log2(D), renyi(v, 3), 1e-5);
    EXPECT_NEAR(std::log2(D), renyi(v, infty), 1e-5);
}

/// BEGIN template <typename Derived> realT tsallis(
///       const Eigen::MatrixBase<Derived>& A, realT q)
TEST(qpp_tsallis, Matrix) {
    // 1 x 1 case
    cmat A(1, 1);
    A << 1.;
    EXPECT_NEAR(0, tsallis(A, 0), 1e-5);
    EXPECT_NEAR(0, tsallis(A, 1 / 2.), 1e-5);
    EXPECT_NEAR(0, tsallis(A, 1), 1e-5);
    EXPECT_NEAR(0, tsallis(A, 2), 1e-5);
    EXPECT_NEAR(0, tsallis(A, infty), 1e-5);

    // 2 x 2 random matrix with fixed eigenvalues
    idx D = 2;
    cmat evals = cmat::Zero(D, 1);
    evals << 0.6, 0.4;
    A = evals.asDiagonal();
    cmat U = randU(D);
    EXPECT_NEAR(1, tsallis(A, 0), 1e-5);
    EXPECT_NEAR(0.81410440255, tsallis(A, 1 / 2.), 1e-5);
    EXPECT_NEAR(entropy(A) * std::log(2), tsallis(A, 1), 1e-5);
    EXPECT_NEAR(0.48, tsallis(A, 2), 1e-5);
    EXPECT_NEAR(0.36, tsallis(A, 3), 1e-5);
    EXPECT_NEAR(0, tsallis(A, infty), 1e-5);

    // 2 x 2 random matrix with fixed equal eigenvalues
    D = 2;
    evals = cmat::Zero(D, 1);
    evals << 0.5, 0.5;
    A = evals.asDiagonal();
    U = randU(D);
    EXPECT_NEAR(1, tsallis(A, 0), 1e-5);
    EXPECT_NEAR(0.828427124746, tsallis(A, 1 / 2.), 1e-5);
    EXPECT_NEAR(entropy(A) * std::log(2), tsallis(A, 1), 1e-5);
    EXPECT_NEAR(0.5, tsallis(A, 2), 1e-5);
    EXPECT_NEAR(0.375, tsallis(A, 3), 1e-5);
    EXPECT_NEAR(0, tsallis(A, infty), 1e-5);

    // 3 x 3 random matrix with fixed only 1 non-zero eigenvalue
    D = 3;
    evals = cmat::Zero(D, 1);
    evals << 1., 0., 0.;
    A = evals.asDiagonal();
    U = randU(D);
    EXPECT_NEAR(0, tsallis(A, 0), 1e-5);
    EXPECT_NEAR(0, tsallis(A, 1 / 2.), 1e-5);
    EXPECT_NEAR(0, tsallis(A, 1), 1e-5);
    EXPECT_NEAR(0, tsallis(A, 2), 1e-5);
    EXPECT_NEAR(0, tsallis(A, infty), 1e-5);

    // 3 x 3 random matrix with fixed equal eigenvalues
    D = 3;
    evals = cmat::Zero(D, 1);
    evals << 1 / 3., 1 / 3., 1 / 3.;
    A = evals.asDiagonal();
    U = randU(D);
    EXPECT_NEAR(2, tsallis(A, 0), 1e-5);
    EXPECT_NEAR(1.46410161514, tsallis(A, 1 / 2.), 1e-5);
    EXPECT_NEAR(entropy(A) * std::log(2), tsallis(A, 1), 1e-5);
    EXPECT_NEAR(2 / 3., tsallis(A, 2), 1e-5);
    EXPECT_NEAR(4 / 9., tsallis(A, 3), 1e-5);
    EXPECT_NEAR(0, tsallis(A, infty), 1e-5);
}

/// BEGIN inline realT tsallis(const std::vector<realT>& prob, realT q)
TEST(qpp_tsallis, Vector) {
    // 1 value
    std::vector<realT> v = {1};
    EXPECT_NEAR(0, tsallis(v, 0), 1e-5);
    EXPECT_NEAR(0, tsallis(v, 1 / 2.), 1e-5);
    EXPECT_NEAR(0, tsallis(v, 1), 1e-5);
    EXPECT_NEAR(0, tsallis(v, 2), 1e-5);
    EXPECT_NEAR(0, tsallis(v, infty), 1e-5);

    // 2 fixed values
    v = {0.6, 0.4};
    EXPECT_NEAR(1, tsallis(v, 0), 1e-5);
    EXPECT_NEAR(0.81410440255, tsallis(v, 1 / 2.), 1e-5);
    EXPECT_NEAR(entropy(v) * std::log(2), tsallis(v, 1), 1e-5);
    EXPECT_NEAR(0.48, tsallis(v, 2), 1e-5);
    EXPECT_NEAR(0.36, tsallis(v, 3), 1e-5);
    EXPECT_NEAR(0, tsallis(v, infty), 1e-5);

    // 2 equal values
    v = {0.5, 0.5};
    EXPECT_NEAR(1, tsallis(v, 0), 1e-5);
    EXPECT_NEAR(0.828427124746, tsallis(v, 1 / 2.), 1e-5);
    EXPECT_NEAR(entropy(v) * std::log(2), tsallis(v, 1), 1e-5);
    EXPECT_NEAR(0.5, tsallis(v, 2), 1e-5);
    EXPECT_NEAR(0.375, tsallis(v, 3), 1e-5);
    EXPECT_NEAR(0, tsallis(v, infty), 1e-5);

    // 3 values, only 1 non-zero
    v = {1, 0, 0};
    EXPECT_NEAR(0, tsallis(v, 0), 1e-5);
    EXPECT_NEAR(0, tsallis(v, 1 / 2.), 1e-5);
    EXPECT_NEAR(0, tsallis(v, 1), 1e-5);
    EXPECT_NEAR(0, tsallis(v, 2), 1e-5);
    EXPECT_NEAR(0, tsallis(v, infty), 1e-5);

    // 3 equal values
    v = {1 / 3., 1 / 3., 1 / 3.};
    EXPECT_NEAR(2, tsallis(v, 0), 1e-5);
    EXPECT_NEAR(1.46410161514, tsallis(v, 1 / 2.), 1e-5);
    EXPECT_NEAR(entropy(v) * std::log(2), tsallis(v, 1), 1e-5);
    EXPECT_NEAR(2 / 3., tsallis(v, 2), 1e-5);
    EXPECT_NEAR(4 / 9., tsallis(v, 3), 1e-5);
    EXPECT_NEAR(0, tsallis(v, infty), 1e-5);
}
