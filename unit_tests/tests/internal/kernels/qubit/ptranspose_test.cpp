#include "gtest/gtest.h"

#include "qpp/qpp.hpp"

#include "qpp/internal/kernels/qubit/ptranspose.hpp"

#include <cmath>
#include <complex>
#include <sstream>
#include <string>
#include <vector>

using namespace qpp;

namespace {

constexpr double tolerance = 1e-11;

ket make_state(idx n, idx seed = 1) {
    ket result(idx{1} << n);
    for (idx j = 0; j < static_cast<idx>(result.size()); ++j) {
        result(j) = cplx{static_cast<double>((7 * j + 3 * seed) % 19) - 9.0,
                         static_cast<double>((5 * j + 2 * seed) % 17) - 8.0};
    }
    result.normalize();
    return result;
}

cmat make_density(idx n) {
    const ket psi = make_state(n, 1);
    const ket phi = make_state(n, 3);
    return 0.4 * (psi * psi.adjoint()) + 0.6 * (phi * phi.adjoint());
}

cmat make_general_matrix(idx n) {
    const idx D = idx{1} << n;
    cmat result(D, D);
    for (idx r = 0; r < D; ++r) {
        for (idx c = 0; c < D; ++c) {
            result(r, c) =
                cplx{0.05 * (r + 1) - 0.03 * (c + 1),
                     0.02 * (2 * r + c + 1)};
        }
    }
    return result;
}

std::string describe(const std::vector<idx>& target) {
    std::ostringstream out;
    out << "transpose targets={";
    for (idx j = 0; j < target.size(); ++j) {
        if (j != 0) {
            out << ",";
        }
        out << target[j];
    }
    out << "}";
    return out.str();
}

// Independent reference for partial transpose: interchange the row and
// column bit belonging to each selected big-endian qubit.
cmat reference_ptranspose(const cmat& rho,
                          const std::vector<idx>& target, idx n) {
    const idx D = idx{1} << n;
    cmat result(D, D);
    for (idx r = 0; r < D; ++r) {
        for (idx c = 0; c < D; ++c) {
            idx source_row = r;
            idx source_col = c;
            for (idx q : target) {
                const idx mask = idx{1} << (n - 1 - q);
                if ((r & mask) != (c & mask)) {
                    source_row ^= mask;
                    source_col ^= mask;
                }
            }
            result(r, c) = rho(source_row, source_col);
        }
    }
    return result;
}

template <typename Derived1, typename Derived2>
void expect_same(const Eigen::MatrixBase<Derived1>& actual,
                 const Eigen::MatrixBase<Derived2>& expected) {
    ASSERT_EQ(actual.rows(), expected.rows());
    ASSERT_EQ(actual.cols(), expected.cols());
    for (Eigen::Index r = 0; r < actual.rows(); ++r) {
        for (Eigen::Index c = 0; c < actual.cols(); ++c) {
            EXPECT_LE(std::abs(actual(r, c) - expected(r, c)), tolerance)
                << "Mismatch at (" << r << ", " << c << ")";
        }
    }
}

} // namespace

TEST(qpp_ptranspose_psi_kq, EmptyTargetReturnsPureDensityMatrix) {
    for (idx n : {idx{1}, idx{2}, idx{4}}) {
        const ket psi = make_state(n);
        const cmat actual =
            internal::kernels::qubit::ptranspose_psi_kq(psi, {}, n);
        const cmat expected = psi * psi.adjoint();
        expect_same(actual, expected);
    }
}

TEST(qpp_ptranspose_psi_kq, FullTargetReturnsOrdinaryTranspose) {
    for (idx n : {idx{1}, idx{2}, idx{4}}) {
        const ket psi = make_state(n);
        const cmat rho = psi * psi.adjoint();
        std::vector<idx> target;
        for (idx q = 0; q < n; ++q) {
            target.push_back(q);
        }
        const cmat actual =
            internal::kernels::qubit::ptranspose_psi_kq(psi, target, n);
        expect_same(actual, rho.transpose());
    }
}

TEST(qpp_ptranspose_psi_kq, BellStatePartialTranspose) {
    ket psi = ket::Zero(4);
    psi(0) = 1.0 / std::sqrt(2.0);
    psi(3) = 1.0 / std::sqrt(2.0);

    cmat expected = cmat::Zero(4, 4);
    expected(0, 0) = 0.5;
    expected(3, 3) = 0.5;
    expected(1, 2) = 0.5;
    expected(2, 1) = 0.5;

    for (idx q : {idx{0}, idx{1}}) {
        const cmat actual =
            internal::kernels::qubit::ptranspose_psi_kq(psi, {q}, 2);
        expect_same(actual, expected);
    }
}

TEST(qpp_ptranspose_psi_kq, ComplexBellStatePreservesPhase) {
    ket psi = ket::Zero(4);
    psi(0) = 1.0 / std::sqrt(2.0);
    psi(3) = cplx{0, 1.0 / std::sqrt(2.0)};
    const cmat rho = psi * psi.adjoint();
    const cmat expected = reference_ptranspose(rho, {1}, 2);
    const cmat actual =
        internal::kernels::qubit::ptranspose_psi_kq(psi, {1}, 2);

    expect_same(actual, expected);
    EXPECT_NEAR(std::abs(actual(1, 2) - cplx{0, -0.5}), 0.0, tolerance);
    EXPECT_NEAR(std::abs(actual(2, 1) - cplx{0, 0.5}), 0.0, tolerance);
}

TEST(qpp_ptranspose_psi_kq, ArbitrarySubsetsAndTargetOrders) {
    const idx n = 5;
    const ket psi = make_state(n);
    const ket before = psi;
    const cmat rho = psi * psi.adjoint();
    const std::vector<std::vector<idx>> targets{
        {}, {0}, {4}, {1, 3}, {3, 1}, {0, 2, 4},
        {4, 0, 2}, {3, 0, 4, 1}, {0, 1, 2, 3, 4}};
    for (const auto& target : targets) {
        SCOPED_TRACE(describe(target));
        const cmat actual =
            internal::kernels::qubit::ptranspose_psi_kq(psi, target, n);
        expect_same(actual, reference_ptranspose(rho, target, n));
        EXPECT_NEAR(actual.trace().real(), 1.0, tolerance);
        expect_same(actual, actual.adjoint());
        expect_same(psi, before);
    }
}

TEST(qpp_ptranspose_psi_kq, ProductStatePartialTranspose) {
    // |0> tensor (|0>+i|1>)/sqrt(2), transpose only the second qubit.
    ket psi = ket::Zero(4);
    psi(0) = 1.0 / std::sqrt(2.0);
    psi(1) = cplx{0, 1.0 / std::sqrt(2.0)};
    const cmat rho = psi * psi.adjoint();
    const cmat actual =
        internal::kernels::qubit::ptranspose_psi_kq(psi, {1}, 2);
    expect_same(actual, reference_ptranspose(rho, {1}, 2));
}

TEST(qpp_ptranspose_rho_kq, EmptyTargetReturnsInput) {
    for (idx n : {idx{1}, idx{2}, idx{4}}) {
        const cmat rho = make_general_matrix(n);
        const cmat actual =
            internal::kernels::qubit::ptranspose_rho_kq(rho, {}, n);
        expect_same(actual, rho);
    }
}

TEST(qpp_ptranspose_rho_kq, FullTargetReturnsOrdinaryTranspose) {
    for (idx n : {idx{1}, idx{2}, idx{4}}) {
        const cmat rho = make_general_matrix(n);
        std::vector<idx> target;
        for (idx q = 0; q < n; ++q) {
            target.push_back(q);
        }
        const cmat actual =
            internal::kernels::qubit::ptranspose_rho_kq(rho, target, n);
        expect_same(actual, rho.transpose());
    }
}

TEST(qpp_ptranspose_rho_kq, BellStateHasNegativeEigenvalue) {
    ket psi = ket::Zero(4);
    psi(0) = 1.0 / std::sqrt(2.0);
    psi(3) = 1.0 / std::sqrt(2.0);
    const cmat rho = psi * psi.adjoint();
    const cmat actual =
        internal::kernels::qubit::ptranspose_rho_kq(rho, {1}, 2);
    const cmat expected = reference_ptranspose(rho, {1}, 2);
    expect_same(actual, expected);

    Eigen::SelfAdjointEigenSolver<cmat> solver(actual);
    ASSERT_EQ(solver.info(), Eigen::Success);
    EXPECT_NEAR(solver.eigenvalues()(0), -0.5, tolerance);
    EXPECT_NEAR(solver.eigenvalues()(1), 0.5, tolerance);
    EXPECT_NEAR(solver.eigenvalues()(2), 0.5, tolerance);
    EXPECT_NEAR(solver.eigenvalues()(3), 0.5, tolerance);
}

TEST(qpp_ptranspose_rho_kq, MixedStateArbitrarySubsets) {
    const idx n = 5;
    const cmat rho = make_density(n);
    const cmat before = rho;
    const std::vector<std::vector<idx>> targets{
        {}, {0}, {4}, {1, 3}, {3, 1}, {0, 2, 4},
        {4, 0, 2}, {3, 0, 4, 1}, {0, 1, 2, 3, 4}};
    for (const auto& target : targets) {
        SCOPED_TRACE(describe(target));
        const cmat actual =
            internal::kernels::qubit::ptranspose_rho_kq(rho, target, n);
        expect_same(actual, reference_ptranspose(rho, target, n));
        EXPECT_NEAR(actual.trace().real(), 1.0, tolerance);
        expect_same(actual, actual.adjoint());
        expect_same(rho, before);
    }
}

TEST(qpp_ptranspose_rho_kq, GeneralNonHermitianMatrix) {
    const idx n = 4;
    const cmat original = make_general_matrix(n);
    for (const auto& target :
         std::vector<std::vector<idx>>{{0}, {3}, {0, 2}, {2, 0}, {1, 2, 3}}) {
        SCOPED_TRACE(describe(target));
        const cmat actual =
            internal::kernels::qubit::ptranspose_rho_kq(original, target, n);
        expect_same(actual, reference_ptranspose(original, target, n));
    }
}

TEST(qpp_ptranspose_rho_kq, ApplyingTwiceRestoresOriginal) {
    const idx n = 4;
    const cmat original = make_general_matrix(n);
    for (const auto& target :
         std::vector<std::vector<idx>>{{0}, {3}, {1, 3}, {2, 0, 3}}) {
        SCOPED_TRACE(describe(target));
        const cmat once =
            internal::kernels::qubit::ptranspose_rho_kq(original, target, n);
        const cmat twice =
            internal::kernels::qubit::ptranspose_rho_kq(once, target, n);
        expect_same(twice, original);
    }
}

TEST(qpp_ptranspose_rho_kq, DisjointPartialTransposesCommute) {
    const idx n = 4;
    const cmat rho = make_general_matrix(n);
    const cmat t0 =
        internal::kernels::qubit::ptranspose_rho_kq(rho, {0}, n);
    const cmat t2 =
        internal::kernels::qubit::ptranspose_rho_kq(rho, {2}, n);
    const cmat first =
        internal::kernels::qubit::ptranspose_rho_kq(t0, {2}, n);
    const cmat second =
        internal::kernels::qubit::ptranspose_rho_kq(t2, {0}, n);
    const cmat both =
        internal::kernels::qubit::ptranspose_rho_kq(rho, {0, 2}, n);
    expect_same(first, second);
    expect_same(first, both);
    expect_same(first, reference_ptranspose(rho, {0, 2}, n));
}

TEST(qpp_ptranspose_rho_kq, Linearity) {
    const idx n = 3;
    const cmat A = make_density(n);
    const cmat B = make_general_matrix(n);
    const cplx alpha{0.7, 0.2};
    const cplx beta{-0.4, 0.3};
    const cmat combination = alpha * A + beta * B;
    const std::vector<idx> target{0, 2};

    const cmat actual =
        internal::kernels::qubit::ptranspose_rho_kq(combination, target, n);
    const cmat expected =
        alpha * internal::kernels::qubit::ptranspose_rho_kq(A, target, n) +
        beta * internal::kernels::qubit::ptranspose_rho_kq(B, target, n);
    expect_same(actual, expected);
    expect_same(actual, reference_ptranspose(combination, target, n));
}

TEST(qpp_ptranspose_rho_kq, PureStateConsistency) {
    const idx n = 5;
    const ket psi = make_state(n);
    const cmat rho = psi * psi.adjoint();
    for (const auto& target :
         std::vector<std::vector<idx>>{{0}, {3, 1}, {4, 0, 2}}) {
        SCOPED_TRACE(describe(target));
        const cmat from_psi =
            internal::kernels::qubit::ptranspose_psi_kq(psi, target, n);
        const cmat from_rho =
            internal::kernels::qubit::ptranspose_rho_kq(rho, target, n);
        expect_same(from_psi, from_rho);
        expect_same(from_rho, reference_ptranspose(rho, target, n));
    }
}
