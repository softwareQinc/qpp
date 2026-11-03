#include "gtest/gtest.h"

#include "qpp/qpp.hpp"

#include "qpp/internal/kernels/qubit/ptrace.hpp"

#include <cmath>
#include <sstream>
#include <string>
#include <vector>

using namespace qpp;

namespace {

constexpr double tolerance = 1e-11;

ket make_state(idx n, idx seed = 1) {
    ket result(idx{1} << n);
    for (idx i = 0; i < static_cast<idx>(result.size()); ++i) {
        result(i) = cplx{static_cast<double>((7 * i + 3 * seed) % 19) - 9.0,
                         static_cast<double>((5 * i + 2 * seed) % 17) - 8.0};
    }
    result.normalize();
    return result;
}

cmat make_density(idx n) {
    const ket psi = make_state(n, 1);
    const ket phi = make_state(n, 3);
    return 0.35 * (psi * psi.adjoint()) + 0.65 * (phi * phi.adjoint());
}

cmat make_general_matrix(idx n) {
    const idx D = idx{1} << n;
    cmat result(D, D);
    for (idx r = 0; r < D; ++r) {
        for (idx c = 0; c < D; ++c) {
            result(r, c) =
                cplx{0.03 * (r + 1) - 0.02 * (c + 1),
                     0.01 * (2 * r + c + 1)};
        }
    }
    return result;
}

std::string describe(const std::vector<idx>& target) {
    std::ostringstream out;
    out << "trace targets={";
    for (idx i = 0; i < target.size(); ++i) {
        if (i != 0) {
            out << ",";
        }
        out << target[i];
    }
    out << "}";
    return out.str();
}

// Independent reference: sum density-matrix elements whose traced-out
// computational-basis bits agree, keeping the remaining qubits in ascending
// big-endian subsystem order.
cmat reference_ptrace(const cmat& rho, const std::vector<idx>& target, idx n) {
    const idx D = idx{1} << n;
    std::vector<bool> traced(n, false);
    idx trace_mask = 0;
    for (idx q : target) {
        traced[q] = true;
        trace_mask |= idx{1} << (n - 1 - q);
    }

    std::vector<idx> kept;
    for (idx q = 0; q < n; ++q) {
        if (!traced[q]) {
            kept.push_back(q);
        }
    }
    const idx kept_dim = idx{1} << kept.size();
    cmat result = cmat::Zero(kept_dim, kept_dim);

    auto compress = [&](idx basis) {
        idx compressed = 0;
        for (idx q : kept) {
            compressed = (compressed << 1) |
                         ((basis >> (n - 1 - q)) & idx{1});
        }
        return compressed;
    };

    for (idx r = 0; r < D; ++r) {
        for (idx c = 0; c < D; ++c) {
            if ((r & trace_mask) == (c & trace_mask)) {
                result(compress(r), compress(c)) += rho(r, c);
            }
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

TEST(qpp_ptrace_psi_kq, EmptyTraceProducesPureDensityMatrix) {
    for (idx n : {idx{1}, idx{3}, idx{5}}) {
        const ket psi = make_state(n);
        const cmat actual =
            internal::kernels::qubit::ptrace_psi_kq(psi, {}, n);
        const cmat expected = psi * psi.adjoint();
        expect_same(actual, expected);
    }
}

TEST(qpp_ptrace_psi_kq, FullTraceProducesNormSquared) {
    for (idx n : {idx{1}, idx{2}, idx{4}}) {
        ket psi = make_state(n);
        psi *= cplx{2.0, -0.5};
        std::vector<idx> target;
        for (idx q = 0; q < n; ++q) {
            target.push_back(q);
        }
        const cmat actual =
            internal::kernels::qubit::ptrace_psi_kq(psi, target, n);
        ASSERT_EQ(actual.rows(), 1);
        ASSERT_EQ(actual.cols(), 1);
        EXPECT_NEAR(actual(0, 0).real(), psi.squaredNorm(), tolerance);
        EXPECT_NEAR(actual(0, 0).imag(), 0.0, tolerance);
    }
}

TEST(qpp_ptrace_psi_kq, BellStateEitherQubit) {
    ket psi = ket::Zero(4);
    psi(0) = 1.0 / std::sqrt(2.0);
    psi(3) = 1.0 / std::sqrt(2.0);
    const cmat expected = 0.5 * cmat::Identity(2, 2);

    for (idx q : {idx{0}, idx{1}}) {
        const cmat actual =
            internal::kernels::qubit::ptrace_psi_kq(psi, {q}, 2);
        expect_same(actual, expected);
    }
}

TEST(qpp_ptrace_psi_kq, GHZStateTwoQubitsTraced) {
    ket psi = ket::Zero(8);
    psi(0) = 1.0 / std::sqrt(2.0);
    psi(7) = 1.0 / std::sqrt(2.0);
    const cmat expected = 0.5 * cmat::Identity(2, 2);
    for (const auto& target :
         std::vector<std::vector<idx>>{{0, 1}, {1, 2}, {2, 0}, {2, 1}}) {
        SCOPED_TRACE(describe(target));
        const cmat actual =
            internal::kernels::qubit::ptrace_psi_kq(psi, target, 3);
        expect_same(actual, expected);
    }
}

TEST(qpp_ptrace_psi_kq, ProductStateTracesMiddleQubit) {
    // |0> tensor |+> tensor |1> -> |01><01| after tracing qubit 1.
    ket psi = ket::Zero(8);
    psi(1) = 1.0 / std::sqrt(2.0);
    psi(3) = 1.0 / std::sqrt(2.0);
    cmat expected = cmat::Zero(4, 4);
    expected(1, 1) = 1;
    const cmat actual =
        internal::kernels::qubit::ptrace_psi_kq(psi, {1}, 3);
    expect_same(actual, expected);
}

TEST(qpp_ptrace_psi_kq, BruteForceAllTargetConfigurations) {
    const idx n = 5;
    const ket psi = make_state(n);
    const ket before = psi;
    const cmat rho = psi * psi.adjoint();
    const std::vector<std::vector<idx>> targets{
        {}, {0}, {4}, {1, 3}, {3, 1}, {0, 2, 4},
        {4, 0, 2}, {4, 1, 3, 0}, {0, 1, 2, 3, 4}};
    for (const auto& target : targets) {
        SCOPED_TRACE(describe(target));
        const cmat actual =
            internal::kernels::qubit::ptrace_psi_kq(psi, target, n);
        const cmat expected = reference_ptrace(rho, target, n);
        expect_same(actual, expected);
        EXPECT_NEAR(actual.trace().real(), 1.0, tolerance);
        expect_same(psi, before);
    }
}

TEST(qpp_ptrace_psi_kq, TracedTargetOrderDoesNotMatter) {
    const idx n = 4;
    const ket psi = make_state(n);
    const cmat first =
        internal::kernels::qubit::ptrace_psi_kq(psi, {0, 2}, n);
    const cmat second =
        internal::kernels::qubit::ptrace_psi_kq(psi, {2, 0}, n);
    expect_same(first, second);
}

TEST(qpp_ptrace_psi_kq, ZeroStateAndZeroQubitSystem) {
    const ket zero = ket::Zero(8);
    const cmat actual =
        internal::kernels::qubit::ptrace_psi_kq(zero, {0, 2}, 3);
    expect_same(actual, cmat::Zero(2, 2));

    ket scalar(1);
    scalar(0) = cplx{1.5, -0.5};
    const cmat scalar_trace =
        internal::kernels::qubit::ptrace_psi_kq(scalar, {}, 0);
    ASSERT_EQ(scalar_trace.size(), 1);
    EXPECT_NEAR(scalar_trace(0, 0).real(), 2.5, tolerance);
    EXPECT_NEAR(scalar_trace(0, 0).imag(), 0.0, tolerance);
}

TEST(qpp_ptrace_rho_kq, EmptyTraceReturnsOriginalMatrix) {
    for (idx n : {idx{1}, idx{3}, idx{5}}) {
        const cmat rho = make_density(n);
        const cmat before = rho;
        const cmat actual =
            internal::kernels::qubit::ptrace_rho_kq(rho, {}, n);
        expect_same(actual, before);
        expect_same(rho, before);
    }
}

TEST(qpp_ptrace_rho_kq, FullTraceReturnsScalarTrace) {
    for (idx n : {idx{1}, idx{2}, idx{4}}) {
        const cmat rho = make_general_matrix(n);
        std::vector<idx> target;
        for (idx q = 0; q < n; ++q) {
            target.push_back(q);
        }
        const cmat actual =
            internal::kernels::qubit::ptrace_rho_kq(rho, target, n);
        ASSERT_EQ(actual.rows(), 1);
        ASSERT_EQ(actual.cols(), 1);
        EXPECT_LE(std::abs(actual(0, 0) - rho.trace()), tolerance);
    }
}

TEST(qpp_ptrace_rho_kq, BellStatePartialTrace) {
    ket psi = ket::Zero(4);
    psi(0) = 1.0 / std::sqrt(2.0);
    psi(3) = 1.0 / std::sqrt(2.0);
    const cmat rho = psi * psi.adjoint();
    const cmat expected = 0.5 * cmat::Identity(2, 2);
    for (idx q : {idx{0}, idx{1}}) {
        const cmat actual =
            internal::kernels::qubit::ptrace_rho_kq(rho, {q}, 2);
        expect_same(actual, expected);
    }
}

TEST(qpp_ptrace_rho_kq, MixedStateAllTargetConfigurations) {
    const idx n = 5;
    const cmat rho = make_density(n);
    const cmat before = rho;
    const std::vector<std::vector<idx>> targets{
        {}, {0}, {4}, {1, 3}, {3, 1}, {0, 2, 4},
        {4, 0, 2}, {4, 1, 3, 0}, {0, 1, 2, 3, 4}};
    for (const auto& target : targets) {
        SCOPED_TRACE(describe(target));
        const cmat actual =
            internal::kernels::qubit::ptrace_rho_kq(rho, target, n);
        const cmat expected = reference_ptrace(rho, target, n);
        expect_same(actual, expected);
        EXPECT_NEAR(actual.trace().real(), 1.0, tolerance);
        expect_same(actual, actual.adjoint());
        expect_same(rho, before);
    }
}

TEST(qpp_ptrace_rho_kq, ArbitraryNonHermitianMatrix) {
    const idx n = 4;
    const cmat original = make_general_matrix(n);
    for (const auto& target :
         std::vector<std::vector<idx>>{{0}, {3}, {1, 3}, {2, 0, 3}}) {
        SCOPED_TRACE(describe(target));
        const cmat actual =
            internal::kernels::qubit::ptrace_rho_kq(original, target, n);
        expect_same(actual, reference_ptrace(original, target, n));
    }
}

TEST(qpp_ptrace_rho_kq, Linearity) {
    const idx n = 4;
    const cmat rho = make_density(n);
    const cmat M = make_general_matrix(n);
    const cplx a{0.7, 0.3};
    const cplx b{-0.2, 0.5};
    const cmat combination = a * rho + b * M;
    const std::vector<idx> target{3, 1};

    const cmat actual =
        internal::kernels::qubit::ptrace_rho_kq(combination, target, n);
    const cmat expected =
        a * internal::kernels::qubit::ptrace_rho_kq(rho, target, n) +
        b * internal::kernels::qubit::ptrace_rho_kq(M, target, n);
    expect_same(actual, expected);
    expect_same(actual, reference_ptrace(combination, target, n));
}

TEST(qpp_ptrace_rho_kq, TracedTargetOrderDoesNotMatter) {
    const idx n = 4;
    const cmat rho = make_density(n);
    const cmat first =
        internal::kernels::qubit::ptrace_rho_kq(rho, {0, 2}, n);
    const cmat second =
        internal::kernels::qubit::ptrace_rho_kq(rho, {2, 0}, n);
    expect_same(first, second);
}

TEST(qpp_ptrace_rho_kq, PureStateConsistency) {
    const idx n = 5;
    const ket psi = make_state(n);
    const cmat rho = psi * psi.adjoint();
    for (const auto& target :
         std::vector<std::vector<idx>>{{0}, {3, 1}, {4, 0, 2}}) {
        SCOPED_TRACE(describe(target));
        const cmat from_psi =
            internal::kernels::qubit::ptrace_psi_kq(psi, target, n);
        const cmat from_rho =
            internal::kernels::qubit::ptrace_rho_kq(rho, target, n);
        expect_same(from_psi, from_rho);
    }
}

TEST(qpp_ptrace_rho_kq, ZeroQubitSystem) {
    cmat scalar(1, 1);
    scalar(0, 0) = cplx{0.7, -0.2};
    const cmat actual =
        internal::kernels::qubit::ptrace_rho_kq(scalar, {}, 0);
    expect_same(actual, scalar);
}
