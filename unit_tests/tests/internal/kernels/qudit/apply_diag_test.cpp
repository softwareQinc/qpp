#include <cmath>
#include <complex>
#include <numeric>
#include <vector>

#include "gtest/gtest.h"

#include "qpp/qpp.hpp"

#include "qpp/internal/kernels/qudit/apply_diag.hpp"

using namespace qpp;

namespace {
constexpr double tolerance = 1e-12;
namespace qudit = qpp::internal::kernels::qudit;

struct Case {
    std::vector<idx> dims;
    std::vector<idx> target;
};

idx dimension(const std::vector<idx>& dims) {
    return std::accumulate(dims.begin(), dims.end(), idx{1},
                           [](idx a, idx b) { return a * b; });
}

idx target_dimension(const Case& c) {
    idx D = 1;
    for (idx k : c.target) {
        D *= c.dims[k];
    }
    return D;
}

ket make_state(idx D) {
    ket psi(D);
    for (idx i = 0; i < D; ++i) {
        const auto x = static_cast<double>(i + 1);
        psi(i) = cplx{std::sin(x), std::cos(2 * x)};
    }
    return psi.normalized();
}

ket complex_gate(idx D) {
    ket A(D);
    for (idx i = 0; i < D; ++i) {
        A(i) = std::exp(cplx{0, 0.37 * static_cast<double>(i + 1)});
    }
    return A;
}

ket nonunitary_gate(idx D) {
    ket A = complex_gate(D);
    for (idx i = 0; i < D; ++i) {
        A(i) *= 0.2 * static_cast<double>(i);
    }
    return A;
}

template <typename Derived1, typename Derived2>
void expect_same(const Eigen::MatrixBase<Derived1>& actual,
                 const Eigen::MatrixBase<Derived2>& expected) {
    ASSERT_EQ(actual.rows(), expected.rows());
    ASSERT_EQ(actual.cols(), expected.cols());
    EXPECT_NEAR((actual - expected).norm(), 0, tolerance);
}

// Independent reference: decode each global basis index and gather
// the target digits in the specified order, without using QPP kernels.
cmat reference_operator(const ket& A, const Case& c) {
    const idx D = dimension(c.dims);
    cmat U = cmat::Zero(D, D);
    for (idx i = 0; i < D; ++i) {
        idx j = i;
        std::vector<idx> digits(c.dims.size());
        for (idx k = c.dims.size(); k > 0; --k) {
            digits[k - 1] = j % c.dims[k - 1];
            j /= c.dims[k - 1];
        }
        idx m = 0;
        for (idx k : c.target) {
            m = m * c.dims[k] + digits[k];
        }
        U(i, i) = A(m);
    }
    return U;
}

const std::vector<Case> cases{
    {{3}, {0}},
    {{3, 3}, {0}},
    {{3, 3}, {1}},
    {{3, 3, 3}, {0, 2}},
    {{3, 3, 3}, {2, 0}},
    {{4, 4}, {1}},
    {{2, 3, 4}, {1}},
    {{2, 3, 4}, {2, 0}},
    {{2, 3, 4}, {0, 2}},
    {{2, 3, 4}, {2, 0, 1}},
    {{3, 3, 3}, {0, 1, 2}},
    {{2, 3}, {}},
    {{1}, {0}},
    {{}, {}},
};
} // namespace

// =============================================================================
// apply_psi_kq_diag_inplace()
// =============================================================================

TEST(qpp_qudit_apply_psi_kq_diag, Reference) {
    for (idx i = 0; i < cases.size(); ++i) {
        SCOPED_TRACE(i);
        const auto& c = cases[i];
        const ket psi = make_state(dimension(c.dims));
        const ket A = complex_gate(target_dimension(c));
        const ket expected = reference_operator(A, c) * psi;
        ket actual = psi;
        qudit::apply_psi_kq_diag_inplace(actual, A, c.target, c.dims);
        expect_same(actual, expected);
    }
}

TEST(qpp_qudit_apply_psi_kq_diag, BasisStates) {
    const Case c{{2, 3, 4}, {2, 0}};
    const idx D = dimension(c.dims);
    const ket A = nonunitary_gate(target_dimension(c));
    const cmat U = reference_operator(A, c);
    for (idx i = 0; i < D; ++i) {
        SCOPED_TRACE(i);
        ket actual = ket::Zero(D);
        actual(i) = 1;
        qudit::apply_psi_kq_diag_inplace(actual, A, c.target, c.dims);
        expect_same(actual, U.col(i));
    }
}

TEST(qpp_qudit_apply_psi_kq_diag, UnsortedTargets) {
    const std::vector<idx> dims{2, 3, 4};
    ket A(8);
    for (idx i = 0; i < 8; ++i) {
        A(i) = cplx{static_cast<double>(i + 1), -0.5};
    }
    const ket psi = mket({1, 2, 2}, dims);
    // Target digits (2, 1) give index 2 * 2 + 1 = 5.
    const ket expected = A(5) * psi;
    ket actual = psi;
    qudit::apply_psi_kq_diag_inplace(actual, A, {2, 0}, dims);
    expect_same(actual, expected);
}

TEST(qpp_qudit_apply_psi_kq_diag, Identity) {
    const ket psi = make_state(24);
    const ket A = ket::Ones(8);
    ket actual = psi;
    qudit::apply_psi_kq_diag_inplace(actual, A, {2, 0}, {2, 3, 4});
    expect_same(actual, psi);
}

TEST(qpp_qudit_apply_psi_kq_diag, NormPreservation) {
    for (idx i = 0; i < cases.size(); ++i) {
        SCOPED_TRACE(i);
        const auto& c = cases[i];
        const ket psi = make_state(dimension(c.dims));
        const ket A = complex_gate(target_dimension(c));
        ket actual = psi;
        qudit::apply_psi_kq_diag_inplace(actual, A, c.target, c.dims);
        EXPECT_NEAR(actual.norm(), psi.norm(), tolerance);
    }
}

TEST(qpp_qudit_apply_psi_kq_diag, NonUnitaryGate) {
    for (idx i = 0; i < cases.size(); ++i) {
        SCOPED_TRACE(i);
        const auto& c = cases[i];
        const ket psi = make_state(dimension(c.dims));
        const ket A = nonunitary_gate(target_dimension(c));
        const ket expected = reference_operator(A, c) * psi;
        ket actual = psi;
        qudit::apply_psi_kq_diag_inplace(actual, A, c.target, c.dims);
        expect_same(actual, expected);
    }
}

TEST(qpp_qudit_apply_psi_kq_diag, ZeroGate) {
    ket actual = make_state(24);
    const ket A = ket::Zero(8);
    qudit::apply_psi_kq_diag_inplace(actual, A, {2, 0}, {2, 3, 4});
    expect_same(actual, ket::Zero(24));
}

TEST(qpp_qudit_apply_psi_kq_diag, ZeroState) {
    ket actual = ket::Zero(24);
    const ket A = complex_gate(8);
    qudit::apply_psi_kq_diag_inplace(actual, A, {2, 0}, {2, 3, 4});
    expect_same(actual, ket::Zero(24));
}

TEST(qpp_qudit_apply_psi_kq_diag, TrivialDimension) {
    ket actual(1), A(1);
    actual(0) = cplx{0.3, 0.7};
    A(0) = cplx{2, -3};
    const ket expected = A(0) * actual;
    qudit::apply_psi_kq_diag_inplace(actual, A, {}, {});
    expect_same(actual, expected);
}

// =============================================================================
// apply_rho_kq_diag_inplace()
// =============================================================================

TEST(qpp_qudit_apply_rho_kq_diag, Reference) {
    for (idx i = 0; i < cases.size(); ++i) {
        SCOPED_TRACE(i);
        const auto& c = cases[i];
        const idx D = dimension(c.dims);
        const ket psi = make_state(D);
        const ket A = complex_gate(target_dimension(c));
        const cmat rho = 0.7 * psi * psi.adjoint() +
                         0.3 * cmat::Identity(D, D) / static_cast<double>(D);
        const cmat U = reference_operator(A, c);
        const cmat expected = U * rho * U.adjoint();
        cmat actual = rho;
        qudit::apply_rho_kq_diag_inplace(actual, A, c.target, c.dims);
        expect_same(actual, expected);
    }
}

TEST(qpp_qudit_apply_rho_kq_diag, PureStateConsistencyAndTrace) {
    for (idx i = 0; i < cases.size(); ++i) {
        SCOPED_TRACE(i);
        const auto& c = cases[i];
        ket psi = make_state(dimension(c.dims));
        cmat rho = psi * psi.adjoint();
        const ket A = complex_gate(target_dimension(c));
        qudit::apply_psi_kq_diag_inplace(psi, A, c.target, c.dims);
        qudit::apply_rho_kq_diag_inplace(rho, A, c.target, c.dims);
        expect_same(rho, psi * psi.adjoint());
        EXPECT_NEAR(rho.trace().real(), 1.0, tolerance);
        EXPECT_NEAR(rho.trace().imag(), 0.0, tolerance);
        expect_same(rho, rho.adjoint());
    }
}

TEST(qpp_qudit_apply_rho_kq_diag, MatrixUnits) {
    const Case c{{2, 3}, {1, 0}};
    const idx D = dimension(c.dims);
    const ket A = nonunitary_gate(D);
    const cmat U = reference_operator(A, c);
    for (idx i = 0; i < D; ++i) {
        for (idx j = 0; j < D; ++j) {
            SCOPED_TRACE(i * D + j);
            cmat actual = cmat::Zero(D, D);
            actual(i, j) = 1;
            const cmat expected = U * actual * U.adjoint();
            qudit::apply_rho_kq_diag_inplace(actual, A, c.target, c.dims);
            expect_same(actual, expected);
        }
    }
}

TEST(qpp_qudit_apply_rho_kq_diag, Identity) {
    const ket psi = make_state(24);
    const cmat rho = psi * psi.adjoint();
    const ket A = ket::Ones(8);
    cmat actual = rho;
    qudit::apply_rho_kq_diag_inplace(actual, A, {2, 0}, {2, 3, 4});
    expect_same(actual, rho);
}

TEST(qpp_qudit_apply_rho_kq_diag, MaximallyMixedState) {
    const cmat rho = cmat::Identity(24, 24) / 24.0;
    const ket A = complex_gate(8);
    cmat actual = rho;
    qudit::apply_rho_kq_diag_inplace(actual, A, {2, 0}, {2, 3, 4});
    expect_same(actual, rho);
}

TEST(qpp_qudit_apply_rho_kq_diag, NonUnitaryGate) {
    for (idx i = 0; i < cases.size(); ++i) {
        SCOPED_TRACE(i);
        const auto& c = cases[i];
        const ket psi = make_state(dimension(c.dims));
        const ket A = nonunitary_gate(target_dimension(c));
        const cmat rho = psi * psi.adjoint();
        const cmat U = reference_operator(A, c);
        const cmat expected = U * rho * U.adjoint();
        cmat actual = rho;
        qudit::apply_rho_kq_diag_inplace(actual, A, c.target, c.dims);
        expect_same(actual, expected);
    }
}

TEST(qpp_qudit_apply_rho_kq_diag, NonHermitianOperator) {
    const Case c{{2, 3, 4}, {2, 0}};
    const ket psi = make_state(24);
    const ket phi = mket({1, 2, 2}, c.dims);
    const cmat rho = psi * phi.adjoint();
    const ket A = complex_gate(8);
    const cmat U = reference_operator(A, c);
    const cmat expected = U * rho * U.adjoint();
    cmat actual = rho;
    qudit::apply_rho_kq_diag_inplace(actual, A, c.target, c.dims);
    expect_same(actual, expected);
}

TEST(qpp_qudit_apply_rho_kq_diag, ZeroGate) {
    const ket psi = make_state(24);
    cmat actual = psi * psi.adjoint();
    const ket A = ket::Zero(8);
    qudit::apply_rho_kq_diag_inplace(actual, A, {2, 0}, {2, 3, 4});
    expect_same(actual, cmat::Zero(24, 24));
}

TEST(qpp_qudit_apply_rho_kq_diag, ZeroState) {
    cmat actual = cmat::Zero(24, 24);
    const ket A = complex_gate(8);
    qudit::apply_rho_kq_diag_inplace(actual, A, {2, 0}, {2, 3, 4});
    expect_same(actual, cmat::Zero(24, 24));
}

TEST(qpp_qudit_apply_rho_kq_diag, TrivialDimension) {
    cmat actual(1, 1);
    actual(0, 0) = 0.7;
    ket A(1);
    A(0) = cplx{2, -3};
    const cmat expected = std::norm(A(0)) * actual;
    qudit::apply_rho_kq_diag_inplace(actual, A, {}, {});
    expect_same(actual, expected);
}
