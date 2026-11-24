#include <cmath>
#include <complex>
#include <numeric>
#include <vector>

#include "gtest/gtest.h"

#include "qpp/qpp.hpp"

#include "qpp/internal/kernels/qudit/apply.hpp"

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

cmat complex_gate(idx D) {
    cmat A(D, D);
    const double pi = std::acos(-1.0);
    for (idx i = 0; i < D; ++i) {
        for (idx j = 0; j < D; ++j) {
            const double theta = 2 * pi * static_cast<double>(i * j) /
                                 static_cast<double>(D);
            A(i, j) = std::exp(cplx{0, theta}) /
                      std::sqrt(static_cast<double>(D));
        }
    }
    return A;
}

cmat nonunitary_gate(idx D) {
    cmat A = complex_gate(D);
    for (idx i = 0; i < D; ++i) {
        A.row(i) *= 0.2 * static_cast<double>(i);
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

// Independent reference: embed the local gate in the full Hilbert space
// using mixed-radix basis indexes, without calling any QPP kernels.
std::vector<idx> decode(idx i, const std::vector<idx>& dims) {
    std::vector<idx> digits(dims.size());
    for (idx k = dims.size(); k > 0; --k) {
        digits[k - 1] = i % dims[k - 1];
        i /= dims[k - 1];
    }
    return digits;
}

cmat reference_operator(const cmat& A, const Case& c) {
    const idx D = dimension(c.dims);
    cmat U = cmat::Zero(D, D);
    std::vector<bool> is_target(c.dims.size(), false);
    for (idx k : c.target) {
        is_target[k] = true;
    }
    for (idx i = 0; i < D; ++i) {
        const auto row = decode(i, c.dims);
        idx m = 0;
        for (idx k : c.target) {
            m = m * c.dims[k] + row[k];
        }
        for (idx j = 0; j < D; ++j) {
            const auto col = decode(j, c.dims);
            bool same_complement = true;
            for (idx k = 0; k < c.dims.size(); ++k) {
                if (!is_target[k] && row[k] != col[k]) {
                    same_complement = false;
                    break;
                }
            }
            if (!same_complement) {
                continue;
            }
            idx n = 0;
            for (idx k : c.target) {
                n = n * c.dims[k] + col[k];
            }
            U(i, j) = A(m, n);
        }
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
    {{2, 3}, {1, 0}},
    {{3, 3}, {0, 1}},
    {{2, 3}, {}},
};
} // namespace

// =============================================================================
// apply_psi_kq_inplace()
// =============================================================================

TEST(qpp_qudit_apply_psi_kq, Reference) {
    for (idx i = 0; i < cases.size(); ++i) {
        SCOPED_TRACE(i);
        const auto& c = cases[i];
        const ket psi = make_state(dimension(c.dims));
        const cmat A = complex_gate(target_dimension(c));
        const ket expected = reference_operator(A, c) * psi;
        ket actual = psi;
        qudit::apply_psi_kq_inplace(actual, A, c.target, c.dims);
        expect_same(actual, expected);
    }
}

TEST(qpp_qudit_apply_psi_kq, BasisStates) {
    const Case c{{2, 3, 4}, {2, 0}};
    const idx D = dimension(c.dims);
    const cmat A = nonunitary_gate(target_dimension(c));
    const cmat U = reference_operator(A, c);
    for (idx i = 0; i < D; ++i) {
        SCOPED_TRACE(i);
        ket actual = ket::Zero(D);
        actual(i) = 1;
        qudit::apply_psi_kq_inplace(actual, A, c.target, c.dims);
        expect_same(actual, U.col(i));
    }
}

TEST(qpp_qudit_apply_psi_kq, UnsortedTargets) {
    const std::vector<idx> dims{2, 3, 4};
    cmat A = cmat::Zero(8, 8);
    for (idx i = 0; i < 8; ++i) {
        A((i + 1) % 8, i) = 1;
    }
    // Target digits (2, 1) encode 5; the gate maps 5 to 6,
    // corresponding to target digits (3, 0).
    ket actual = mket({1, 2, 2}, dims);
    const ket expected = mket({0, 2, 3}, dims);
    qudit::apply_psi_kq_inplace(actual, A, {2, 0}, dims);
    expect_same(actual, expected);
}

TEST(qpp_qudit_apply_psi_kq, EntanglingGate) {
    const std::vector<idx> dims{3, 3};
    cmat A = cmat::Zero(9, 9);
    for (idx i = 0; i < 3; ++i) {
        for (idx j = 0; j < 3; ++j) {
            A(3 * i + (i + j) % 3, 3 * i + j) = 1;
        }
    }
    ket actual = ket::Zero(9);
    ket expected = ket::Zero(9);
    for (idx i = 0; i < 3; ++i) {
        actual(3 * i) = 1 / std::sqrt(3.0);
        expected(4 * i) = 1 / std::sqrt(3.0);
    }
    qudit::apply_psi_kq_inplace(actual, A, {0, 1}, dims);
    expect_same(actual, expected);
}

TEST(qpp_qudit_apply_psi_kq, Identity) {
    const ket psi = make_state(24);
    const cmat A = cmat::Identity(8, 8);
    ket actual = psi;
    qudit::apply_psi_kq_inplace(actual, A, {2, 0}, {2, 3, 4});
    expect_same(actual, psi);
}

TEST(qpp_qudit_apply_psi_kq, NormPreservation) {
    for (idx i = 0; i < cases.size(); ++i) {
        SCOPED_TRACE(i);
        const auto& c = cases[i];
        const ket psi = make_state(dimension(c.dims));
        const cmat A = complex_gate(target_dimension(c));
        ket actual = psi;
        qudit::apply_psi_kq_inplace(actual, A, c.target, c.dims);
        EXPECT_NEAR(actual.norm(), psi.norm(), tolerance);
    }
}

TEST(qpp_qudit_apply_psi_kq, NonUnitaryGate) {
    for (idx i = 0; i < cases.size(); ++i) {
        SCOPED_TRACE(i);
        const auto& c = cases[i];
        const ket psi = make_state(dimension(c.dims));
        const cmat A = nonunitary_gate(target_dimension(c));
        const ket expected = reference_operator(A, c) * psi;
        ket actual = psi;
        qudit::apply_psi_kq_inplace(actual, A, c.target, c.dims);
        expect_same(actual, expected);
    }
}

TEST(qpp_qudit_apply_psi_kq, ZeroGate) {
    ket actual = make_state(24);
    const cmat A = cmat::Zero(8, 8);
    qudit::apply_psi_kq_inplace(actual, A, {2, 0}, {2, 3, 4});
    expect_same(actual, ket::Zero(24));
}

TEST(qpp_qudit_apply_psi_kq, ZeroState) {
    ket actual = ket::Zero(24);
    const cmat A = complex_gate(8);
    qudit::apply_psi_kq_inplace(actual, A, {2, 0}, {2, 3, 4});
    expect_same(actual, ket::Zero(24));
}

TEST(qpp_qudit_apply_psi_kq, TrivialDimension) {
    // The supplied kernel returns unchanged when D == 1.
    ket actual(1);
    actual(0) = cplx{0.3, 0.7};
    cmat A(1, 1);
    A(0, 0) = cplx{2, -3};
    const ket expected = actual;
    qudit::apply_psi_kq_inplace(actual, A, {}, {});
    expect_same(actual, expected);
}

// =============================================================================
// apply_rho_kq_inplace()
// =============================================================================

TEST(qpp_qudit_apply_rho_kq, Reference) {
    for (idx i = 0; i < cases.size(); ++i) {
        SCOPED_TRACE(i);
        const auto& c = cases[i];
        const idx D = dimension(c.dims);
        const ket psi = make_state(D);
        const cmat A = complex_gate(target_dimension(c));
        const cmat rho = 0.7 * psi * psi.adjoint() +
                         0.3 * cmat::Identity(D, D) / static_cast<double>(D);
        const cmat U = reference_operator(A, c);
        const cmat expected = U * rho * U.adjoint();
        cmat actual = rho;
        qudit::apply_rho_kq_inplace(actual, A, c.target, c.dims);
        expect_same(actual, expected);
    }
}

TEST(qpp_qudit_apply_rho_kq, PureStateConsistencyAndTrace) {
    for (idx i = 0; i < cases.size(); ++i) {
        SCOPED_TRACE(i);
        const auto& c = cases[i];
        ket psi = make_state(dimension(c.dims));
        cmat rho = psi * psi.adjoint();
        const cmat A = complex_gate(target_dimension(c));
        qudit::apply_psi_kq_inplace(psi, A, c.target, c.dims);
        qudit::apply_rho_kq_inplace(rho, A, c.target, c.dims);
        expect_same(rho, psi * psi.adjoint());
        EXPECT_NEAR(rho.trace().real(), 1.0, tolerance);
        EXPECT_NEAR(rho.trace().imag(), 0.0, tolerance);
        expect_same(rho, rho.adjoint());
    }
}

TEST(qpp_qudit_apply_rho_kq, MatrixUnits) {
    const Case c{{2, 3}, {1, 0}};
    const idx D = dimension(c.dims);
    const cmat A = nonunitary_gate(D);
    const cmat U = reference_operator(A, c);
    for (idx i = 0; i < D; ++i) {
        for (idx j = 0; j < D; ++j) {
            SCOPED_TRACE(i * D + j);
            cmat actual = cmat::Zero(D, D);
            actual(i, j) = 1;
            const cmat expected = U * actual * U.adjoint();
            qudit::apply_rho_kq_inplace(actual, A, c.target, c.dims);
            expect_same(actual, expected);
        }
    }
}

TEST(qpp_qudit_apply_rho_kq, Identity) {
    const ket psi = make_state(24);
    const cmat rho = psi * psi.adjoint();
    const cmat A = cmat::Identity(8, 8);
    cmat actual = rho;
    qudit::apply_rho_kq_inplace(actual, A, {2, 0}, {2, 3, 4});
    expect_same(actual, rho);
}

TEST(qpp_qudit_apply_rho_kq, MaximallyMixedState) {
    const cmat rho = cmat::Identity(24, 24) / 24.0;
    const cmat A = complex_gate(8);
    cmat actual = rho;
    qudit::apply_rho_kq_inplace(actual, A, {2, 0}, {2, 3, 4});
    expect_same(actual, rho);
}

TEST(qpp_qudit_apply_rho_kq, NonUnitaryGate) {
    for (idx i = 0; i < cases.size(); ++i) {
        SCOPED_TRACE(i);
        const auto& c = cases[i];
        const ket psi = make_state(dimension(c.dims));
        const cmat A = nonunitary_gate(target_dimension(c));
        const cmat rho = psi * psi.adjoint();
        const cmat U = reference_operator(A, c);
        const cmat expected = U * rho * U.adjoint();
        cmat actual = rho;
        qudit::apply_rho_kq_inplace(actual, A, c.target, c.dims);
        expect_same(actual, expected);
    }
}

TEST(qpp_qudit_apply_rho_kq, NonHermitianOperator) {
    const Case c{{2, 3, 4}, {2, 0}};
    const ket psi = make_state(24);
    const ket phi = mket({1, 2, 2}, c.dims);
    const cmat rho = psi * phi.adjoint();
    const cmat A = complex_gate(8);
    const cmat U = reference_operator(A, c);
    const cmat expected = U * rho * U.adjoint();
    cmat actual = rho;
    qudit::apply_rho_kq_inplace(actual, A, c.target, c.dims);
    expect_same(actual, expected);
}

TEST(qpp_qudit_apply_rho_kq, ZeroGate) {
    const ket psi = make_state(24);
    cmat actual = psi * psi.adjoint();
    const cmat A = cmat::Zero(8, 8);
    qudit::apply_rho_kq_inplace(actual, A, {2, 0}, {2, 3, 4});
    expect_same(actual, cmat::Zero(24, 24));
}

TEST(qpp_qudit_apply_rho_kq, ZeroState) {
    cmat actual = cmat::Zero(24, 24);
    const cmat A = complex_gate(8);
    qudit::apply_rho_kq_inplace(actual, A, {2, 0}, {2, 3, 4});
    expect_same(actual, cmat::Zero(24, 24));
}

TEST(qpp_qudit_apply_rho_kq, TrivialDimension) {
    // The supplied kernel returns unchanged when D == 1.
    cmat actual(1, 1);
    actual(0, 0) = 0.7;
    cmat A(1, 1);
    A(0, 0) = cplx{2, -3};
    const cmat expected = actual;
    qudit::apply_rho_kq_inplace(actual, A, {}, {});
    expect_same(actual, expected);
}
