#include <cmath>
#include <complex>
#include <numeric>
#include <vector>

#include "gtest/gtest.h"

#include "qpp/qpp.hpp"

#include "qpp/internal/kernels/qudit/syspermute.hpp"

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

std::vector<idx> decode(idx i, const std::vector<idx>& dims) {
    std::vector<idx> digits(dims.size());
    for (idx k = dims.size(); k > 0; --k) {
        digits[k - 1] = i % dims[k - 1];
        i /= dims[k - 1];
    }
    return digits;
}

ket make_state(idx D) {
    ket psi(D);
    for (idx i = 0; i < D; ++i) {
        const double x = static_cast<double>(i + 1);
        psi(i) = cplx{std::sin(x), std::cos(2 * x)};
    }
    return psi.normalized();
}

cmat make_operator(idx D) {
    cmat A(D, D);
    for (idx i = 0; i < D; ++i) {
        for (idx j = 0; j < D; ++j) {
            const double x = static_cast<double>(i + 1);
            const double y = static_cast<double>(j + 1);
            A(i, j) = cplx{std::sin(x + 2 * y), std::cos(3 * x - y)};
        }
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

// Independent reference: build a permutation matrix in the basis ordering.
cmat permutation_matrix(const Case& c) {
    const idx D = dimension(c.dims);
    cmat P = cmat::Zero(D, D);
    for (idx i = 0; i < D; ++i) {
        const auto digits = decode(i, c.dims);
        idx pos = 0;
        for (idx k : c.target) {
            pos = pos * c.dims[k] + digits[k];
        }
        P(pos, i) = 1;
    }
    return P;
}

std::vector<idx> permuted_dims(const Case& c) {
    std::vector<idx> result;
    for (idx k : c.target) {
        result.push_back(c.dims[k]);
    }
    return result;
}

const std::vector<Case> cases{
    {{3}, {0}}, {{3, 3}, {1, 0}}, {{4, 4}, {1, 0}},
    {{2, 3, 4}, {0, 1, 2}}, {{2, 3, 4}, {2, 0, 1}},
    {{2, 3, 4}, {1, 2, 0}}, {{2, 3, 4}, {2, 1, 0}},
    {{3, 2, 3, 2}, {3, 1, 0, 2}}, {{1}, {0}},
};
} // namespace

TEST(qpp_qudit_syspermute_psi_kq, Reference) {
    for (idx i = 0; i < cases.size(); ++i) {
        SCOPED_TRACE(i);
        const auto& c = cases[i];
        const ket psi = make_state(dimension(c.dims));
        const cmat expected = permutation_matrix(c) * psi;
        expect_same(qudit::syspermute_psi_kq(psi, c.target, c.dims), expected);
    }
}

TEST(qpp_qudit_syspermute_rho_kq, Reference) {
    for (idx i = 0; i < cases.size(); ++i) {
        SCOPED_TRACE(i);
        const auto& c = cases[i];
        const cmat A = make_operator(dimension(c.dims));
        const cmat P = permutation_matrix(c);
        const cmat expected = P * A * P.adjoint();
        expect_same(qudit::syspermute_rho_kq(A, c.target, c.dims), expected);
    }
}

TEST(qpp_qudit_syspermute_rho_kq, PureStateConsistency) {
    for (idx i = 0; i < cases.size(); ++i) {
        SCOPED_TRACE(i);
        const auto& c = cases[i];
        const ket psi = make_state(dimension(c.dims));
        const cmat A = psi * psi.adjoint();
        const cmat actual = qudit::syspermute_rho_kq(A, c.target, c.dims);
        const cmat evolved = qudit::syspermute_psi_kq(psi, c.target, c.dims);
        expect_same(actual, evolved * evolved.adjoint());
        EXPECT_NEAR(actual.trace().real(), 1.0, tolerance);
        EXPECT_NEAR(actual.trace().imag(), 0.0, tolerance);
        expect_same(actual, actual.adjoint());
    }
}

TEST(qpp_qudit_syspermute_rho_kq, MatrixUnits) {
    const Case c{{2, 3, 2}, {2, 0, 1}};
    const idx D = dimension(c.dims);
    const cmat P = permutation_matrix(c);
    for (idx i = 0; i < D; ++i) {
        for (idx j = 0; j < D; ++j) {
            SCOPED_TRACE(i * D + j);
            cmat A = cmat::Zero(D, D);
            A(i, j) = cplx{0.5, -0.75};
            const cmat expected = P * A * P.adjoint();
            expect_same(qudit::syspermute_rho_kq(A, c.target, c.dims), expected);
        }
    }
}

TEST(qpp_qudit_syspermute_psi_kq, ZeroState) {
    const Case c{{2, 3, 4}, {2, 0, 1}};
    const ket psi = ket::Zero(24);
    const cmat actual = qudit::syspermute_psi_kq(psi, c.target, c.dims);
    expect_same(actual, cmat::Zero(24, 1));
}

TEST(qpp_qudit_syspermute_rho_kq, ZeroOperator) {
    const Case c{{2, 3, 4}, {2, 0, 1}};
    const cmat A = cmat::Zero(24, 24);
    const cmat actual = qudit::syspermute_rho_kq(A, c.target, c.dims);
    expect_same(actual, cmat::Zero(24, 24));
}

TEST(qpp_qudit_syspermute_psi_kq, BasisStates) {
    const Case c{{2, 3, 4}, {2, 0, 1}};
    const cmat P = permutation_matrix(c);
    for (idx i = 0; i < 24; ++i) {
        SCOPED_TRACE(i);
        ket psi = ket::Zero(24);
        psi(i) = 1;
        expect_same(qudit::syspermute_psi_kq(psi, c.target, c.dims), P.col(i));
    }
}

TEST(qpp_qudit_syspermute_psi_kq, InverseAndNorm) {
    for (const auto& c : cases) {
        const ket psi = make_state(dimension(c.dims));
        const ket evolved = qudit::syspermute_psi_kq(psi, c.target, c.dims);
        std::vector<idx> inverse(c.target.size());
        for (idx k = 0; k < c.target.size(); ++k) {
            inverse[c.target[k]] = k;
        }
        expect_same(qudit::syspermute_psi_kq(evolved, inverse, permuted_dims(c)), psi);
        EXPECT_NEAR(evolved.norm(), psi.norm(), tolerance);
    }
}

TEST(qpp_qudit_syspermute_rho_kq, InverseAndNorm) {
    for (const auto& c : cases) {
        const cmat A = make_operator(dimension(c.dims));
        const cmat B = qudit::syspermute_rho_kq(A, c.target, c.dims);
        std::vector<idx> inverse(c.target.size());
        for (idx k = 0; k < c.target.size(); ++k) {
            inverse[c.target[k]] = k;
        }
        expect_same(qudit::syspermute_rho_kq(B, inverse, permuted_dims(c)), A);
        EXPECT_NEAR(B.norm(), A.norm(), tolerance);
        EXPECT_NEAR(std::abs(B.trace() - A.trace()), 0, tolerance);
    }
}

TEST(qpp_qudit_syspermute_rho_kq, MaximallyMixedState) {
    const cmat A = cmat::Identity(24, 24) / 24.0;
    expect_same(qudit::syspermute_rho_kq(A, {2, 0, 1}, {2, 3, 4}), A);
}
