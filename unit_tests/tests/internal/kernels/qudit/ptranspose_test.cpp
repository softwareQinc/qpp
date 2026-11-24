#include <cmath>
#include <complex>
#include <numeric>
#include <vector>

#include "gtest/gtest.h"

#include "qpp/qpp.hpp"

#include "qpp/internal/kernels/qudit/ptranspose.hpp"

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

// Independent reference: map each matrix unit to its transposed location.
cmat reference(const cmat& A, const Case& c) {
    const idx D = dimension(c.dims);
    std::vector<bool> selected(c.dims.size(), false);
    for (idx t : c.target) {
        selected[t] = true;
    }
    cmat result = cmat::Zero(D, D);
    for (idx i = 0; i < D; ++i) {
        const auto a = decode(i, c.dims);
        for (idx j = 0; j < D; ++j) {
            const auto b = decode(j, c.dims);
            idx row = 0, col = 0;
            for (idx k = 0; k < c.dims.size(); ++k) {
                row = row * c.dims[k] + (selected[k] ? b[k] : a[k]);
                col = col * c.dims[k] + (selected[k] ? a[k] : b[k]);
            }
            result(row, col) = A(i, j);
        }
    }
    return result;
}

const std::vector<Case> cases{
    {{3}, {}}, {{3}, {0}}, {{3, 3}, {0}}, {{3, 3}, {1}},
    {{4, 4}, {1}}, {{2, 3, 4}, {1}}, {{2, 3, 4}, {2, 0}},
    {{3, 2, 3, 2}, {3, 1}}, {{2, 3, 4}, {}},
    {{2, 3, 4}, {2, 0, 1}}, {{1}, {0}}, {{1}, {}}, {{}, {}},
};
} // namespace

TEST(qpp_qudit_ptranspose_psi_kq, Reference) {
    for (idx i = 0; i < cases.size(); ++i) {
        SCOPED_TRACE(i);
        const auto& c = cases[i];
        const ket psi = make_state(dimension(c.dims));
        const cmat expected = reference(psi * psi.adjoint(), c);
        expect_same(qudit::ptranspose_psi_kq(psi, c.target, c.dims), expected);
    }
}

TEST(qpp_qudit_ptranspose_rho_kq, Reference) {
    for (idx i = 0; i < cases.size(); ++i) {
        SCOPED_TRACE(i);
        const auto& c = cases[i];
        const cmat A = make_operator(dimension(c.dims));
        const cmat expected = reference(A, c);
        expect_same(qudit::ptranspose_rho_kq(A, c.target, c.dims), expected);
    }
}

TEST(qpp_qudit_ptranspose_rho_kq, PureStateConsistency) {
    for (idx i = 0; i < cases.size(); ++i) {
        SCOPED_TRACE(i);
        const auto& c = cases[i];
        const ket psi = make_state(dimension(c.dims));
        const cmat A = psi * psi.adjoint();
        const cmat actual = qudit::ptranspose_rho_kq(A, c.target, c.dims);
        const cmat evolved = qudit::ptranspose_psi_kq(psi, c.target, c.dims);
        expect_same(actual, evolved);
        EXPECT_NEAR(actual.trace().real(), 1.0, tolerance);
        EXPECT_NEAR(actual.trace().imag(), 0.0, tolerance);
        expect_same(actual, actual.adjoint());
    }
}

TEST(qpp_qudit_ptranspose_rho_kq, MatrixUnits) {
    const Case c{{2, 3, 2}, {2, 0}};
    const idx D = dimension(c.dims);
    for (idx i = 0; i < D; ++i) {
        for (idx j = 0; j < D; ++j) {
            SCOPED_TRACE(i * D + j);
            cmat A = cmat::Zero(D, D);
            A(i, j) = cplx{0.5, -0.75};
            const cmat expected = reference(A, c);
            expect_same(qudit::ptranspose_rho_kq(A, c.target, c.dims), expected);
        }
    }
}

TEST(qpp_qudit_ptranspose_psi_kq, ZeroState) {
    const Case c{{2, 3, 4}, {2, 0}};
    const ket psi = ket::Zero(24);
    const cmat actual = qudit::ptranspose_psi_kq(psi, c.target, c.dims);
    expect_same(actual, cmat::Zero(24, 24));
}

TEST(qpp_qudit_ptranspose_rho_kq, ZeroOperator) {
    const Case c{{2, 3, 4}, {2, 0}};
    const cmat A = cmat::Zero(24, 24);
    const cmat actual = qudit::ptranspose_rho_kq(A, c.target, c.dims);
    expect_same(actual, cmat::Zero(24, 24));
}

TEST(qpp_qudit_ptranspose_rho_kq, InvolutionAndNorm) {
    for (const auto& c : cases) {
        const cmat A = make_operator(dimension(c.dims));
        const cmat B = qudit::ptranspose_rho_kq(A, c.target, c.dims);
        expect_same(qudit::ptranspose_rho_kq(B, c.target, c.dims), A);
        EXPECT_NEAR(B.norm(), A.norm(), tolerance);
        EXPECT_NEAR(std::abs(B.trace() - A.trace()), 0, tolerance);
    }
}

TEST(qpp_qudit_ptranspose_psi_kq, MaximallyEntangledQutrits) {
    ket psi = ket::Zero(9);
    psi(0) = psi(4) = psi(8) = 1 / std::sqrt(3.0);
    cmat expected = cmat::Zero(9, 9);
    for (idx i = 0; i < 3; ++i) {
        for (idx j = 0; j < 3; ++j) {
            expected(3 * i + j, 3 * j + i) = 1 / 3.0;
        }
    }
    expect_same(qudit::ptranspose_psi_kq(psi, {0}, {3, 3}), expected);
    expect_same(qudit::ptranspose_psi_kq(psi, {1}, {3, 3}), expected);
}

TEST(qpp_qudit_ptranspose_rho_kq, FullTransposeNotAdjoint) {
    const cmat A = make_operator(24);
    expect_same(qudit::ptranspose_rho_kq(A, {2, 0, 1}, {2, 3, 4}),
                A.transpose());
}

TEST(qpp_qudit_ptranspose_rho_kq, DisjointComposition) {
    const cmat A = make_operator(24);
    const cmat B = qudit::ptranspose_rho_kq(A, {0}, {2, 3, 4});
    expect_same(qudit::ptranspose_rho_kq(B, {2}, {2, 3, 4}),
                qudit::ptranspose_rho_kq(A, {2, 0}, {2, 3, 4}));
}
