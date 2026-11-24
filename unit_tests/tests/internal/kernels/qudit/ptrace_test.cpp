#include <cmath>
#include <complex>
#include <numeric>
#include <vector>

#include "gtest/gtest.h"

#include "qpp/qpp.hpp"

#include "qpp/internal/kernels/qudit/ptrace.hpp"

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

// Independent reference: sum matrix entries whose traced digits agree.
cmat reference(const cmat& A, const Case& c) {
    std::vector<bool> traced(c.dims.size(), false);
    for (idx t : c.target) {
        traced[t] = true;
    }
    idx R = 1;
    for (idx k = 0; k < c.dims.size(); ++k) {
        if (!traced[k]) {
            R *= c.dims[k];
        }
    }
    cmat result = cmat::Zero(R, R);
    const idx D = dimension(c.dims);
    for (idx i = 0; i < D; ++i) {
        const auto a = decode(i, c.dims);
        for (idx j = 0; j < D; ++j) {
            const auto b = decode(j, c.dims);
            idx row = 0, col = 0;
            bool match = true;
            for (idx k = 0; k < c.dims.size(); ++k) {
                if (traced[k]) {
                    if (a[k] != b[k]) {
                        match = false;
                        break;
                    }
                } else {
                    row = row * c.dims[k] + a[k];
                    col = col * c.dims[k] + b[k];
                }
            }
            if (match) {
                result(row, col) += A(i, j);
            }
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

TEST(qpp_qudit_ptrace_psi_kq, Reference) {
    for (idx i = 0; i < cases.size(); ++i) {
        SCOPED_TRACE(i);
        const auto& c = cases[i];
        const ket psi = make_state(dimension(c.dims));
        const cmat expected = reference(psi * psi.adjoint(), c);
        expect_same(qudit::ptrace_psi_kq(psi, c.target, c.dims), expected);
    }
}

TEST(qpp_qudit_ptrace_rho_kq, Reference) {
    for (idx i = 0; i < cases.size(); ++i) {
        SCOPED_TRACE(i);
        const auto& c = cases[i];
        const cmat A = make_operator(dimension(c.dims));
        const cmat expected = reference(A, c);
        expect_same(qudit::ptrace_rho_kq(A, c.target, c.dims), expected);
    }
}

TEST(qpp_qudit_ptrace_rho_kq, PureStateConsistency) {
    for (idx i = 0; i < cases.size(); ++i) {
        SCOPED_TRACE(i);
        const auto& c = cases[i];
        const ket psi = make_state(dimension(c.dims));
        const cmat A = psi * psi.adjoint();
        const cmat actual = qudit::ptrace_rho_kq(A, c.target, c.dims);
        const cmat evolved = qudit::ptrace_psi_kq(psi, c.target, c.dims);
        expect_same(actual, evolved);
        EXPECT_NEAR(actual.trace().real(), 1.0, tolerance);
        EXPECT_NEAR(actual.trace().imag(), 0.0, tolerance);
        expect_same(actual, actual.adjoint());
    }
}

TEST(qpp_qudit_ptrace_rho_kq, MatrixUnits) {
    const Case c{{2, 3, 2}, {2, 0}};
    const idx D = dimension(c.dims);
    for (idx i = 0; i < D; ++i) {
        for (idx j = 0; j < D; ++j) {
            SCOPED_TRACE(i * D + j);
            cmat A = cmat::Zero(D, D);
            A(i, j) = cplx{0.5, -0.75};
            const cmat expected = reference(A, c);
            expect_same(qudit::ptrace_rho_kq(A, c.target, c.dims), expected);
        }
    }
}

TEST(qpp_qudit_ptrace_psi_kq, ZeroState) {
    const Case c{{2, 3, 4}, {2, 0}};
    const ket psi = ket::Zero(24);
    const cmat actual = qudit::ptrace_psi_kq(psi, c.target, c.dims);
    expect_same(actual, cmat::Zero(3, 3));
}

TEST(qpp_qudit_ptrace_rho_kq, ZeroOperator) {
    const Case c{{2, 3, 4}, {2, 0}};
    const cmat A = cmat::Zero(24, 24);
    const cmat actual = qudit::ptrace_rho_kq(A, c.target, c.dims);
    expect_same(actual, cmat::Zero(3, 3));
}

TEST(qpp_qudit_ptrace_psi_kq, MaximallyEntangledQutrits) {
    ket psi = ket::Zero(9);
    psi(0) = psi(4) = psi(8) = 1 / std::sqrt(3.0);
    const cmat expected = cmat::Identity(3, 3) / 3.0;
    expect_same(qudit::ptrace_psi_kq(psi, {0}, {3, 3}), expected);
    expect_same(qudit::ptrace_psi_kq(psi, {1}, {3, 3}), expected);
}

TEST(qpp_qudit_ptrace_psi_kq, UnnormalizedFullTrace) {
    const ket psi = cplx{2, 3} * make_state(24);
    cmat expected(1, 1);
    expected(0, 0) = 13;
    expect_same(qudit::ptrace_psi_kq(psi, {2, 0, 1}, {2, 3, 4}), expected);
}

TEST(qpp_qudit_ptrace_rho_kq, MaximallyMixedState) {
    const cmat A = cmat::Identity(24, 24) / 24.0;
    const cmat expected = cmat::Identity(3, 3) / 3.0;
    expect_same(qudit::ptrace_rho_kq(A, {2, 0}, {2, 3, 4}), expected);
}

TEST(qpp_qudit_ptrace_rho_kq, TargetOrder) {
    const cmat A = make_operator(24);
    expect_same(qudit::ptrace_rho_kq(A, {0, 2}, {2, 3, 4}),
                qudit::ptrace_rho_kq(A, {2, 0}, {2, 3, 4}));
}

TEST(qpp_qudit_ptrace_rho_kq, SequentialTrace) {
    const cmat A = make_operator(24);
    const cmat first = qudit::ptrace_rho_kq(A, {0}, {2, 3, 4});
    expect_same(qudit::ptrace_rho_kq(first, {1}, {3, 4}),
                qudit::ptrace_rho_kq(A, {0, 2}, {2, 3, 4}));
}
