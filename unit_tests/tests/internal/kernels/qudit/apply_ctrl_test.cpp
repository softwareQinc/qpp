#include <cmath>
#include <complex>
#include <numeric>
#include <vector>

#include "gtest/gtest.h"

#include "qpp/qpp.hpp"

#include "qpp/internal/kernels/qudit/apply_ctrl.hpp"

using namespace qpp;

namespace {
constexpr double tolerance = 1e-12;

namespace qudit = qpp::internal::kernels::qudit;

struct Case {
    std::vector<idx> dims;
    std::vector<idx> ctrl;
    std::vector<idx> target;
    std::vector<idx> shift;
};

idx dimension(const std::vector<idx>& dims) {
    return std::accumulate(dims.begin(), dims.end(), idx{1},
                           [](idx a, idx b) { return a * b; });
}

std::vector<idx> decode(idx i, const std::vector<idx>& dims) {
    std::vector<idx> result(dims.size());

    for (idx k = dims.size(); k > 0; --k) {
        result[k - 1] = i % dims[k - 1];
        i /= dims[k - 1];
    }

    return result;
}

idx target_dimension(const Case& test_case) {
    idx D = 1;
    for (idx k : test_case.target) {
        D *= test_case.dims[k];
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

cmat complex_gate(idx d) {
    cmat A(d, d);

    const double pi = std::acos(-1.0);

    for (idx i = 0; i < d; ++i) {
        for (idx j = 0; j < d; ++j) {
            const double theta =
                2 * pi * static_cast<double>(i * j) /
                static_cast<double>(d);

            A(i, j) = std::exp(cplx{0, theta}) /
                      std::sqrt(static_cast<double>(d));
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

// Independent reference implementation.
//
// The kernel accepts inverted control shifts, so the effective
// control value is (digit + d - shift) % d.
cmat reference_operator(const cmat& A, const Case& test_case) {
    const auto& dims = test_case.dims;
    const auto& ctrl = test_case.ctrl;
    const auto& target = test_case.target;
    const auto& shift = test_case.shift;

    const idx D = dimension(dims);
    const idx d = dims[ctrl[0]];
    const idx d_target = static_cast<idx>(A.rows());

    std::vector<cmat> Ak(d);
    Ak[0] = cmat::Identity(d_target, d_target);

    for (idx k = 1; k < d; ++k) {
        Ak[k] = Ak[k - 1] * A;
    }

    std::vector<bool> is_target(dims.size(), false);

    for (idx elem : target) {
        is_target[elem] = true;
    }

    cmat U = cmat::Zero(D, D);

    for (idx j = 0; j < D; ++j) {
        const auto input = decode(j, dims);

        const idx r =
            (input[ctrl[0]] + d - shift[0]) % d;

        bool active = true;

        for (idx k = 1; k < ctrl.size(); ++k) {
            if ((input[ctrl[k]] + d - shift[k]) % d != r) {
                active = false;
                break;
            }
        }

        if (!active || r == 0) {
            U(j, j) = 1;
            continue;
        }

        for (idx i = 0; i < D; ++i) {
            const auto output = decode(i, dims);

            bool same_spectators = true;
            for (idx k = 0; k < dims.size(); ++k) {
                if (!is_target[k] && output[k] != input[k]) {
                    same_spectators = false;
                    break;
                }
            }

            if (!same_spectators) {
                continue;
            }

            idx row = 0;
            idx col = 0;
            for (idx k : target) {
                row = row * dims[k] + output[k];
                col = col * dims[k] + input[k];
            }
            U(i, j) = Ak[r](row, col);
        }
    }

    return U;
}

const std::vector<Case> cases{
    // Single control, single target
    {{3, 3}, {0}, {1}, {0}},

    // Single control, multiple targets
    {{3, 3, 3}, {1}, {2, 0}, {0}},

    // Multiple controls, single target
    {{3, 3, 3}, {2, 0}, {1}, {0, 0}},

    // Multiple controls, multiple targets
    {{3, 3, 3, 3}, {2, 0}, {3, 1}, {0, 0}},

    // Shifted single control
    {{3, 3}, {0}, {1}, {2}},

    // Shifted multiple controls
    {{3, 3, 3, 3}, {2, 0}, {3, 1}, {1, 2}},

    // Ququarts
    {{4, 4, 4}, {1}, {2, 0}, {3}},

    // Heterogeneous dimensions
    {{3, 2, 3, 2}, {2, 0}, {3, 1}, {1, 2}},

    // Control dimension smaller than target dimension
    {{2, 3, 3}, {0}, {2, 1}, {1}},

    // Targets with different dimensions, in non-canonical order
    {{3, 2, 4, 3}, {3, 0}, {2, 1}, {1, 2}},

    // Untouched spectator subsystem
    {{3, 2, 3}, {0}, {2}, {1}},
};
} // namespace

// =============================================================================
// apply_ctrl_psi_kq_inplace()
// =============================================================================

TEST(qpp_qudit_apply_ctrl_psi_kq, Reference) {
    for (idx i = 0; i < cases.size(); ++i) {
        SCOPED_TRACE(i);

        const auto& test_case = cases[i];

        const idx D = dimension(test_case.dims);
        const ket psi = make_state(D);
        const cmat A =
            complex_gate(target_dimension(test_case));

        const ket expected =
            reference_operator(A, test_case) * psi;

        ket actual = psi;

        qudit::apply_ctrl_psi_kq_inplace(
            actual, A, test_case.ctrl, test_case.target,
            test_case.dims, test_case.shift);

        expect_same(actual, expected);
    }
}

TEST(qpp_qudit_apply_ctrl_psi_kq, BasisStates) {
    const Case test_case{{3, 3, 3, 3},
                         {2, 0},
                         {3, 1},
                         {0, 0}};

    const idx D = dimension(test_case.dims);
    const cmat A = kron(gt.Xd(3), gt.Xd(3));
    const cmat U = reference_operator(A, test_case);

    for (idx i = 0; i < D; ++i) {
        SCOPED_TRACE(i);

        ket actual = ket::Zero(D);
        actual(i) = 1;

        qudit::apply_ctrl_psi_kq_inplace(
            actual, A, test_case.ctrl, test_case.target,
            test_case.dims, test_case.shift);

        expect_same(actual, U.col(i));
    }
}

TEST(qpp_qudit_apply_ctrl_psi_kq, InactiveControls) {
    const std::vector<idx> dims{3, 3, 3, 3};
    const cmat A = kron(gt.Xd(3), gt.Xd(3));

    // The control qudits have different values.
    const ket psi = mket({1, 0, 2, 0}, dims);
    ket actual = psi;

    qudit::apply_ctrl_psi_kq_inplace(
        actual, A, {0, 2}, {1, 3}, dims, {0, 0});

    expect_same(actual, psi);
}

TEST(qpp_qudit_apply_ctrl_psi_kq, ZeroControls) {
    const std::vector<idx> dims{3, 3, 3, 3};
    const cmat A = kron(gt.Xd(3), gt.Xd(3));

    const ket psi = mket({0, 1, 0, 2}, dims);
    ket actual = psi;

    qudit::apply_ctrl_psi_kq_inplace(
        actual, A, {0, 2}, {1, 3}, dims, {0, 0});

    expect_same(actual, psi);
}

TEST(qpp_qudit_apply_ctrl_psi_kq, NonzeroControls) {
    const std::vector<idx> dims{3, 3, 3, 3};
    const cmat A = kron(gt.Xd(3), gt.Xd(3));

    // Both controls are 2, so X^2 is applied to each target.
    ket actual = mket({2, 0, 2, 1}, dims);
    const ket expected = mket({2, 2, 2, 0}, dims);

    qudit::apply_ctrl_psi_kq_inplace(
        actual, A, {0, 2}, {1, 3}, dims, {0, 0});

    expect_same(actual, expected);
}

TEST(qpp_qudit_apply_ctrl_psi_kq, ShiftedControls) {
    const std::vector<idx> dims{3, 3, 3, 3};
    const cmat A = kron(gt.Xd(3), gt.Xd(3));

    // Inverted shifts:
    // (2 - 1) % 3 = (0 - 2) % 3 = 1 (modulo 3).
    // Therefore, X is applied to both targets.
    ket actual = mket({2, 0, 0, 1}, dims);
    const ket expected = mket({2, 1, 0, 2}, dims);

    qudit::apply_ctrl_psi_kq_inplace(
        actual, A, {0, 2}, {1, 3}, dims, {1, 2});

    expect_same(actual, expected);
}

TEST(qpp_qudit_apply_ctrl_psi_kq, Identity) {
    const ket psi = make_state(81);
    const cmat A = cmat::Identity(9, 9);

    ket actual = psi;

    qudit::apply_ctrl_psi_kq_inplace(
        actual, A, {2, 0}, {3, 1}, {3, 3, 3, 3}, {1, 2});

    expect_same(actual, psi);
}

TEST(qpp_qudit_apply_ctrl_psi_kq, NormPreservation) {
    for (idx i = 0; i < cases.size(); ++i) {
        SCOPED_TRACE(i);

        const auto& test_case = cases[i];

        const ket psi = make_state(dimension(test_case.dims));
        const cmat A =
            complex_gate(target_dimension(test_case));

        ket actual = psi;

        qudit::apply_ctrl_psi_kq_inplace(
            actual, A, test_case.ctrl, test_case.target,
            test_case.dims, test_case.shift);

        EXPECT_NEAR(actual.norm(), psi.norm(), tolerance);
    }
}

TEST(qpp_qudit_apply_ctrl_psi_kq, NonUnitaryGate) {
    const Case test_case{{3, 3, 3},
                         {1},
                         {2, 0},
                         {1}};

    const ket psi = make_state(27);
    const cmat A = 0.7 * complex_gate(9);

    const ket expected =
        reference_operator(A, test_case) * psi;

    ket actual = psi;

    qudit::apply_ctrl_psi_kq_inplace(
        actual, A, test_case.ctrl, test_case.target,
        test_case.dims, test_case.shift);

    expect_same(actual, expected);
}

TEST(qpp_qudit_apply_ctrl_psi_kq, ZeroState) {
    ket actual = ket::Zero(27);
    const cmat A = complex_gate(9);

    qudit::apply_ctrl_psi_kq_inplace(
        actual, A, {0}, {1, 2}, {3, 3, 3}, {1});

    expect_same(actual, ket::Zero(27));
}


// =============================================================================
// apply_ctrl_rho_kq_inplace()
// =============================================================================

TEST(qpp_qudit_apply_ctrl_rho_kq, Reference) {
    for (idx i = 0; i < cases.size(); ++i) {
        SCOPED_TRACE(i);

        const auto& test_case = cases[i];

        const idx D = dimension(test_case.dims);
        const ket psi = make_state(D);
        const cmat A =
            complex_gate(target_dimension(test_case));

        // Mixed state
        const cmat rho =
            0.7 * psi * psi.adjoint() +
            0.3 * cmat::Identity(D, D) / static_cast<double>(D);

        const cmat U = reference_operator(A, test_case);
        const cmat expected = U * rho * U.adjoint();

        cmat actual = rho;

        qudit::apply_ctrl_rho_kq_inplace(
            actual, A, test_case.ctrl, test_case.target,
            test_case.dims, test_case.shift);

        expect_same(actual, expected);
    }
}

TEST(qpp_qudit_apply_ctrl_rho_kq,
     PureStateConsistencyAndTrace) {
    for (idx i = 0; i < cases.size(); ++i) {
        SCOPED_TRACE(i);

        const auto& test_case = cases[i];

        const idx D = dimension(test_case.dims);
        const ket psi = make_state(D);
        const cmat rho = psi * psi.adjoint();
        const cmat A =
            complex_gate(target_dimension(test_case));

        ket evolved_psi = psi;
        cmat evolved_rho = rho;

        qudit::apply_ctrl_psi_kq_inplace(
            evolved_psi, A, test_case.ctrl, test_case.target,
            test_case.dims, test_case.shift);

        qudit::apply_ctrl_rho_kq_inplace(
            evolved_rho, A, test_case.ctrl, test_case.target,
            test_case.dims, test_case.shift);

        const cmat expected =
            evolved_psi * evolved_psi.adjoint();

        expect_same(evolved_rho, expected);
        EXPECT_NEAR(evolved_rho.trace().real(), 1.0, tolerance);
        EXPECT_NEAR(evolved_rho.trace().imag(), 0.0, tolerance);
        expect_same(evolved_rho, evolved_rho.adjoint());
    }
}

TEST(qpp_qudit_apply_ctrl_rho_kq, InactiveControls) {
    const std::vector<idx> dims{3, 3, 3, 3};
    const cmat A = kron(gt.Xd(3), gt.Xd(3));

    const ket psi = mket({1, 0, 2, 0}, dims);
    const cmat rho = psi * psi.adjoint();

    cmat actual = rho;

    qudit::apply_ctrl_rho_kq_inplace(
        actual, A, {0, 2}, {1, 3}, dims, {0, 0});

    expect_same(actual, rho);
}

TEST(qpp_qudit_apply_ctrl_rho_kq, ShiftedControls) {
    const std::vector<idx> dims{3, 3, 3, 3};
    const cmat A = kron(gt.Xd(3), gt.Xd(3));

    const ket psi = mket({2, 0, 0, 1}, dims);
    const ket expected_psi = mket({2, 1, 0, 2}, dims);

    cmat actual = psi * psi.adjoint();
    const cmat expected =
        expected_psi * expected_psi.adjoint();

    qudit::apply_ctrl_rho_kq_inplace(
        actual, A, {0, 2}, {1, 3}, dims, {1, 2});

    expect_same(actual, expected);
}

TEST(qpp_qudit_apply_ctrl_rho_kq, Identity) {
    const ket psi = make_state(27);
    const cmat rho = psi * psi.adjoint();
    const cmat A = cmat::Identity(9, 9);

    cmat actual = rho;

    qudit::apply_ctrl_rho_kq_inplace(
        actual, A, {0}, {2, 1}, {3, 3, 3}, {2});

    expect_same(actual, rho);
}

TEST(qpp_qudit_apply_ctrl_rho_kq, MaximallyMixedState) {
    const idx D = 27;
    const cmat rho =
        cmat::Identity(D, D) / static_cast<double>(D);

    const cmat A = complex_gate(9);

    cmat actual = rho;

    qudit::apply_ctrl_rho_kq_inplace(
        actual, A, {1}, {2, 0}, {3, 3, 3}, {1});

    expect_same(actual, rho);
}

TEST(qpp_qudit_apply_ctrl_rho_kq, NonUnitaryGate) {
    const Case test_case{{3, 3, 3},
                         {1},
                         {2, 0},
                         {1}};

    const ket psi = make_state(27);
    const cmat rho = psi * psi.adjoint();
    const cmat A = 0.7 * complex_gate(9);

    const cmat U = reference_operator(A, test_case);
    const cmat expected = U * rho * U.adjoint();

    cmat actual = rho;

    qudit::apply_ctrl_rho_kq_inplace(
        actual, A, test_case.ctrl, test_case.target,
        test_case.dims, test_case.shift);

    expect_same(actual, expected);
}

TEST(qpp_qudit_apply_ctrl_rho_kq, NonHermitianOperator) {
    const Case test_case{{3, 3},
                         {0},
                         {1},
                         {2}};

    const ket psi = make_state(9);
    const ket phi = mket({2, 1}, test_case.dims);

    const cmat rho = psi * phi.adjoint();
    const cmat A = complex_gate(3);

    const cmat U = reference_operator(A, test_case);
    const cmat expected = U * rho * U.adjoint();

    cmat actual = rho;

    qudit::apply_ctrl_rho_kq_inplace(
        actual, A, test_case.ctrl, test_case.target,
        test_case.dims, test_case.shift);

    expect_same(actual, expected);
}

TEST(qpp_qudit_apply_ctrl_rho_kq, ZeroState) {
    cmat actual = cmat::Zero(27, 27);
    const cmat A = complex_gate(9);

    qudit::apply_ctrl_rho_kq_inplace(
        actual, A, {0}, {1, 2}, {3, 3, 3}, {1});

    expect_same(actual, cmat::Zero(27, 27));
}

TEST(qpp_qudit_apply_ctrl_psi_kq, OrderedHeterogeneousTargets) {
    const std::vector<idx> dims{2, 3, 4};
    // Target order {2, 0}: local index = 2 * digit[2] + digit[0].
    const cmat A = gt.Xd(8);
    ket actual = mket({1, 1, 2}, dims);
    const ket expected = mket({0, 1, 3}, dims);

    qudit::apply_ctrl_psi_kq_inplace(actual, A, {1}, {2, 0}, dims, {0});

    expect_same(actual, expected);
}

TEST(qpp_qudit_apply_ctrl_rho_kq, OrderedHeterogeneousTargets) {
    const std::vector<idx> dims{2, 3, 4};
    const cmat A = gt.Xd(8);
    const ket psi = mket({1, 1, 2}, dims);
    const ket expected_psi = mket({0, 1, 3}, dims);
    cmat actual = psi * psi.adjoint();
    const cmat expected = expected_psi * expected_psi.adjoint();

    qudit::apply_ctrl_rho_kq_inplace(actual, A, {1}, {2, 0}, dims, {0});

    expect_same(actual, expected);
}
