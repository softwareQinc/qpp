#include "gtest/gtest.h"

#include "qpp/qpp.hpp"

#include "qpp/internal/kernels/qubit/apply_ctrl_diag.hpp"

#include <cmath>
#include <complex>
#include <sstream>
#include <string>
#include <vector>

using namespace qpp;

namespace {

constexpr double tolerance = 1e-11;

struct Case {
    std::vector<idx> ctrl;
    std::vector<idx> target;
    std::vector<idx> shift;
};

std::string describe(const Case& test_case) {
    std::ostringstream out;
    auto print = [&](const char* name, const std::vector<idx>& values) {
        out << name << "={";
        for (idx j = 0; j < values.size(); ++j) {
            if (j != 0) {
                out << ",";
            }
            out << values[j];
        }
        out << "} ";
    };
    print("ctrl", test_case.ctrl);
    print("target", test_case.target);
    print("shift", test_case.shift);
    return out.str();
}

ket make_state(idx n, idx seed = 1) {
    ket result(idx{1} << n);
    for (idx j = 0; j < static_cast<idx>(result.size()); ++j) {
        result(j) = cplx{
            static_cast<double>((7 * j + 3 * seed) % 19) - 9.0,
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

// Modes: 0 = identity, 1 = phases with A(0)=1, 2 = general complex
// non-unitary values, 3 = only A(0) differs from one.
ket diagonal_gate(idx d, idx mode) {
    ket A = ket::Ones(d);
    for (idx j = 0; j < d; ++j) {
        if (mode == 1 && j != 0) {
            A(j) = std::polar(1.0, 0.19 * static_cast<double>(j));
        } else if (mode == 2) {
            A(j) = cplx{0.7 + 0.04 * static_cast<double>(j),
                        0.06 * static_cast<double>(j + 1)};
        } else if (mode == 3 && j == 0) {
            A(j) = cplx{0.8, 0.2};
        }
    }
    return A;
}

// Independent per-basis reference. Controls and targets use big-endian indices;
// shift=0 requires |1>, while shift=1 requires |0>.
ket controlled_diagonal(const ket& A, const Case& test_case, idx n) {
    const idx D = idx{1} << n;
    ket result = ket::Ones(D);
    for (idx basis = 0; basis < D; ++basis) {
        bool enabled = true;
        for (idx j = 0; j < test_case.ctrl.size(); ++j) {
            const idx bit =
                (basis >> (n - 1 - test_case.ctrl[j])) & idx{1};
            if (bit != (idx{1} - test_case.shift[j])) {
                enabled = false;
                break;
            }
        }
        if (!enabled) {
            continue;
        }

        idx local_basis = 0;
        for (idx target : test_case.target) {
            local_basis = (local_basis << 1) |
                          ((basis >> (n - 1 - target)) & idx{1});
        }
        result(basis) = A(local_basis);
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

void expect_psi(const ket& actual, const ket& original, const ket& A,
                const Case& test_case, idx n) {
    const ket factors = controlled_diagonal(A, test_case, n);
    const ket expected = original.cwiseProduct(factors);
    expect_same(actual, expected);
}

void expect_rho(const cmat& actual, const cmat& original, const ket& A,
                const Case& test_case, idx n) {
    const ket factors = controlled_diagonal(A, test_case, n);
    const idx D = idx{1} << n;
    cmat expected(D, D);
    for (idx r = 0; r < D; ++r) {
        for (idx c = 0; c < D; ++c) {
            expected(r, c) =
                factors(r) * original(r, c) * std::conj(factors(c));
        }
    }
    expect_same(actual, expected);
}

} // namespace

TEST(qpp_apply_ctrl_psi_1q_diag_inplace, ControlPatternsAndDiagonalShapes) {
    const idx n = 5;
    const ket original = make_state(n);
    const std::vector<Case> cases{
        {{}, {0}, {}},
        {{4}, {0}, {0}},
        {{4}, {0}, {1}},
        {{0}, {4}, {0}},
        {{0}, {4}, {1}},
        {{1, 3}, {0}, {0, 1}},
        {{1, 3}, {0}, {1, 0}},
        {{4, 0, 2, 3}, {1}, {1, 0, 1, 0}}};

    for (const auto& test_case : cases) {
        SCOPED_TRACE(describe(test_case));
        for (idx mode : {idx{0}, idx{1}, idx{2}, idx{3}}) {
            SCOPED_TRACE(::testing::Message() << "mode=" << mode);
            const ket A = diagonal_gate(2, mode);
            ket actual = original;
            internal::kernels::qubit::apply_ctrl_psi_1q_diag_inplace(
                actual, A, test_case.ctrl, test_case.target[0],
                test_case.shift, n);
            expect_psi(actual, original, A, test_case, n);
        }
    }
}

TEST(qpp_apply_ctrl_psi_1q_diag_inplace, OptimizedPhasePathsLargeState) {
    const idx n = 16;
    const ket original = make_state(n);
    const ket A = diagonal_gate(2, 1);
    // Target 0/control 15: ctrl_bit < step (contiguous fast path).
    // Target 15/control 0: ctrl_bit > step (bit-deposit fallback).
    const std::vector<Case> cases{
        {{15}, {0}, {0}}, {{0}, {15}, {0}}};
    for (const auto& test_case : cases) {
        SCOPED_TRACE(describe(test_case));
        ket actual = original;
        internal::kernels::qubit::apply_ctrl_psi_1q_diag_inplace(
            actual, A, test_case.ctrl, test_case.target[0],
            test_case.shift, n);
        expect_psi(actual, original, A, test_case, n);
    }
}

TEST(qpp_apply_ctrl_psi_1q_diag_inplace, PositiveAndNegativeControlBasis) {
    const ket A = diagonal_gate(2, 1);
    for (idx shift : {idx{0}, idx{1}}) {
        const Case test_case{{0}, {1}, {shift}};
        ket original = ket::Zero(4);
        original(shift == 0 ? 3 : 1) = 1;
        ket actual = original;
        internal::kernels::qubit::apply_ctrl_psi_1q_diag_inplace(
            actual, A, {0}, 1, {shift}, 2);
        expect_psi(actual, original, A, test_case, 2);
    }
}

TEST(qpp_apply_ctrl_psi_1q_diag, ReturnsResultWithoutChangingInput) {
    const idx n = 4;
    ket original = make_state(n);
    const ket before = original;
    const std::vector<Case> cases{
        {{}, {0}, {}},
        {{3}, {0}, {0}},
        {{3}, {0}, {1}},
        {{0}, {3}, {0}},
        {{1, 2}, {0}, {1, 0}}};
    for (const auto& test_case : cases) {
        SCOPED_TRACE(describe(test_case));
        const ket A = diagonal_gate(2, 2);
        const ket actual = internal::kernels::qubit::apply_ctrl_psi_1q_diag(
            original, A, test_case.ctrl, test_case.target[0],
            test_case.shift, n);
        expect_psi(actual, before, A, test_case, n);
        expect_same(original, before);
    }
}

TEST(qpp_apply_ctrl_psi_2q_diag_inplace, ReversedTargetsAndControls) {
    const idx n = 5;
    const ket original = make_state(n);
    const std::vector<Case> cases{
        {{}, {0, 4}, {}},
        {{1}, {0, 4}, {0}},
        {{1}, {0, 4}, {1}},
        {{1, 3}, {0, 4}, {0, 1}},
        {{3, 1}, {4, 0}, {1, 0}},
        {{0, 2, 3}, {4, 1}, {1, 0, 1}},
        {{4}, {2, 3}, {0}}};
    for (const auto& test_case : cases) {
        SCOPED_TRACE(describe(test_case));
        for (idx mode : {idx{0}, idx{1}, idx{2}}) {
            const ket A = diagonal_gate(4, mode);
            ket actual = original;
            internal::kernels::qubit::apply_ctrl_psi_2q_diag_inplace(
                actual, A, test_case.ctrl, test_case.target[0],
                test_case.target[1], test_case.shift, n);
            expect_psi(actual, original, A, test_case, n);
        }
    }
}

TEST(qpp_apply_ctrl_psi_2q_diag_inplace, NoControlsTwoQubits) {
    const ket original = make_state(2);
    const Case test_case{{}, {1, 0}, {}};
    const ket A = diagonal_gate(4, 2);
    ket actual = original;
    internal::kernels::qubit::apply_ctrl_psi_2q_diag_inplace(
        actual, A, {}, 1, 0, {}, 2);
    expect_psi(actual, original, A, test_case, 2);
}

TEST(qpp_apply_ctrl_psi_2q_diag, ReturnsResultWithoutChangingInput) {
    const idx n = 4;
    ket original = make_state(n);
    const ket before = original;
    const std::vector<Case> cases{
        {{}, {0, 3}, {}},
        {{1}, {3, 0}, {0}},
        {{1, 2}, {3, 0}, {1, 0}}};
    for (const auto& test_case : cases) {
        SCOPED_TRACE(describe(test_case));
        const ket A = diagonal_gate(4, 2);
        const ket actual = internal::kernels::qubit::apply_ctrl_psi_2q_diag(
            original, A, test_case.ctrl, test_case.target[0],
            test_case.target[1], test_case.shift, n);
        expect_psi(actual, before, A, test_case, n);
        expect_same(original, before);
    }
}

TEST(qpp_apply_ctrl_psi_kq_diag_inplace, ZeroToFiveTargets) {
    const idx n = 5;
    const ket original = make_state(n);
    const std::vector<Case> cases{
        {{}, {}, {}},
        {{0}, {}, {0}},
        {{0}, {}, {1}},
        {{4}, {0}, {0}},
        {{0}, {4}, {1}},
        {{1, 3}, {0, 4}, {0, 1}},
        {{3, 1}, {4, 0}, {1, 0}},
        {{1, 3}, {4, 2, 0}, {0, 1}},
        {{1}, {4, 2, 0, 3}, {1}},
        {{}, {4, 2, 0, 3, 1}, {}}};
    for (const auto& test_case : cases) {
        SCOPED_TRACE(describe(test_case));
        const idx d = idx{1} << test_case.target.size();
        for (idx mode : {idx{0}, idx{1}, idx{2}}) {
            SCOPED_TRACE(::testing::Message() << "mode=" << mode);
            const ket A = diagonal_gate(d, mode);
            ket actual = original;
            internal::kernels::qubit::apply_ctrl_psi_kq_diag_inplace(
                actual, A, test_case.ctrl, test_case.target,
                test_case.shift, n);
            expect_psi(actual, original, A, test_case, n);
        }
    }
}

TEST(qpp_apply_ctrl_psi_kq_diag, ReturnsResultWithoutChangingInput) {
    const idx n = 5;
    ket original = make_state(n);
    const ket before = original;
    const std::vector<Case> cases{
        {{0}, {}, {0}},
        {{0}, {3}, {1}},
        {{4, 1}, {0, 3}, {0, 1}},
        {{2}, {4, 0, 1}, {1}},
        {{3}, {4, 2, 0, 1}, {0}},
        {{}, {2, 4, 1, 3, 0}, {}}};
    for (const auto& test_case : cases) {
        SCOPED_TRACE(describe(test_case));
        const ket A = diagonal_gate(idx{1} << test_case.target.size(), 2);
        const ket actual = internal::kernels::qubit::apply_ctrl_psi_kq_diag(
            original, A, test_case.ctrl, test_case.target,
            test_case.shift, n);
        expect_psi(actual, before, A, test_case, n);
        expect_same(original, before);
    }
}

TEST(qpp_apply_ctrl_rho_1q_diag_inplace, MixedStateControlBranches) {
    const idx n = 4;
    const cmat original = make_density(n);
    const std::vector<Case> cases{
        {{}, {0}, {}},
        {{3}, {0}, {0}},
        {{3}, {0}, {1}},
        {{0}, {3}, {0}},
        {{0}, {3}, {1}},
        {{1, 3}, {0}, {0, 1}},
        {{1, 2, 3}, {0}, {1, 0, 1}}};
    for (const auto& test_case : cases) {
        SCOPED_TRACE(describe(test_case));
        for (idx mode : {idx{0}, idx{1}, idx{2}, idx{3}}) {
            const ket A = diagonal_gate(2, mode);
            cmat actual = original;
            internal::kernels::qubit::apply_ctrl_rho_1q_diag_inplace(
                actual, A, test_case.ctrl, test_case.target[0],
                test_case.shift, n);
            expect_rho(actual, original, A, test_case, n);
        }
    }
}

TEST(qpp_apply_ctrl_rho_1q_diag_inplace, OffDiagonalCoherence) {
    const idx n = 2;
    ket psi = ket::Zero(4);
    psi(0) = 1.0 / std::sqrt(2.0);
    psi(3) = 1.0 / std::sqrt(2.0);
    const cmat original = psi * psi.adjoint();
    const Case test_case{{0}, {1}, {0}};
    const ket A = diagonal_gate(2, 1);
    cmat actual = original;

    internal::kernels::qubit::apply_ctrl_rho_1q_diag_inplace(
        actual, A, test_case.ctrl, 1, test_case.shift, n);

    expect_rho(actual, original, A, test_case, n);
    EXPECT_NEAR(std::abs(actual(0, 3)), 0.5, tolerance);
    EXPECT_NEAR(std::abs(actual(0, 3) - original(0, 3)), 
                std::abs((std::conj(A(1)) - cplx{1}) * original(0, 3)),
                tolerance);
}

TEST(qpp_apply_ctrl_rho_1q_diag, ReturnsResultWithoutChangingInput) {
    const idx n = 4;
    cmat original = make_density(n);
    const cmat before = original;
    const std::vector<Case> cases{
        {{}, {0}, {}},
        {{0}, {3}, {0}},
        {{0}, {3}, {1}},
        {{1, 2}, {0}, {0, 1}}};
    for (const auto& test_case : cases) {
        SCOPED_TRACE(describe(test_case));
        const ket A = diagonal_gate(2, 2);
        const cmat actual = internal::kernels::qubit::apply_ctrl_rho_1q_diag(
            original, A, test_case.ctrl, test_case.target[0],
            test_case.shift, n);
        expect_rho(actual, before, A, test_case, n);
        expect_same(original, before);
    }
}

TEST(qpp_apply_ctrl_rho_2q_diag_inplace, ReversedTargetsAndMixedControls) {
    const idx n = 5;
    const cmat original = make_density(n);
    const std::vector<Case> cases{
        {{}, {0, 4}, {}},
        {{1}, {0, 4}, {0}},
        {{1}, {0, 4}, {1}},
        {{1, 3}, {0, 4}, {0, 1}},
        {{3, 1}, {4, 0}, {1, 0}},
        {{0, 2, 3}, {4, 1}, {1, 0, 1}},
        {{4}, {2, 3}, {0}}};
    for (const auto& test_case : cases) {
        SCOPED_TRACE(describe(test_case));
        for (idx mode : {idx{0}, idx{1}, idx{2}}) {
            const ket A = diagonal_gate(4, mode);
            cmat actual = original;
            internal::kernels::qubit::apply_ctrl_rho_2q_diag_inplace(
                actual, A, test_case.ctrl, test_case.target[0],
                test_case.target[1], test_case.shift, n);
            expect_rho(actual, original, A, test_case, n);
        }
    }
}

TEST(qpp_apply_ctrl_rho_2q_diag_inplace, NoControlsTwoQubits) {
    const cmat original = make_density(2);
    const Case test_case{{}, {1, 0}, {}};
    const ket A = diagonal_gate(4, 2);
    cmat actual = original;
    internal::kernels::qubit::apply_ctrl_rho_2q_diag_inplace(
        actual, A, {}, 1, 0, {}, 2);
    expect_rho(actual, original, A, test_case, 2);
}

TEST(qpp_apply_ctrl_rho_2q_diag, ReturnsResultWithoutChangingInput) {
    const idx n = 4;
    cmat original = make_density(n);
    const cmat before = original;
    const std::vector<Case> cases{
        {{}, {0, 3}, {}},
        {{1}, {3, 0}, {0}},
        {{1, 2}, {3, 0}, {1, 0}}};
    for (const auto& test_case : cases) {
        SCOPED_TRACE(describe(test_case));
        const ket A = diagonal_gate(4, 2);
        const cmat actual = internal::kernels::qubit::apply_ctrl_rho_2q_diag(
            original, A, test_case.ctrl, test_case.target[0],
            test_case.target[1], test_case.shift, n);
        expect_rho(actual, before, A, test_case, n);
        expect_same(original, before);
    }
}

TEST(qpp_apply_ctrl_rho_kq_diag_inplace, ZeroToFiveTargets) {
    const idx n = 5;
    const cmat original = make_density(n);
    const std::vector<Case> cases{
        {{}, {}, {}},
        {{0}, {}, {0}},
        {{0}, {}, {1}},
        {{4}, {0}, {0}},
        {{0}, {4}, {1}},
        {{1, 3}, {0, 4}, {0, 1}},
        {{3, 1}, {4, 0}, {1, 0}},
        {{1, 3}, {4, 2, 0}, {0, 1}},
        {{1}, {4, 2, 0, 3}, {1}},
        {{}, {4, 2, 0, 3, 1}, {}}};
    for (const auto& test_case : cases) {
        SCOPED_TRACE(describe(test_case));
        const idx d = idx{1} << test_case.target.size();
        for (idx mode : {idx{0}, idx{1}, idx{2}}) {
            const ket A = diagonal_gate(d, mode);
            cmat actual = original;
            internal::kernels::qubit::apply_ctrl_rho_kq_diag_inplace(
                actual, A, test_case.ctrl, test_case.target,
                test_case.shift, n);
            expect_rho(actual, original, A, test_case, n);
        }
    }
}

TEST(qpp_apply_ctrl_rho_kq_diag, ReturnsResultWithoutChangingInput) {
    const idx n = 5;
    cmat original = make_density(n);
    const cmat before = original;
    const std::vector<Case> cases{
        {{0}, {}, {0}},
        {{0}, {3}, {1}},
        {{4, 1}, {0, 3}, {0, 1}},
        {{2}, {4, 0, 1}, {1}},
        {{3}, {4, 2, 0, 1}, {0}},
        {{}, {2, 4, 1, 3, 0}, {}}};
    for (const auto& test_case : cases) {
        SCOPED_TRACE(describe(test_case));
        const ket A =
            diagonal_gate(idx{1} << test_case.target.size(), 2);
        const cmat actual = internal::kernels::qubit::apply_ctrl_rho_kq_diag(
            original, A, test_case.ctrl, test_case.target,
            test_case.shift, n);
        expect_rho(actual, before, A, test_case, n);
        expect_same(original, before);
    }
}

TEST(qpp_apply_ctrl_rho_kq_diag, PureStateConsistencyAndTrace) {
    const idx n = 5;
    const ket psi = make_state(n);
    const cmat rho = psi * psi.adjoint();
    const Case test_case{{1, 4}, {3, 0, 2}, {0, 1}};
    const ket A = diagonal_gate(8, 1);

    const ket transformed_psi =
        internal::kernels::qubit::apply_ctrl_psi_kq_diag(
            psi, A, test_case.ctrl, test_case.target, test_case.shift, n);
    const cmat transformed_rho =
        internal::kernels::qubit::apply_ctrl_rho_kq_diag(
            rho, A, test_case.ctrl, test_case.target, test_case.shift, n);
    const cmat expected = transformed_psi * transformed_psi.adjoint();

    expect_same(transformed_rho, expected);
    EXPECT_NEAR(transformed_rho.trace().real(), 1.0, tolerance);
    expect_same(transformed_rho, transformed_rho.adjoint());
}
