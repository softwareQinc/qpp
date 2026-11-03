#include "gtest/gtest.h"

#include "qpp/qpp.hpp"

#include "qpp/internal/kernels/qubit/apply_ctrl.hpp"

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
    const idx D = idx{1} << n;
    ket result(D);
    for (idx i = 0; i < D; ++i) {
        result(i) = cplx{static_cast<double>((7 * i + 3 * seed) % 19) - 9.0,
                         static_cast<double>((5 * i + 2 * seed) % 17) - 8.0};
    }
    result.normalize();
    return result;
}

cmat make_density(idx n) {
    const ket psi = make_state(n, 1);
    const ket phi = make_state(n, 3);
    return 0.4 * (psi * psi.adjoint()) + 0.6 * (phi * phi.adjoint());
}

cmat fourier_gate(idx d) {
    cmat result(d, d);
    const double pi = std::acos(-1.0);
    const double scale = 1.0 / std::sqrt(static_cast<double>(d));
    for (idx r = 0; r < d; ++r) {
        for (idx c = 0; c < d; ++c) {
            const double angle =
                2.0 * pi * static_cast<double>(r * c) / static_cast<double>(d) +
                0.17 * static_cast<double>(r) - 0.23 * static_cast<double>(c);
            result(r, c) = std::polar(scale, angle);
        }
    }
    return result;
}

cmat dense_gate(idx d) {
    cmat result(d, d);
    for (idx r = 0; r < d; ++r) {
        for (idx c = 0; c < d; ++c) {
            result(r, c) =
                cplx{(r == c ? 0.8 : 0.0) + 0.07 * (r + 1) - 0.03 * (c + 1),
                     0.04 * (r + 1) + 0.02 * (c + 1)};
        }
    }
    return result;
}

// Independent reference: start with the full identity and replace only the
// columns whose computational-basis control bits satisfy the requested shifts.
// A shift of 0 requires |1>; a shift of 1 requires |0>.
cmat controlled_operator(const cmat& A, const Case& test_case, idx n) {
    const idx D = idx{1} << n;
    const idx k = test_case.target.size();
    cmat result = cmat::Identity(D, D);

    for (idx input = 0; input < D; ++input) {
        bool enabled = true;
        for (idx j = 0; j < test_case.ctrl.size(); ++j) {
            const idx bit = (input >> (n - 1 - test_case.ctrl[j])) & idx{1};
            if (bit != (idx{1} - test_case.shift[j])) {
                enabled = false;
                break;
            }
        }
        if (!enabled) {
            continue;
        }

        idx local_input = 0;
        for (idx q : test_case.target) {
            local_input =
                (local_input << 1) | ((input >> (n - 1 - q)) & idx{1});
        }

        result(input, input) = 0;
        for (idx local_output = 0; local_output < (idx{1} << k);
             ++local_output) {
            idx output = input;
            for (idx j = 0; j < k; ++j) {
                const idx mask = idx{1} << (n - 1 - test_case.target[j]);
                const idx bit = (local_output >> (k - 1 - j)) & idx{1};
                output = bit ? (output | mask) : (output & ~mask);
            }
            result(output, input) = A(local_output, local_input);
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

void expect_psi(const ket& actual, const ket& original, const cmat& gate,
                const Case& test_case, idx n) {
    const ket expected = controlled_operator(gate, test_case, n) * original;
    expect_same(actual, expected);
}

void expect_rho(const cmat& actual, const cmat& original, const cmat& gate,
                const Case& test_case, idx n) {
    const cmat U = controlled_operator(gate, test_case, n);
    const cmat expected = U * original * U.adjoint();
    expect_same(actual, expected);
}

} // namespace

TEST(qpp_qubit_apply_ctrl_psi_1q_inplace, PositiveNegativeAndMixedControls) {
    const idx n = 5;
    const ket original = make_state(n);
    const std::vector<Case> cases{
        {{}, {0}, {}},         {{4}, {0}, {0}},
        {{4}, {0}, {1}},       {{0}, {4}, {0}},
        {{0}, {4}, {1}},       {{1, 3}, {0}, {0, 1}},
        {{1, 3}, {0}, {1, 0}}, {{4, 0, 2, 3}, {1}, {1, 0, 1, 0}}};

    for (const auto& test_case : cases) {
        SCOPED_TRACE(describe(test_case));
        for (const cmat& gate : {cmat(gt.X), fourier_gate(2), dense_gate(2)}) {
            ket actual = original;
            internal::kernels::qubit::apply_ctrl_psi_1q_inplace(
                actual, gate, test_case.ctrl, test_case.target[0],
                test_case.shift, n);
            expect_psi(actual, original, gate, test_case, n);
        }
    }
}

TEST(qpp_qubit_apply_ctrl_psi_1q_inplace, PositiveControlCreatesBellState) {
    ket state = ket::Zero(4);
    state(0) = 1.0 / std::sqrt(2.0);
    state(2) = 1.0 / std::sqrt(2.0);

    internal::kernels::qubit::apply_ctrl_psi_1q_inplace(state, gt.X, {0}, 1,
                                                        {0}, 2);

    ket expected = ket::Zero(4);
    expected(0) = 1.0 / std::sqrt(2.0);
    expected(3) = 1.0 / std::sqrt(2.0);
    expect_same(state, expected);
}

TEST(qpp_qubit_apply_ctrl_psi_1q_inplace, NegativeControlOnZeroBranch) {
    ket state = ket::Zero(4);
    state(0) = 1;
    internal::kernels::qubit::apply_ctrl_psi_1q_inplace(state, gt.X, {0}, 1,
                                                        {1}, 2);

    ket expected = ket::Zero(4);
    expected(1) = 1;
    expect_same(state, expected);
}

TEST(qpp_qubit_apply_ctrl_psi_1q_inplace, NoControlsOnSingleQubit) {
    const ket original = make_state(1);
    const Case test_case{{}, {0}, {}};
    const cmat gate = dense_gate(2);
    ket actual = original;
    internal::kernels::qubit::apply_ctrl_psi_1q_inplace(actual, gate, {}, 0, {},
                                                        1);
    expect_psi(actual, original, gate, test_case, 1);
}

TEST(qpp_qubit_apply_ctrl_psi_1q, ReturnsResultWithoutChangingInput) {
    const idx n = 4;
    ket original = make_state(n);
    const ket before = original;
    const cmat gate = fourier_gate(2);
    const std::vector<Case> cases{
        {{}, {0}, {}}, {{0}, {3}, {0}}, {{0}, {3}, {1}}, {{1, 2}, {0}, {0, 1}}};
    for (const auto& test_case : cases) {
        SCOPED_TRACE(describe(test_case));
        const ket actual = internal::kernels::qubit::apply_ctrl_psi_1q(
            original, gate, test_case.ctrl, test_case.target[0],
            test_case.shift, n);
        expect_psi(actual, before, gate, test_case, n);
        expect_same(original, before);
    }
}

TEST(qpp_qubit_apply_ctrl_psi_2q_inplace, ReversedTargetsAndMixedControls) {
    const idx n = 5;
    const ket original = make_state(n);
    const std::vector<Case> cases{
        {{}, {0, 4}, {}},         {{1}, {0, 4}, {0}},
        {{1}, {0, 4}, {1}},       {{1, 3}, {0, 4}, {0, 1}},
        {{3, 1}, {4, 0}, {1, 0}}, {{0, 2, 3}, {4, 1}, {1, 0, 1}},
        {{4}, {2, 3}, {0}}};

    for (const auto& test_case : cases) {
        SCOPED_TRACE(describe(test_case));
        for (const cmat& gate : {fourier_gate(4), dense_gate(4)}) {
            ket actual = original;
            internal::kernels::qubit::apply_ctrl_psi_2q_inplace(
                actual, gate, test_case.ctrl, test_case.target[0],
                test_case.target[1], test_case.shift, n);
            expect_psi(actual, original, gate, test_case, n);
        }
    }
}

TEST(qpp_qubit_apply_ctrl_psi_2q_inplace, NoControlsOnTwoQubits) {
    const idx n = 2;
    const ket original = make_state(n);
    const Case test_case{{}, {1, 0}, {}};
    const cmat gate = dense_gate(4);
    ket actual = original;
    internal::kernels::qubit::apply_ctrl_psi_2q_inplace(actual, gate, {}, 1, 0,
                                                        {}, n);
    expect_psi(actual, original, gate, test_case, n);
}

TEST(qpp_qubit_apply_ctrl_psi_2q, ReturnsResultWithoutChangingInput) {
    const idx n = 4;
    ket original = make_state(n);
    const ket before = original;
    const cmat gate = fourier_gate(4);
    const std::vector<Case> cases{
        {{}, {0, 3}, {}}, {{1}, {3, 0}, {0}}, {{1, 2}, {3, 0}, {1, 0}}};
    for (const auto& test_case : cases) {
        SCOPED_TRACE(describe(test_case));
        const ket actual = internal::kernels::qubit::apply_ctrl_psi_2q(
            original, gate, test_case.ctrl, test_case.target[0],
            test_case.target[1], test_case.shift, n);
        expect_psi(actual, before, gate, test_case, n);
        expect_same(original, before);
    }
}

TEST(qpp_qubit_apply_ctrl_psi_kq_inplace, ZeroToFiveTargets) {
    const idx n = 5;
    const ket original = make_state(n);
    const std::vector<Case> cases{{{}, {}, {}},
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
        for (const cmat& gate : {fourier_gate(d), dense_gate(d)}) {
            ket actual = original;
            internal::kernels::qubit::apply_ctrl_psi_kq_inplace(
                actual, gate, test_case.ctrl, test_case.target, test_case.shift,
                n);
            expect_psi(actual, original, gate, test_case, n);
        }
    }
}

TEST(qpp_qubit_apply_ctrl_psi_kq, ReturnsResultWithoutChangingInput) {
    const idx n = 5;
    ket original = make_state(n);
    const ket before = original;
    const std::vector<Case> cases{
        {{0}, {}, {0}},           {{0}, {3}, {1}},
        {{4, 1}, {0, 3}, {0, 1}}, {{2}, {4, 0, 1}, {1}},
        {{3}, {4, 2, 0, 1}, {0}}, {{}, {2, 4, 1, 3, 0}, {}}};
    for (const auto& test_case : cases) {
        SCOPED_TRACE(describe(test_case));
        const cmat gate = fourier_gate(idx{1} << test_case.target.size());
        const ket actual = internal::kernels::qubit::apply_ctrl_psi_kq(
            original, gate, test_case.ctrl, test_case.target, test_case.shift,
            n);
        expect_psi(actual, before, gate, test_case, n);
        expect_same(original, before);
    }
}

TEST(qpp_qubit_apply_ctrl_rho_1q_inplace, MixedStateControlBranches) {
    const idx n = 4;
    const cmat original = make_density(n);
    const std::vector<Case> cases{{{}, {0}, {}},
                                  {{3}, {0}, {0}},
                                  {{3}, {0}, {1}},
                                  {{0}, {3}, {0}},
                                  {{0}, {3}, {1}},
                                  {{1, 3}, {0}, {0, 1}},
                                  {{1, 2, 3}, {0}, {1, 0, 1}}};
    for (const auto& test_case : cases) {
        SCOPED_TRACE(describe(test_case));
        for (const cmat& gate : {fourier_gate(2), dense_gate(2)}) {
            cmat actual = original;
            internal::kernels::qubit::apply_ctrl_rho_1q_inplace(
                actual, gate, test_case.ctrl, test_case.target[0],
                test_case.shift, n);
            expect_rho(actual, original, gate, test_case, n);
        }
    }
}

TEST(qpp_qubit_apply_ctrl_rho_1q_inplace, CoherenceAcrossControlBranches) {
    const idx n = 2;
    ket psi = ket::Zero(4);
    psi(0) = 1.0 / std::sqrt(2.0);
    psi(2) = 1.0 / std::sqrt(2.0);
    cmat actual = psi * psi.adjoint();
    const Case test_case{{0}, {1}, {0}};

    internal::kernels::qubit::apply_ctrl_rho_1q_inplace(
        actual, gt.X, test_case.ctrl, 1, test_case.shift, n);

    ket bell = ket::Zero(4);
    bell(0) = 1.0 / std::sqrt(2.0);
    bell(3) = 1.0 / std::sqrt(2.0);
    const cmat expected = bell * bell.adjoint();
    expect_same(actual, expected);
    EXPECT_NEAR(std::abs(actual(0, 3)), 0.5, tolerance);
}

TEST(qpp_qubit_apply_ctrl_rho_1q_inplace, NoControlsOnSingleQubit) {
    const cmat original = make_density(1);
    const Case test_case{{}, {0}, {}};
    const cmat gate = dense_gate(2);
    cmat actual = original;
    internal::kernels::qubit::apply_ctrl_rho_1q_inplace(actual, gate, {}, 0, {},
                                                        1);
    expect_rho(actual, original, gate, test_case, 1);
}

TEST(qpp_qubit_apply_ctrl_rho_1q, ReturnsResultWithoutChangingInput) {
    const idx n = 4;
    cmat original = make_density(n);
    const cmat before = original;
    const cmat gate = fourier_gate(2);
    const std::vector<Case> cases{
        {{}, {0}, {}}, {{0}, {3}, {0}}, {{0}, {3}, {1}}, {{1, 2}, {0}, {0, 1}}};
    for (const auto& test_case : cases) {
        SCOPED_TRACE(describe(test_case));
        const cmat actual = internal::kernels::qubit::apply_ctrl_rho_1q(
            original, gate, test_case.ctrl, test_case.target[0],
            test_case.shift, n);
        expect_rho(actual, before, gate, test_case, n);
        expect_same(original, before);
    }
}

TEST(qpp_qubit_apply_ctrl_rho_2q_inplace, ReversedTargetsAndMixedControls) {
    const idx n = 5;
    const cmat original = make_density(n);
    const std::vector<Case> cases{
        {{}, {0, 4}, {}},         {{1}, {0, 4}, {0}},
        {{1}, {0, 4}, {1}},       {{1, 3}, {0, 4}, {0, 1}},
        {{3, 1}, {4, 0}, {1, 0}}, {{0, 2, 3}, {4, 1}, {1, 0, 1}},
        {{4}, {2, 3}, {0}}};
    for (const auto& test_case : cases) {
        SCOPED_TRACE(describe(test_case));
        for (const cmat& gate : {fourier_gate(4), dense_gate(4)}) {
            cmat actual = original;
            internal::kernels::qubit::apply_ctrl_rho_2q_inplace(
                actual, gate, test_case.ctrl, test_case.target[0],
                test_case.target[1], test_case.shift, n);
            expect_rho(actual, original, gate, test_case, n);
        }
    }
}

TEST(qpp_qubit_apply_ctrl_rho_2q_inplace, NoControlsOnTwoQubits) {
    const idx n = 2;
    const cmat original = make_density(n);
    const Case test_case{{}, {1, 0}, {}};
    const cmat gate = dense_gate(4);
    cmat actual = original;
    internal::kernels::qubit::apply_ctrl_rho_2q_inplace(actual, gate, {}, 1, 0,
                                                        {}, n);
    expect_rho(actual, original, gate, test_case, n);
}

TEST(qpp_qubit_apply_ctrl_rho_2q, ReturnsResultWithoutChangingInput) {
    const idx n = 4;
    cmat original = make_density(n);
    const cmat before = original;
    const cmat gate = fourier_gate(4);
    const std::vector<Case> cases{
        {{}, {0, 3}, {}}, {{1}, {3, 0}, {0}}, {{1, 2}, {3, 0}, {1, 0}}};
    for (const auto& test_case : cases) {
        SCOPED_TRACE(describe(test_case));
        const cmat actual = internal::kernels::qubit::apply_ctrl_rho_2q(
            original, gate, test_case.ctrl, test_case.target[0],
            test_case.target[1], test_case.shift, n);
        expect_rho(actual, before, gate, test_case, n);
        expect_same(original, before);
    }
}

TEST(qpp_qubit_apply_ctrl_rho_kq_inplace, ZeroToFiveTargets) {
    const idx n = 5;
    const cmat original = make_density(n);
    const std::vector<Case> cases{{{}, {}, {}},
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
        for (const cmat& gate : {fourier_gate(d), dense_gate(d)}) {
            cmat actual = original;
            internal::kernels::qubit::apply_ctrl_rho_kq_inplace(
                actual, gate, test_case.ctrl, test_case.target, test_case.shift,
                n);
            expect_rho(actual, original, gate, test_case, n);
        }
    }
}

TEST(qpp_qubit_apply_ctrl_rho_kq, ReturnsResultWithoutChangingInput) {
    const idx n = 5;
    cmat original = make_density(n);
    const cmat before = original;
    const std::vector<Case> cases{
        {{0}, {}, {0}},           {{0}, {3}, {1}},
        {{4, 1}, {0, 3}, {0, 1}}, {{2}, {4, 0, 1}, {1}},
        {{3}, {4, 2, 0, 1}, {0}}, {{}, {2, 4, 1, 3, 0}, {}}};
    for (const auto& test_case : cases) {
        SCOPED_TRACE(describe(test_case));
        const cmat gate = fourier_gate(idx{1} << test_case.target.size());
        const cmat actual = internal::kernels::qubit::apply_ctrl_rho_kq(
            original, gate, test_case.ctrl, test_case.target, test_case.shift,
            n);
        expect_rho(actual, before, gate, test_case, n);
        expect_same(original, before);
    }
}

TEST(qpp_qubit_apply_ctrl_rho_kq, PureStateConsistencyAndTrace) {
    const idx n = 5;
    const ket psi = make_state(n);
    const cmat rho = psi * psi.adjoint();
    const Case test_case{{1, 4}, {3, 0, 2}, {0, 1}};
    const cmat gate = fourier_gate(8);

    const ket transformed_psi = internal::kernels::qubit::apply_ctrl_psi_kq(
        psi, gate, test_case.ctrl, test_case.target, test_case.shift, n);
    const cmat transformed_rho = internal::kernels::qubit::apply_ctrl_rho_kq(
        rho, gate, test_case.ctrl, test_case.target, test_case.shift, n);
    const cmat expected = transformed_psi * transformed_psi.adjoint();

    expect_same(transformed_rho, expected);
    EXPECT_NEAR(transformed_rho.trace().real(), 1.0, tolerance);
    expect_same(transformed_rho, transformed_rho.adjoint());
}
