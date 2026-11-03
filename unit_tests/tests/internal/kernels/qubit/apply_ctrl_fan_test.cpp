#include "gtest/gtest.h"

#include "qpp/qpp.hpp"

#include "qpp/internal/kernels/qubit/apply_ctrl_fan.hpp"

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
    auto print = [&](const char* label, const std::vector<idx>& indices) {
        out << label << "={";
        for (idx j = 0; j < indices.size(); ++j) {
            if (j != 0) {
                out << ",";
            }
            out << indices[j];
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
    return 0.4 * (psi * psi.adjoint()) + 0.6 * (phi * phi.adjoint());
}

cmat complex_gate() {
    cmat A(2, 2);
    const double scale = 1.0 / std::sqrt(2.0);
    const cplx p0 = std::polar(1.0, 0.27);
    const cplx p1 = std::polar(1.0, -0.31);
    A << scale * p0, scale * p1, scale * p0, -scale * p1;
    return A;
}

cmat dense_gate() {
    cmat A(2, 2);
    A << cplx{0.8, 0.2}, cplx{0.3, -0.1}, cplx{-0.2, 0.4},
        cplx{0.7, -0.3};
    return A;
}

// Independently construct one controlled single-target gate in the full
// computational basis. This does not use any of the internal kernels.
cmat single_target_operator(const cmat& A, const Case& test_case, idx target,
                            idx n) {
    const idx D = idx{1} << n;
    const idx target_bit = idx{1} << (n - 1 - target);
    cmat U = cmat::Identity(D, D);

    for (idx input = 0; input < D; ++input) {
        bool enabled = true;
        for (idx j = 0; j < test_case.ctrl.size(); ++j) {
            const idx control_bit =
                (input >> (n - 1 - test_case.ctrl[j])) & idx{1};
            if (control_bit != (idx{1} - test_case.shift[j])) {
                enabled = false;
                break;
            }
        }
        if (!enabled) {
            continue;
        }

        const idx local_input = (input & target_bit) ? 1 : 0;
        const idx base = input & ~target_bit;
        U(input, input) = 0;
        U(base, input) = A(0, local_input);
        U(base | target_bit, input) = A(1, local_input);
    }
    return U;
}

// A fan applies the same controlled A to every target.
cmat full_fan_operator(const cmat& A, const Case& test_case, idx n) {
    const idx D = idx{1} << n;
    cmat U = cmat::Identity(D, D);
    for (idx target : test_case.target) {
        U = single_target_operator(A, test_case, target, n) * U;
    }
    return U;
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

void expect_psi(const ket& actual, const ket& original, const cmat& A,
                const Case& test_case, idx n) {
    const ket expected = full_fan_operator(A, test_case, n) * original;
    expect_same(actual, expected);
}

void expect_rho(const cmat& actual, const cmat& original, const cmat& A,
                const Case& test_case, idx n) {
    const cmat U = full_fan_operator(A, test_case, n);
    const cmat expected = U * original * U.adjoint();
    expect_same(actual, expected);
}

} // namespace

TEST(qpp_qubit_apply_ctrl_fan_psi_inplace, ControlPatternsAndMultipleTargets) {
    const idx n = 5;
    const ket original = make_state(n);
    const std::vector<Case> cases{
        {{}, {}, {}},
        {{0}, {}, {0}},
        {{}, {0}, {}},
        {{}, {4, 1, 3}, {}},
        {{0}, {1, 3, 4}, {0}},
        {{0}, {1, 3, 4}, {1}},
        {{1, 4}, {0, 2, 3}, {0, 1}},
        {{1, 4}, {3, 2, 0}, {1, 0}},
        {{0, 2, 4}, {1, 3}, {1, 0, 1}},
        {{0, 1, 2, 3}, {4}, {0, 0, 0, 0}},
        {{}, {4, 3, 2, 1, 0}, {}}};

    for (const auto& test_case : cases) {
        SCOPED_TRACE(describe(test_case));
        for (const cmat& A : {cmat(gt.X), cmat(gt.H), complex_gate(),
                              dense_gate()}) {
            ket actual = original;
            internal::kernels::qubit::apply_ctrl_fan_psi_inplace(
                actual, A, test_case.ctrl, test_case.target,
                test_case.shift, n);
            expect_psi(actual, original, A, test_case, n);
        }
    }
}

TEST(qpp_qubit_apply_ctrl_fan_psi_inplace, TwoControlsTwoTargets) {
    // With both controls on, X must flip each target: |1100> -> |1111>.
    const idx n = 4;
    ket actual = ket::Zero(16);
    actual(12) = 1;

    internal::kernels::qubit::apply_ctrl_fan_psi_inplace(
        actual, gt.X, {0, 1}, {2, 3}, {0, 0}, n);

    ket expected = ket::Zero(16);
    expected(15) = 1;
    expect_same(actual, expected);
}

TEST(qpp_qubit_apply_ctrl_fan_psi_inplace, ControlsNotSatisfied) {
    const idx n = 4;
    ket actual = ket::Zero(16);
    actual(8) = 1; // |1000>: controls {0,1} are not both |1>.
    const ket original = actual;

    internal::kernels::qubit::apply_ctrl_fan_psi_inplace(
        actual, gt.X, {0, 1}, {2, 3}, {0, 0}, n);

    expect_same(actual, original);
}

TEST(qpp_qubit_apply_ctrl_fan_psi_inplace, NegativeControls) {
    const idx n = 4;
    ket actual = ket::Zero(16);
    actual(8) = 1; // |1000>: q0=1 and q1=0.

    internal::kernels::qubit::apply_ctrl_fan_psi_inplace(
        actual, gt.X, {0, 1}, {2, 3}, {0, 1}, n);

    ket expected = ket::Zero(16);
    expected(11) = 1; // |1011>
    expect_same(actual, expected);
}

TEST(qpp_qubit_apply_ctrl_fan_psi_inplace, EmptyAndSingleQubitSystems) {
    const cmat A = complex_gate();
    for (idx n : {idx{0}, idx{1}}) {
        const ket original = make_state(n);
        const Case test_case =
            (n == 0) ? Case{{}, {}, {}} : Case{{}, {0}, {}};
        ket actual = original;
        internal::kernels::qubit::apply_ctrl_fan_psi_inplace(
            actual, A, test_case.ctrl, test_case.target, test_case.shift, n);
        expect_psi(actual, original, A, test_case, n);
    }
}

TEST(qpp_qubit_apply_ctrl_fan_psi, ReturnsResultWithoutChangingInput) {
    const idx n = 5;
    ket original = make_state(n);
    const ket before = original;
    const cmat A = complex_gate();
    const std::vector<Case> cases{
        {{}, {}, {}},
        {{0}, {2, 4}, {0}},
        {{0}, {4, 2}, {1}},
        {{0, 3}, {1, 2, 4}, {0, 1}},
        {{1, 4}, {3, 0, 2}, {1, 0}}};
    for (const auto& test_case : cases) {
        SCOPED_TRACE(describe(test_case));
        const ket actual = internal::kernels::qubit::apply_ctrl_fan_psi(
            original, A, test_case.ctrl, test_case.target,
            test_case.shift, n);
        expect_psi(actual, before, A, test_case, n);
        expect_same(original, before);
    }
}

TEST(qpp_qubit_apply_ctrl_fan_rho_inplace, MixedStateControlBranches) {
    const idx n = 5;
    const cmat original = make_density(n);
    const std::vector<Case> cases{
        {{}, {}, {}},
        {{0}, {}, {0}},
        {{}, {0}, {}},
        {{}, {4, 1, 3}, {}},
        {{0}, {1, 3, 4}, {0}},
        {{0}, {1, 3, 4}, {1}},
        {{1, 4}, {0, 2, 3}, {0, 1}},
        {{1, 4}, {3, 2, 0}, {1, 0}},
        {{0, 2, 4}, {1, 3}, {1, 0, 1}},
        {{0, 1, 2, 3}, {4}, {0, 0, 0, 0}},
        {{}, {4, 3, 2, 1, 0}, {}}};

    for (const auto& test_case : cases) {
        SCOPED_TRACE(describe(test_case));
        for (const cmat& A : {cmat(gt.X), complex_gate(), dense_gate()}) {
            cmat actual = original;
            internal::kernels::qubit::apply_ctrl_fan_rho_inplace(
                actual, A, test_case.ctrl, test_case.target,
                test_case.shift, n);
            expect_rho(actual, original, A, test_case, n);
        }
    }
}

TEST(qpp_qubit_apply_ctrl_fan_rho_inplace, CoherenceAcrossControlBranches) {
    const idx n = 4;
    ket psi = ket::Zero(16);
    psi(0) = 1.0 / std::sqrt(2.0);
    psi(12) = 1.0 / std::sqrt(2.0);
    cmat actual = psi * psi.adjoint();

    internal::kernels::qubit::apply_ctrl_fan_rho_inplace(
        actual, gt.X, {0, 1}, {2, 3}, {0, 0}, n);

    ket expected_psi = ket::Zero(16);
    expected_psi(0) = 1.0 / std::sqrt(2.0);
    expected_psi(15) = 1.0 / std::sqrt(2.0);
    const cmat expected = expected_psi * expected_psi.adjoint();

    expect_same(actual, expected);
    EXPECT_NEAR(std::abs(actual(0, 15)), 0.5, tolerance);
}

TEST(qpp_qubit_apply_ctrl_fan_rho_inplace, EmptyAndSingleQubitSystems) {
    const cmat A = complex_gate();
    for (idx n : {idx{0}, idx{1}}) {
        const cmat original = make_density(n);
        const Case test_case =
            (n == 0) ? Case{{}, {}, {}} : Case{{}, {0}, {}};
        cmat actual = original;
        internal::kernels::qubit::apply_ctrl_fan_rho_inplace(
            actual, A, test_case.ctrl, test_case.target, test_case.shift, n);
        expect_rho(actual, original, A, test_case, n);
    }
}

TEST(qpp_qubit_apply_ctrl_fan_rho, ReturnsResultWithoutChangingInput) {
    const idx n = 5;
    cmat original = make_density(n);
    const cmat before = original;
    const cmat A = complex_gate();
    const std::vector<Case> cases{
        {{}, {}, {}},
        {{0}, {2, 4}, {0}},
        {{0}, {4, 2}, {1}},
        {{0, 3}, {1, 2, 4}, {0, 1}},
        {{1, 4}, {3, 0, 2}, {1, 0}}};
    for (const auto& test_case : cases) {
        SCOPED_TRACE(describe(test_case));
        const cmat actual = internal::kernels::qubit::apply_ctrl_fan_rho(
            original, A, test_case.ctrl, test_case.target,
            test_case.shift, n);
        expect_rho(actual, before, A, test_case, n);
        expect_same(original, before);
    }
}

TEST(qpp_qubit_apply_ctrl_fan_rho, PureStateConsistencyAndTrace) {
    const idx n = 5;
    const ket psi = make_state(n);
    const cmat rho = psi * psi.adjoint();
    const Case test_case{{1, 4}, {3, 0, 2}, {0, 1}};
    const cmat A = complex_gate();

    const ket evolved_psi = internal::kernels::qubit::apply_ctrl_fan_psi(
        psi, A, test_case.ctrl, test_case.target, test_case.shift, n);
    const cmat evolved_rho = internal::kernels::qubit::apply_ctrl_fan_rho(
        rho, A, test_case.ctrl, test_case.target, test_case.shift, n);
    const cmat expected = evolved_psi * evolved_psi.adjoint();

    expect_same(evolved_rho, expected);
    EXPECT_NEAR(evolved_rho.trace().real(), 1.0, tolerance);
    expect_same(evolved_rho, evolved_rho.adjoint());
}
