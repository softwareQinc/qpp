#include "gtest/gtest.h"

#include "qpp/qpp.hpp"

#include "qpp/internal/kernels/qudit/apply_ctrl_diag.hpp"

#include <cmath>
#include <complex>
#include <sstream>
#include <vector>

using namespace qpp;

namespace {

constexpr double tolerance = 1e-11;

struct TestCase {
    idx d;
    idx n;
    std::vector<idx> ctrl;
    std::vector<idx> target;
    std::vector<idx> shift;
};

idx power(idx d, idx n) {
    idx result = 1;
    for (idx k = 0; k < n; ++k) {
        result *= d;
    }
    return result;
}

idx digit(idx basis, idx subsystem, idx d, idx n) {
    return (basis / power(d, n - 1 - subsystem)) % d;
}

// Oracle: for each computational-basis state |x>, the diagonal entry is
// A[target(x)] raised to the collective shifted control value.
// This deliberately does not reconstruct spectator/control base offsets.
cplx diagonal_factor(idx basis, const dyn_col_vect<cplx>& A,
                     const TestCase& tc) {
    idx control_value = 0;
    for (idx j = 0; j < tc.ctrl.size(); ++j) {
        control_value *= tc.d;
        control_value +=
            (digit(basis, tc.ctrl[j], tc.d, tc.n) + tc.shift[j]) % tc.d;
    }

    idx target_value = 0;
    for (idx q : tc.target) {
        target_value = target_value * tc.d + digit(basis, q, tc.d, tc.n);
    }
    return std::pow(A(target_value), control_value);
}

dyn_col_vect<cplx> make_diagonal(const TestCase& tc) {
    const idx dim = power(tc.d, tc.target.size());
    dyn_col_vect<cplx> A(dim);
    for (idx i = 0; i < dim; ++i) {
        A(i) = cplx{0.67 + 0.09 * static_cast<double>(i % 5),
                    -0.22 + 0.05 * static_cast<double>((i * 3) % 7)};
    }
    return A;
}

ket make_state(idx dim) {
    ket result(dim);
    for (idx i = 0; i < dim; ++i) {
        result(i) = cplx{static_cast<double>((i * 7 + 3) % 17) - 8.,
                         static_cast<double>((i * 5 + 1) % 13) - 6.};
    }
    result.normalize();
    return result;
}

cmat make_density(const ket& psi) {
    const idx dim = static_cast<idx>(psi.size());
    const ket other = make_state(dim).conjugate();
    return 0.3 * (psi * psi.adjoint()) +
           0.7 * (other * other.adjoint());
}

ket reference_psi(const ket& original, const dyn_col_vect<cplx>& A,
                  const TestCase& tc) {
    ket expected = original;
    for (idx basis = 0; basis < static_cast<idx>(original.size()); ++basis) {
        expected(basis) *= diagonal_factor(basis, A, tc);
    }
    return expected;
}

cmat reference_rho(const cmat& original, const dyn_col_vect<cplx>& A,
                  const TestCase& tc) {
    cmat expected = original;
    for (idx row = 0; row < static_cast<idx>(original.rows()); ++row) {
        const cplx row_factor = diagonal_factor(row, A, tc);
        for (idx col = 0; col < static_cast<idx>(original.cols()); ++col) {
            expected(row, col) *=
                row_factor * std::conj(diagonal_factor(col, A, tc));
        }
    }
    return expected;
}

template <typename Actual, typename Expected>
void expect_same(const Eigen::MatrixBase<Actual>& actual,
                 const Eigen::MatrixBase<Expected>& expected) {
    ASSERT_EQ(actual.rows(), expected.rows());
    ASSERT_EQ(actual.cols(), expected.cols());
    for (Eigen::Index i = 0; i < actual.rows(); ++i) {
        for (Eigen::Index j = 0; j < actual.cols(); ++j) {
            EXPECT_LE(std::abs(actual(i, j) - expected(i, j)), tolerance)
                << "Mismatch at (" << i << ", " << j << ")";
        }
    }
}

std::vector<TestCase> cases_1q() {
    return {
        {2, 1, {}, {0}, {}},
        {2, 3, {0}, {2}, {0}},
        {2, 3, {2}, {0}, {1}},
        {3, 2, {0}, {1}, {0}},
        {3, 3, {2, 0}, {1}, {1, 2}},
        {3, 4, {3, 1}, {0}, {2, 1}},
        {3, 4, {0, 2, 3}, {1}, {1, 0, 2}},
    };
}

std::vector<TestCase> cases_2q() {
    return {
        {2, 2, {}, {0, 1}, {}},
        {2, 3, {0}, {1, 2}, {0}},
        {2, 4, {3}, {2, 0}, {1}},
        {3, 2, {}, {1, 0}, {}},
        {3, 3, {2}, {1, 0}, {2}},
        {3, 4, {0, 2}, {3, 1}, {1, 2}},
    };
}

std::vector<TestCase> cases_kq() {
    return {
        {2, 2, {0}, {1}, {1}},
        {2, 4, {0}, {3, 1, 2}, {0}},
        {2, 4, {3, 0}, {2, 1}, {1, 0}},
        {3, 3, {2}, {0, 1}, {1}},
        {3, 4, {1}, {3, 0, 2}, {2}},
        {3, 4, {}, {2, 0, 3, 1}, {}},
    };
}

std::string describe(const TestCase& tc) {
    std::ostringstream out;
    out << "d=" << tc.d << " n=" << tc.n << " ctrl={";
    for (idx q : tc.ctrl) out << q << ',';
    out << "} target={";
    for (idx q : tc.target) out << q << ',';
    out << "} shift={";
    for (idx q : tc.shift) out << q << ',';
    out << '}';
    return out.str();
}

} // namespace

TEST(qpp_qudit_apply_ctrl_psi_1q_diag_inplace, ControlPatternsAndShifts) {
    for (const auto& tc : cases_1q()) {
        SCOPED_TRACE(describe(tc));
        const auto A = make_diagonal(tc);
        const ket original = make_state(power(tc.d, tc.n));
        ket actual = original;
        internal::kernels::qudit::apply_ctrl_psi_1q_diag_inplace(
            actual, A, tc.ctrl, tc.target[0], tc.shift, tc.d, tc.n);
        expect_same(actual, reference_psi(original, A, tc));
    }
}

TEST(qpp_qudit_apply_ctrl_psi_2q_diag_inplace, TargetOrderAndShifts) {
    for (const auto& tc : cases_2q()) {
        SCOPED_TRACE(describe(tc));
        const auto A = make_diagonal(tc);
        const ket original = make_state(power(tc.d, tc.n));
        ket actual = original;
        internal::kernels::qudit::apply_ctrl_psi_2q_diag_inplace(
            actual, A, tc.ctrl, tc.target[0], tc.target[1], tc.shift,
            tc.d, tc.n);
        expect_same(actual, reference_psi(original, A, tc));
    }
}

TEST(qpp_qudit_apply_ctrl_psi_kq_diag_inplace, MultipleTargetsAndControls) {
    for (const auto& tc : cases_kq()) {
        SCOPED_TRACE(describe(tc));
        const auto A = make_diagonal(tc);
        const ket original = make_state(power(tc.d, tc.n));
        ket actual = original;
        internal::kernels::qudit::apply_ctrl_psi_kq_diag_inplace(
            actual, A, tc.ctrl, tc.target, tc.shift, tc.d, tc.n);
        expect_same(actual, reference_psi(original, A, tc));
    }
}

TEST(qpp_qudit_apply_ctrl_rho_1q_diag_inplace, ControlPatternsAndShifts) {
    for (const auto& tc : cases_1q()) {
        SCOPED_TRACE(describe(tc));
        const auto A = make_diagonal(tc);
        const cmat original = make_density(make_state(power(tc.d, tc.n)));
        cmat actual = original;
        internal::kernels::qudit::apply_ctrl_rho_1q_diag_inplace(
            actual, A, tc.ctrl, tc.target[0], tc.shift, tc.d, tc.n);
        expect_same(actual, reference_rho(original, A, tc));
    }
}

TEST(qpp_qudit_apply_ctrl_rho_2q_diag_inplace, TargetOrderAndShifts) {
    for (const auto& tc : cases_2q()) {
        SCOPED_TRACE(describe(tc));
        const auto A = make_diagonal(tc);
        const cmat original = make_density(make_state(power(tc.d, tc.n)));
        cmat actual = original;
        internal::kernels::qudit::apply_ctrl_rho_2q_diag_inplace(
            actual, A, tc.ctrl, tc.target[0], tc.target[1], tc.shift,
            tc.d, tc.n);
        expect_same(actual, reference_rho(original, A, tc));
    }
}

TEST(qpp_qudit_apply_ctrl_rho_kq_diag_inplace, MultipleTargetsAndControls) {
    for (const auto& tc : cases_kq()) {
        SCOPED_TRACE(describe(tc));
        const auto A = make_diagonal(tc);
        const cmat original = make_density(make_state(power(tc.d, tc.n)));
        cmat actual = original;
        internal::kernels::qudit::apply_ctrl_rho_kq_diag_inplace(
            actual, A, tc.ctrl, tc.target, tc.shift, tc.d, tc.n);
        expect_same(actual, reference_rho(original, A, tc));
    }
}

TEST(qpp_qudit_apply_ctrl_psi_diag, FunctionalWrappers) {
    for (const auto& tc : cases_1q()) {
        SCOPED_TRACE(describe(tc));
        const auto A = make_diagonal(tc);
        const ket original = make_state(power(tc.d, tc.n));
        const auto actual = internal::kernels::qudit::apply_ctrl_psi_1q_diag(
            original, A, tc.ctrl, tc.target[0], tc.shift, tc.d, tc.n);
        expect_same(actual, reference_psi(original, A, tc));
    }
    for (const auto& tc : cases_2q()) {
        SCOPED_TRACE(describe(tc));
        const auto A = make_diagonal(tc);
        const ket original = make_state(power(tc.d, tc.n));
        const auto actual = internal::kernels::qudit::apply_ctrl_psi_2q_diag(
            original, A, tc.ctrl, tc.target[0], tc.target[1], tc.shift,
            tc.d, tc.n);
        expect_same(actual, reference_psi(original, A, tc));
    }
    for (const auto& tc : cases_kq()) {
        SCOPED_TRACE(describe(tc));
        const auto A = make_diagonal(tc);
        const ket original = make_state(power(tc.d, tc.n));
        const auto actual = internal::kernels::qudit::apply_ctrl_psi_kq_diag(
            original, A, tc.ctrl, tc.target, tc.shift, tc.d, tc.n);
        expect_same(actual, reference_psi(original, A, tc));
    }
}

TEST(qpp_qudit_apply_ctrl_rho_diag, FunctionalWrappers) {
    for (const auto& tc : cases_1q()) {
        SCOPED_TRACE(describe(tc));
        const auto A = make_diagonal(tc);
        const cmat original = make_density(make_state(power(tc.d, tc.n)));
        const auto actual = internal::kernels::qudit::apply_ctrl_rho_1q_diag(
            original, A, tc.ctrl, tc.target[0], tc.shift, tc.d, tc.n);
        expect_same(actual, reference_rho(original, A, tc));
    }
    for (const auto& tc : cases_2q()) {
        SCOPED_TRACE(describe(tc));
        const auto A = make_diagonal(tc);
        const cmat original = make_density(make_state(power(tc.d, tc.n)));
        const auto actual = internal::kernels::qudit::apply_ctrl_rho_2q_diag(
            original, A, tc.ctrl, tc.target[0], tc.target[1], tc.shift,
            tc.d, tc.n);
        expect_same(actual, reference_rho(original, A, tc));
    }
    for (const auto& tc : cases_kq()) {
        SCOPED_TRACE(describe(tc));
        const auto A = make_diagonal(tc);
        const cmat original = make_density(make_state(power(tc.d, tc.n)));
        const auto actual = internal::kernels::qudit::apply_ctrl_rho_kq_diag(
            original, A, tc.ctrl, tc.target, tc.shift, tc.d, tc.n);
        expect_same(actual, reference_rho(original, A, tc));
    }
}

TEST(qpp_qudit_apply_ctrl_diag, PureAndDensityMatrixConsistency) {
    for (const auto& tc : cases_kq()) {
        SCOPED_TRACE(describe(tc));
        const auto A = make_diagonal(tc);
        const ket psi = make_state(power(tc.d, tc.n));
        const cmat rho = psi * psi.adjoint();
        ket psi_actual = psi;
        cmat rho_actual = rho;
        internal::kernels::qudit::apply_ctrl_psi_kq_diag_inplace(
            psi_actual, A, tc.ctrl, tc.target, tc.shift, tc.d, tc.n);
        internal::kernels::qudit::apply_ctrl_rho_kq_diag_inplace(
            rho_actual, A, tc.ctrl, tc.target, tc.shift, tc.d, tc.n);
        expect_same(rho_actual, psi_actual * psi_actual.adjoint());
    }
}
