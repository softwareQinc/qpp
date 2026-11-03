#include "gtest/gtest.h"

#include "qpp/qpp.hpp"

#include "qpp/internal/kernels/qubit/apply_diag.hpp"

#include <cmath>
#include <complex>
#include <sstream>
#include <string>
#include <vector>

using namespace qpp;

namespace {

constexpr double tolerance = 1e-11;

std::string describe(const std::vector<idx>& target) {
    std::ostringstream out;
    out << "targets={";
    for (idx j = 0; j < target.size(); ++j) {
        if (j != 0) {
            out << ",";
        }
        out << target[j];
    }
    out << "}";
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

// Modes: 0 = identity, 1 = unit-modulus phases, 2 = non-unitary complex
// values, 3 = some zero diagonal entries.
ket diagonal_gate(idx d, idx mode) {
    ket A = ket::Ones(d);
    for (idx j = 0; j < d; ++j) {
        if (mode == 1) {
            A(j) = std::polar(1.0, 0.19 * static_cast<double>(j));
        } else if (mode == 2) {
            A(j) = cplx{0.7 + 0.04 * static_cast<double>(j),
                        0.06 * static_cast<double>(j + 1)};
        } else if (mode == 3 && j % 3 == 0) {
            A(j) = 0;
        }
    }
    return A;
}

// Independent big-endian reference: compute the local basis index for each
// global computational-basis state, then look up the diagonal coefficient.
ket diagonal_factors(const ket& A, const std::vector<idx>& target, idx n) {
    const idx D = idx{1} << n;
    ket factors(D);
    for (idx basis = 0; basis < D; ++basis) {
        idx local_basis = 0;
        for (idx q : target) {
            local_basis = (local_basis << 1) |
                          ((basis >> (n - 1 - q)) & idx{1});
        }
        factors(basis) = A(local_basis);
    }
    return factors;
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
                const std::vector<idx>& target, idx n) {
    const ket expected =
        original.cwiseProduct(diagonal_factors(A, target, n));
    expect_same(actual, expected);
}

void expect_rho(const cmat& actual, const cmat& original, const ket& A,
                const std::vector<idx>& target, idx n) {
    const ket factors = diagonal_factors(A, target, n);
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

TEST(qpp_apply_psi_1q_diag_inplace, AllPositionsAndDiagonalShapes) {
    for (idx n : {idx{1}, idx{2}, idx{5}}) {
        const ket original = make_state(n);
        for (idx i = 0; i < n; ++i) {
            SCOPED_TRACE(::testing::Message() << "n=" << n << ", i=" << i);
            for (idx mode : {idx{0}, idx{1}, idx{2}, idx{3}}) {
                const ket A = diagonal_gate(2, mode);
                ket actual = original;
                internal::kernels::qubit::apply_psi_1q_diag_inplace(
                    actual, A, i, n);
                expect_psi(actual, original, A, {i}, n);
            }
        }
    }
}

TEST(qpp_apply_psi_1q_diag_inplace, BigEndianBasisState) {
    ket actual = ket::Zero(8);
    actual(0) = 1;
    actual(4) = 1;
    ket A(2);
    A << cplx{2, 1}, cplx{-1, 3};
    internal::kernels::qubit::apply_psi_1q_diag_inplace(actual, A, 0, 3);

    ket expected = ket::Zero(8);
    expected(0) = A(0);
    expected(4) = A(1);
    expect_same(actual, expected);
}

TEST(qpp_apply_psi_1q_diag_inplace, LargeStateOpenMPThreshold) {
    const idx n = 16;
    const ket original = make_state(n);
    const ket A = diagonal_gate(2, 2);
    for (idx i : {idx{0}, idx{15}}) {
        ket actual = original;
        internal::kernels::qubit::apply_psi_1q_diag_inplace(actual, A, i, n);
        expect_psi(actual, original, A, {i}, n);
    }
}

TEST(qpp_apply_psi_1q_diag_inplace, RowVectorDiagonal) {
    const idx n = 3;
    const ket original = make_state(n);
    const ket A = diagonal_gate(2, 2);
    const Eigen::Matrix<cplx, 1, Eigen::Dynamic> row_A = A.transpose();
    ket actual = original;
    internal::kernels::qubit::apply_psi_1q_diag_inplace(actual, row_A, 1, n);
    expect_psi(actual, original, A, {1}, n);
}

TEST(qpp_apply_psi_1q_diag, ReturnsResultWithoutChangingInput) {
    const idx n = 4;
    ket original = make_state(n);
    const ket before = original;
    const ket A = diagonal_gate(2, 2);
    for (idx i = 0; i < n; ++i) {
        const ket actual = internal::kernels::qubit::apply_psi_1q_diag(
            original, A, i, n);
        expect_psi(actual, before, A, {i}, n);
        expect_same(original, before);
    }
}

TEST(qpp_apply_psi_2q_diag_inplace, TargetOrdersAndSpectators) {
    for (idx n : {idx{2}, idx{5}}) {
        const ket original = make_state(n);
        std::vector<std::vector<idx>> targets;
        if (n == 2) {
            targets = {{0, 1}, {1, 0}};
        } else {
            targets = {{0, 1}, {0, 4}, {4, 0}, {1, 3}, {3, 2}};
        }
        for (const auto& t : targets) {
            SCOPED_TRACE(describe(t));
            for (idx mode : {idx{0}, idx{1}, idx{2}, idx{3}}) {
                const ket A = diagonal_gate(4, mode);
                ket actual = original;
                internal::kernels::qubit::apply_psi_2q_diag_inplace(
                    actual, A, t[0], t[1], n);
                expect_psi(actual, original, A, t, n);
            }
        }
    }
}

TEST(qpp_apply_psi_2q_diag, ReturnsResultWithoutChangingInput) {
    const idx n = 4;
    ket original = make_state(n);
    const ket before = original;
    const ket A = diagonal_gate(4, 2);
    for (const auto& t :
         std::vector<std::vector<idx>>{{0, 3}, {3, 0}, {2, 1}}) {
        const ket actual = internal::kernels::qubit::apply_psi_2q_diag(
            original, A, t[0], t[1], n);
        expect_psi(actual, before, A, t, n);
        expect_same(original, before);
    }
}

TEST(qpp_apply_psi_3q_diag_inplace, PermutationsAndSpectators) {
    for (idx n : {idx{3}, idx{5}}) {
        const ket original = make_state(n);
        std::vector<std::vector<idx>> targets;
        if (n == 3) {
            targets = {{0, 1, 2}, {0, 2, 1}, {1, 0, 2},
                       {1, 2, 0}, {2, 0, 1}, {2, 1, 0}};
        } else {
            targets = {{0, 2, 4}, {4, 0, 2}, {1, 3, 0}, {3, 4, 1}};
        }
        for (const auto& t : targets) {
            SCOPED_TRACE(describe(t));
            for (idx mode : {idx{0}, idx{1}, idx{2}, idx{3}}) {
                const ket A = diagonal_gate(8, mode);
                ket actual = original;
                internal::kernels::qubit::apply_psi_3q_diag_inplace(
                    actual, A, t[0], t[1], t[2], n);
                expect_psi(actual, original, A, t, n);
            }
        }
    }
}

TEST(qpp_apply_psi_3q_diag, ReturnsResultWithoutChangingInput) {
    const idx n = 5;
    ket original = make_state(n);
    const ket before = original;
    const ket A = diagonal_gate(8, 2);
    for (const auto& t :
         std::vector<std::vector<idx>>{{4, 2, 0}, {1, 3, 4}}) {
        const ket actual = internal::kernels::qubit::apply_psi_3q_diag(
            original, A, t[0], t[1], t[2], n);
        expect_psi(actual, before, A, t, n);
        expect_same(original, before);
    }
}

TEST(qpp_apply_psi_kq_diag_inplace, ZeroToFiveTargets) {
    const idx n = 5;
    const ket original = make_state(n);
    const std::vector<std::vector<idx>> targets{
        {}, {0}, {4}, {0, 4}, {4, 0}, {3, 1, 4},
        {4, 2, 0}, {0, 2, 4, 1}, {4, 3, 2, 1, 0}};
    for (const auto& t : targets) {
        SCOPED_TRACE(describe(t));
        const idx d = idx{1} << t.size();
        for (idx mode : {idx{0}, idx{1}, idx{2}, idx{3}}) {
            const ket A = diagonal_gate(d, mode);
            ket actual = original;
            internal::kernels::qubit::apply_psi_kq_diag_inplace(actual, A,
                                                                 t, n);
            expect_psi(actual, original, A, t, n);
        }
    }
}

TEST(qpp_apply_psi_kq_diag_inplace, ZeroQubits) {
    const ket original = make_state(0);
    const ket A = diagonal_gate(1, 2);
    ket actual = original;
    internal::kernels::qubit::apply_psi_kq_diag_inplace(actual, A, {}, 0);
    expect_psi(actual, original, A, {}, 0);
}

TEST(qpp_apply_psi_kq_diag, ReturnsResultWithoutChangingInput) {
    const idx n = 5;
    ket original = make_state(n);
    const ket before = original;
    const std::vector<std::vector<idx>> targets{
        {}, {2}, {3, 0}, {4, 1, 0}, {0, 3, 4, 2}, {2, 0, 4, 1, 3}};
    for (const auto& t : targets) {
        SCOPED_TRACE(describe(t));
        const ket A = diagonal_gate(idx{1} << t.size(), 2);
        const ket actual =
            internal::kernels::qubit::apply_psi_kq_diag(original, A, t, n);
        expect_psi(actual, before, A, t, n);
        expect_same(original, before);
    }
}

TEST(qpp_apply_rho_1q_diag_inplace, MixedStateAllPositions) {
    for (idx n : {idx{1}, idx{4}}) {
        const cmat original = make_density(n);
        for (idx i = 0; i < n; ++i) {
            SCOPED_TRACE(::testing::Message() << "n=" << n << ", i=" << i);
            for (idx mode : {idx{0}, idx{1}, idx{2}, idx{3}}) {
                const ket A = diagonal_gate(2, mode);
                cmat actual = original;
                internal::kernels::qubit::apply_rho_1q_diag_inplace(
                    actual, A, i, n);
                expect_rho(actual, original, A, {i}, n);
            }
        }
    }
}

TEST(qpp_apply_rho_1q_diag, ReturnsResultWithoutChangingInput) {
    const idx n = 4;
    cmat original = make_density(n);
    const cmat before = original;
    const ket A = diagonal_gate(2, 2);
    for (idx i = 0; i < n; ++i) {
        const cmat actual = internal::kernels::qubit::apply_rho_1q_diag(
            original, A, i, n);
        expect_rho(actual, before, A, {i}, n);
        expect_same(original, before);
    }
}

TEST(qpp_apply_rho_2q_diag_inplace, MixedStateTargetOrders) {
    for (idx n : {idx{2}, idx{5}}) {
        const cmat original = make_density(n);
        std::vector<std::vector<idx>> targets;
        if (n == 2) {
            targets = {{0, 1}, {1, 0}};
        } else {
            targets = {{0, 4}, {4, 0}, {1, 3}, {3, 2}};
        }
        for (const auto& t : targets) {
            SCOPED_TRACE(describe(t));
            for (idx mode : {idx{0}, idx{1}, idx{2}, idx{3}}) {
                const ket A = diagonal_gate(4, mode);
                cmat actual = original;
                internal::kernels::qubit::apply_rho_2q_diag_inplace(
                    actual, A, t[0], t[1], n);
                expect_rho(actual, original, A, t, n);
            }
        }
    }
}

TEST(qpp_apply_rho_2q_diag, ReturnsResultWithoutChangingInput) {
    const idx n = 4;
    cmat original = make_density(n);
    const cmat before = original;
    const ket A = diagonal_gate(4, 2);
    for (const auto& t :
         std::vector<std::vector<idx>>{{0, 3}, {3, 0}, {2, 1}}) {
        const cmat actual = internal::kernels::qubit::apply_rho_2q_diag(
            original, A, t[0], t[1], n);
        expect_rho(actual, before, A, t, n);
        expect_same(original, before);
    }
}

TEST(qpp_apply_rho_3q_diag_inplace, MixedStatePermutations) {
    for (idx n : {idx{3}, idx{5}}) {
        const cmat original = make_density(n);
        std::vector<std::vector<idx>> targets;
        if (n == 3) {
            targets = {{0, 1, 2}, {2, 1, 0}, {1, 0, 2}};
        } else {
            targets = {{0, 2, 4}, {4, 0, 2}, {1, 3, 0}, {3, 4, 1}};
        }
        for (const auto& t : targets) {
            SCOPED_TRACE(describe(t));
            for (idx mode : {idx{0}, idx{1}, idx{2}, idx{3}}) {
                const ket A = diagonal_gate(8, mode);
                cmat actual = original;
                internal::kernels::qubit::apply_rho_3q_diag_inplace(
                    actual, A, t[0], t[1], t[2], n);
                expect_rho(actual, original, A, t, n);
            }
        }
    }
}

TEST(qpp_apply_rho_3q_diag, ReturnsResultWithoutChangingInput) {
    const idx n = 5;
    cmat original = make_density(n);
    const cmat before = original;
    const ket A = diagonal_gate(8, 2);
    for (const auto& t :
         std::vector<std::vector<idx>>{{4, 2, 0}, {1, 3, 4}}) {
        const cmat actual = internal::kernels::qubit::apply_rho_3q_diag(
            original, A, t[0], t[1], t[2], n);
        expect_rho(actual, before, A, t, n);
        expect_same(original, before);
    }
}

TEST(qpp_apply_rho_kq_diag_inplace, ZeroToFiveTargets) {
    const idx n = 5;
    const cmat original = make_density(n);
    const std::vector<std::vector<idx>> targets{
        {}, {0}, {4}, {0, 4}, {4, 0}, {3, 1, 4},
        {4, 2, 0}, {0, 2, 4, 1}, {4, 3, 2, 1, 0}};
    for (const auto& t : targets) {
        SCOPED_TRACE(describe(t));
        const idx d = idx{1} << t.size();
        for (idx mode : {idx{0}, idx{1}, idx{2}, idx{3}}) {
            const ket A = diagonal_gate(d, mode);
            cmat actual = original;
            internal::kernels::qubit::apply_rho_kq_diag_inplace(actual, A,
                                                                 t, n);
            expect_rho(actual, original, A, t, n);
        }
    }
}

TEST(qpp_apply_rho_kq_diag_inplace, RowVectorDiagonal) {
    const idx n = 4;
    const cmat original = make_density(n);
    const std::vector<idx> targets{3, 0, 2};
    const ket A = diagonal_gate(8, 2);
    const Eigen::Matrix<cplx, 1, Eigen::Dynamic> row_A = A.transpose();
    cmat actual = original;
    internal::kernels::qubit::apply_rho_kq_diag_inplace(
        actual, row_A, targets, n);
    expect_rho(actual, original, A, targets, n);
}

TEST(qpp_apply_rho_kq_diag, ReturnsResultWithoutChangingInput) {
    const idx n = 5;
    cmat original = make_density(n);
    const cmat before = original;
    const std::vector<std::vector<idx>> targets{
        {}, {2}, {3, 0}, {4, 1, 0}, {0, 3, 4, 2}, {2, 0, 4, 1, 3}};
    for (const auto& t : targets) {
        SCOPED_TRACE(describe(t));
        const ket A = diagonal_gate(idx{1} << t.size(), 2);
        const cmat actual =
            internal::kernels::qubit::apply_rho_kq_diag(original, A, t, n);
        expect_rho(actual, before, A, t, n);
        expect_same(original, before);
    }
}

TEST(qpp_apply_rho_kq_diag, PureStateConsistencyAndUnitaryTrace) {
    const idx n = 5;
    const ket psi = make_state(n);
    const cmat rho = psi * psi.adjoint();
    const std::vector<idx> targets{4, 0, 2};
    const ket A = diagonal_gate(8, 1);

    const ket transformed_psi =
        internal::kernels::qubit::apply_psi_kq_diag(psi, A, targets, n);
    const cmat transformed_rho =
        internal::kernels::qubit::apply_rho_kq_diag(rho, A, targets, n);
    const cmat expected = transformed_psi * transformed_psi.adjoint();

    expect_same(transformed_rho, expected);
    EXPECT_NEAR(transformed_rho.trace().real(), 1.0, tolerance);
    expect_same(transformed_rho, transformed_rho.adjoint());
}
