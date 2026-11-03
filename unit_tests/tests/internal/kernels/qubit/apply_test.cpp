#include "gtest/gtest.h"

#include "qpp/qpp.hpp"

#include "qpp/internal/kernels/qubit/apply.hpp"

#include <cmath>
#include <complex>
#include <vector>

using namespace qpp;

namespace {

constexpr double tolerance = 1e-11;

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

// A dense unitary with nontrivial complex phases in every row and column.
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

// The kernels also accept arbitrary, not necessarily unitary, matrices.
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

// Independent reference: explicitly embed the local operator in the entire
// 2^n-dimensional space. The order of 'target' defines local big-endian order.
cmat embedded_gate(const cmat& A, const std::vector<idx>& target, idx n) {
    const idx D = idx{1} << n;
    const idx k = target.size();
    cmat result = cmat::Zero(D, D);

    for (idx input = 0; input < D; ++input) {
        idx local_input = 0;
        for (idx q : target) {
            local_input =
                (local_input << 1) | ((input >> (n - 1 - q)) & idx{1});
        }

        for (idx local_output = 0; local_output < (idx{1} << k);
             ++local_output) {
            idx output = input;
            for (idx j = 0; j < k; ++j) {
                const idx mask = idx{1} << (n - 1 - target[j]);
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
                const std::vector<idx>& target, idx n) {
    const ket expected = embedded_gate(gate, target, n) * original;
    expect_same(actual, expected);
}

void expect_rho(const cmat& actual, const cmat& original, const cmat& gate,
                const std::vector<idx>& target, idx n) {
    const cmat U = embedded_gate(gate, target, n);
    const cmat expected = U * original * U.adjoint();
    expect_same(actual, expected);
}

} // namespace

TEST(qpp_apply_psi_1q_inplace, HadamardFastPath) {
    for (idx n : {idx{1}, idx{2}, idx{5}}) {
        const ket original = make_state(n);
        for (idx i = 0; i < n; ++i) {
            SCOPED_TRACE(::testing::Message() << "n=" << n << ", i=" << i);
            ket actual = original;
            internal::kernels::qubit::apply_psi_1q_inplace(actual, gt.H, i, n);
            expect_psi(actual, original, gt.H, {i}, n);
        }
    }
}

TEST(qpp_apply_psi_1q_inplace, ComplexAndNonunitaryGates) {
    for (idx n : {idx{1}, idx{4}}) {
        const ket original = make_state(n);
        for (const cmat& gate : {fourier_gate(2), dense_gate(2)}) {
            for (idx i = 0; i < n; ++i) {
                SCOPED_TRACE(::testing::Message() << "n=" << n << ", i=" << i);
                ket actual = original;
                internal::kernels::qubit::apply_psi_1q_inplace(actual, gate, i,
                                                               n);
                expect_psi(actual, original, gate, {i}, n);
            }
        }
    }
}

TEST(qpp_apply_psi_1q_inplace, BigEndianBasisState) {
    ket state = ket::Zero(8);
    state(0) = 1;
    internal::kernels::qubit::apply_psi_1q_inplace(state, gt.H, 0, 3);

    ket expected = ket::Zero(8);
    expected(0) = 1.0 / std::sqrt(2.0);
    expected(4) = 1.0 / std::sqrt(2.0);
    expect_same(state, expected);
}

TEST(qpp_apply_psi_1q_inplace, LargeStateBothPaths) {
    // Exercises the OpenMP threshold (2^16 amplitudes) when OpenMP is enabled.
    const idx n = 16;
    const ket original = make_state(n);
    for (const cmat& gate : {cmat(gt.H), fourier_gate(2)}) {
        for (idx i : {idx{0}, idx{15}}) {
            SCOPED_TRACE(::testing::Message() << "i=" << i);
            ket expected = original;
            const idx mask = idx{1} << (n - 1 - i);
            for (idx p = 0; p < static_cast<idx>(original.size()); ++p) {
                if ((p & mask) != 0) {
                    continue;
                }
                const cplx a = original(p);
                const cplx b = original(p | mask);
                expected(p) = gate(0, 0) * a + gate(0, 1) * b;
                expected(p | mask) = gate(1, 0) * a + gate(1, 1) * b;
            }
            ket actual = original;
            internal::kernels::qubit::apply_psi_1q_inplace(actual, gate, i, n);
            expect_same(actual, expected);
        }
    }
}

TEST(qpp_apply_psi_1q, ReturnsResultWithoutChangingInput) {
    const idx n = 4;
    ket original = make_state(n);
    const ket before = original;
    for (idx i = 0; i < n; ++i) {
        const cmat gate = fourier_gate(2);
        const ket actual =
            internal::kernels::qubit::apply_psi_1q(original, gate, i, n);
        expect_psi(actual, before, gate, {i}, n);
        expect_same(original, before);
    }
}

TEST(qpp_apply_psi_2q_inplace, AllTargetOrdersAndSpectators) {
    for (idx n : {idx{2}, idx{5}}) {
        const ket original = make_state(n);
        std::vector<std::vector<idx>> targets;
        if (n == 2) {
            targets = {{0, 1}, {1, 0}};
        } else {
            targets = {{0, 1}, {0, 4}, {4, 0}, {1, 3}, {3, 2}};
        }
        for (const auto& t : targets) {
            SCOPED_TRACE(::testing::Message()
                         << "n=" << n << ", targets=" << t[0] << "," << t[1]);
            for (const cmat& gate : {fourier_gate(4), dense_gate(4)}) {
                ket actual = original;
                internal::kernels::qubit::apply_psi_2q_inplace(actual, gate,
                                                               t[0], t[1], n);
                expect_psi(actual, original, gate, t, n);
            }
        }
    }
}

TEST(qpp_apply_psi_2q, ReturnsResultWithoutChangingInput) {
    const idx n = 4;
    ket original = make_state(n);
    const ket before = original;
    const cmat gate = fourier_gate(4);
    for (const auto& t :
         std::vector<std::vector<idx>>{{0, 3}, {3, 0}, {2, 1}}) {
        const ket actual = internal::kernels::qubit::apply_psi_2q(
            original, gate, t[0], t[1], n);
        expect_psi(actual, before, gate, t, n);
        expect_same(original, before);
    }
}

TEST(qpp_apply_psi_3q_inplace, PermutationsAndSpectators) {
    for (idx n : {idx{3}, idx{5}}) {
        const ket original = make_state(n);
        std::vector<std::vector<idx>> targets;
        if (n == 3) {
            targets = {{0, 1, 2}, {2, 0, 1}, {2, 1, 0}};
        } else {
            targets = {{0, 2, 4}, {4, 0, 2}, {1, 3, 0}, {3, 4, 1}};
        }
        for (const auto& t : targets) {
            SCOPED_TRACE(::testing::Message()
                         << "n=" << n << ", first=" << t[0]);
            for (const cmat& gate : {fourier_gate(8), dense_gate(8)}) {
                ket actual = original;
                internal::kernels::qubit::apply_psi_3q_inplace(
                    actual, gate, t[0], t[1], t[2], n);
                expect_psi(actual, original, gate, t, n);
            }
        }
    }
}

TEST(qpp_apply_psi_3q, ReturnsResultWithoutChangingInput) {
    const idx n = 5;
    ket original = make_state(n);
    const ket before = original;
    const cmat gate = fourier_gate(8);
    for (const auto& t : std::vector<std::vector<idx>>{{4, 2, 0}, {1, 3, 4}}) {
        const ket actual = internal::kernels::qubit::apply_psi_3q(
            original, gate, t[0], t[1], t[2], n);
        expect_psi(actual, before, gate, t, n);
        expect_same(original, before);
    }
}

TEST(qpp_apply_psi_kq_inplace, ZeroToFourTargets) {
    const idx n = 5;
    const ket original = make_state(n);
    const std::vector<std::vector<idx>> targets{
        {},        {0},       {4},          {0, 4},      {4, 0},
        {3, 1, 4}, {4, 2, 0}, {0, 2, 4, 1}, {4, 3, 2, 1}};
    for (const auto& t : targets) {
        SCOPED_TRACE(::testing::Message() << "target count=" << t.size());
        const idx d = idx{1} << t.size();
        for (const cmat& gate : {fourier_gate(d), dense_gate(d)}) {
            ket actual = original;
            internal::kernels::qubit::apply_psi_kq_inplace(actual, gate, t, n);
            expect_psi(actual, original, gate, t, n);
        }
    }
}

TEST(qpp_apply_psi_kq_inplace, AllQubitsInArbitraryOrder) {
    const idx n = 4;
    const ket original = make_state(n);
    const std::vector<idx> targets{2, 0, 3, 1};
    const cmat gate = fourier_gate(16);
    ket actual = original;
    internal::kernels::qubit::apply_psi_kq_inplace(actual, gate, targets, n);
    expect_psi(actual, original, gate, targets, n);
}

TEST(qpp_apply_psi_kq, ReturnsResultWithoutChangingInput) {
    const idx n = 5;
    ket original = make_state(n);
    const ket before = original;
    for (const auto& t : std::vector<std::vector<idx>>{
             {}, {2}, {3, 0}, {4, 1, 0}, {0, 3, 4, 2}}) {
        const cmat gate = fourier_gate(idx{1} << t.size());
        const ket actual =
            internal::kernels::qubit::apply_psi_kq(original, gate, t, n);
        expect_psi(actual, before, gate, t, n);
        expect_same(original, before);
    }
}

TEST(qpp_apply_rho_1q_inplace, MixedStateAllPositions) {
    for (idx n : {idx{1}, idx{4}}) {
        const cmat original = make_density(n);
        for (idx i = 0; i < n; ++i) {
            SCOPED_TRACE(::testing::Message() << "n=" << n << ", i=" << i);
            for (const cmat& gate :
                 {cmat(gt.H), fourier_gate(2), dense_gate(2)}) {
                cmat actual = original;
                internal::kernels::qubit::apply_rho_1q_inplace(actual, gate, i,
                                                               n);
                expect_rho(actual, original, gate, {i}, n);
            }
        }
    }
}

TEST(qpp_apply_rho_1q, ReturnsResultWithoutChangingInput) {
    const idx n = 3;
    cmat original = make_density(n);
    const cmat before = original;
    const cmat gate = fourier_gate(2);
    for (idx i = 0; i < n; ++i) {
        const cmat actual =
            internal::kernels::qubit::apply_rho_1q(original, gate, i, n);
        expect_rho(actual, before, gate, {i}, n);
        expect_same(original, before);
    }
}

TEST(qpp_apply_rho_2q_inplace, BothOrdersAndSpectators) {
    for (idx n : {idx{2}, idx{4}}) {
        const cmat original = make_density(n);
        std::vector<std::vector<idx>> targets;
        if (n == 2) {
            targets = {{0, 1}, {1, 0}};
        } else {
            targets = {{0, 3}, {3, 0}, {1, 2}, {2, 0}};
        }
        for (const auto& t : targets) {
            SCOPED_TRACE(::testing::Message()
                         << "n=" << n << ", first=" << t[0]);
            for (const cmat& gate : {fourier_gate(4), dense_gate(4)}) {
                cmat actual = original;
                internal::kernels::qubit::apply_rho_2q_inplace(actual, gate,
                                                               t[0], t[1], n);
                expect_rho(actual, original, gate, t, n);
            }
        }
    }
}

TEST(qpp_apply_rho_2q, ReturnsResultWithoutChangingInput) {
    const idx n = 4;
    cmat original = make_density(n);
    const cmat before = original;
    const cmat gate = fourier_gate(4);
    for (const auto& t : std::vector<std::vector<idx>>{{0, 3}, {3, 1}}) {
        const cmat actual = internal::kernels::qubit::apply_rho_2q(
            original, gate, t[0], t[1], n);
        expect_rho(actual, before, gate, t, n);
        expect_same(original, before);
    }
}

TEST(qpp_apply_rho_3q_inplace, PermutationsAndSpectators) {
    for (idx n : {idx{3}, idx{4}}) {
        const cmat original = make_density(n);
        std::vector<std::vector<idx>> targets;
        if (n == 3) {
            targets = {{0, 1, 2}, {2, 1, 0}};
        } else {
            targets = {{0, 2, 3}, {3, 0, 1}, {2, 1, 3}};
        }
        for (const auto& t : targets) {
            SCOPED_TRACE(::testing::Message()
                         << "n=" << n << ", first=" << t[0]);
            for (const cmat& gate : {fourier_gate(8), dense_gate(8)}) {
                cmat actual = original;
                internal::kernels::qubit::apply_rho_3q_inplace(
                    actual, gate, t[0], t[1], t[2], n);
                expect_rho(actual, original, gate, t, n);
            }
        }
    }
}

TEST(qpp_apply_rho_3q, ReturnsResultWithoutChangingInput) {
    const idx n = 4;
    cmat original = make_density(n);
    const cmat before = original;
    const cmat gate = fourier_gate(8);
    for (const auto& t : std::vector<std::vector<idx>>{{3, 2, 0}, {1, 0, 2}}) {
        const cmat actual = internal::kernels::qubit::apply_rho_3q(
            original, gate, t[0], t[1], t[2], n);
        expect_rho(actual, before, gate, t, n);
        expect_same(original, before);
    }
}

TEST(qpp_apply_rho_kq_inplace, ZeroToFourTargets) {
    const idx n = 4;
    const cmat original = make_density(n);
    const std::vector<std::vector<idx>> targets{
        {}, {0}, {3}, {0, 3}, {3, 0}, {2, 0, 3}, {3, 1, 0}, {2, 0, 3, 1}};
    for (const auto& t : targets) {
        SCOPED_TRACE(::testing::Message() << "target count=" << t.size());
        const idx d = idx{1} << t.size();
        for (const cmat& gate : {fourier_gate(d), dense_gate(d)}) {
            cmat actual = original;
            internal::kernels::qubit::apply_rho_kq_inplace(actual, gate, t, n);
            expect_rho(actual, original, gate, t, n);
        }
    }
}

TEST(qpp_apply_rho_kq, ReturnsResultWithoutChangingInput) {
    const idx n = 4;
    cmat original = make_density(n);
    const cmat before = original;
    for (const auto& t : std::vector<std::vector<idx>>{
             {}, {1}, {3, 0}, {2, 0, 1}, {3, 1, 0, 2}}) {
        const cmat gate = fourier_gate(idx{1} << t.size());
        const cmat actual =
            internal::kernels::qubit::apply_rho_kq(original, gate, t, n);
        expect_rho(actual, before, gate, t, n);
        expect_same(original, before);
    }
}

TEST(qpp_apply_rho_kq, PureStateConsistency) {
    const idx n = 4;
    const ket psi = make_state(n);
    const cmat rho = psi * psi.adjoint();
    const std::vector<idx> target{3, 0, 2};
    const cmat gate = fourier_gate(8);

    const ket transformed_psi =
        internal::kernels::qubit::apply_psi_kq(psi, gate, target, n);
    const cmat transformed_rho =
        internal::kernels::qubit::apply_rho_kq(rho, gate, target, n);
    const cmat expected = transformed_psi * transformed_psi.adjoint();

    expect_same(transformed_rho, expected);
    EXPECT_NEAR(transformed_rho.trace().real(), 1.0, tolerance);
    expect_same(transformed_rho, transformed_rho.adjoint());
}
