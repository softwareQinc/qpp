#include "gtest/gtest.h"

#include "qpp/qpp.hpp"

#include "qpp/internal/kernels/qubit/syspermute.hpp"

#include <complex>
#include <sstream>
#include <string>
#include <vector>

using namespace qpp;

namespace {

constexpr double tolerance = 1e-11;

std::string describe(const std::vector<idx>& perm) {
    std::ostringstream out;
    out << "perm={";
    for (idx j = 0; j < perm.size(); ++j) {
        if (j != 0) {
            out << ",";
        }
        out << perm[j];
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

cmat make_general_matrix(idx n) {
    const idx D = idx{1} << n;
    cmat result(D, D);
    for (idx r = 0; r < D; ++r) {
        for (idx c = 0; c < D; ++c) {
            result(r, c) =
                cplx{0.07 * (r + 1) - 0.04 * (c + 1),
                     0.02 * (2 * r + c + 1)};
        }
    }
    return result;
}

// Independent reference: the new qubit at position j contains the old qubit
// at position perm[j], with qubit 0 being the most-significant bit.
idx permuted_index(idx old_index, const std::vector<idx>& perm, idx n) {
    idx new_index = 0;
    for (idx new_qubit = 0; new_qubit < n; ++new_qubit) {
        const idx old_qubit = perm[new_qubit];
        const idx old_bit =
            (old_index >> (n - 1 - old_qubit)) & idx{1};
        new_index |= old_bit << (n - 1 - new_qubit);
    }
    return new_index;
}

ket reference_psi(const ket& psi, const std::vector<idx>& perm, idx n) {
    const idx D = idx{1} << n;
    ket result = ket::Zero(D);
    for (idx old_index = 0; old_index < D; ++old_index) {
        result(permuted_index(old_index, perm, n)) = psi(old_index);
    }
    return result;
}

cmat reference_rho(const cmat& rho, const std::vector<idx>& perm, idx n) {
    const idx D = idx{1} << n;
    cmat result = cmat::Zero(D, D);
    for (idx old_row = 0; old_row < D; ++old_row) {
        for (idx old_col = 0; old_col < D; ++old_col) {
            result(permuted_index(old_row, perm, n),
                   permuted_index(old_col, perm, n)) =
                rho(old_row, old_col);
        }
    }
    return result;
}

std::vector<idx> inverse_permutation(const std::vector<idx>& perm) {
    std::vector<idx> inverse(perm.size());
    for (idx j = 0; j < perm.size(); ++j) {
        inverse[perm[j]] = j;
    }
    return inverse;
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

} // namespace

TEST(qpp_qubit_syspermute_psi_kq, IdentityPermutation) {
    for (idx n : {idx{1}, idx{2}, idx{4}}) {
        const ket original = make_state(n);
        std::vector<idx> perm;
        for (idx q = 0; q < n; ++q) {
            perm.push_back(q);
        }
        const ket actual =
            internal::kernels::qubit::syspermute_psi_kq(original, perm, n);
        expect_same(actual, original);
    }
}

TEST(qpp_qubit_syspermute_psi_kq, TwoQubitSwap) {
    ket original(4);
    original << 1, 2, 3, 4;
    ket expected(4);
    expected << 1, 3, 2, 4;
    const ket actual =
        internal::kernels::qubit::syspermute_psi_kq(original, {1, 0}, 2);
    expect_same(actual, expected);
}

TEST(qpp_qubit_syspermute_psi_kq, ThreeQubitAdjacentSwap) {
    // perm={0,2,1} maps |001> to |010> in big-endian order.
    ket original = ket::Zero(8);
    original(1) = 1;
    const ket actual =
        internal::kernels::qubit::syspermute_psi_kq(original, {0, 2, 1}, 3);
    ket expected = ket::Zero(8);
    expected(2) = 1;
    expect_same(actual, expected);
}

TEST(qpp_qubit_syspermute_psi_kq, ArbitraryThreeAndFourQubitPermutations) {
    for (idx n : {idx{3}, idx{4}}) {
        const ket original = make_state(n);
        std::vector<std::vector<idx>> permutations;
        if (n == 3) {
            permutations = {{0, 1, 2}, {0, 2, 1}, {1, 0, 2},
                            {1, 2, 0}, {2, 0, 1}, {2, 1, 0}};
        } else {
            permutations = {{0, 1, 2, 3}, {1, 0, 2, 3}, {0, 2, 3, 1},
                            {3, 0, 2, 1}, {2, 3, 1, 0}, {3, 2, 1, 0}};
        }
        for (const auto& perm : permutations) {
            SCOPED_TRACE(describe(perm));
            const ket actual =
                internal::kernels::qubit::syspermute_psi_kq(original, perm,
                                                              n);
            expect_same(actual, reference_psi(original, perm, n));
            EXPECT_NEAR(actual.norm(), original.norm(), tolerance);
        }
    }
}

TEST(qpp_qubit_syspermute_psi_kq, PermutationAndInverseRestoreState) {
    const idx n = 5;
    const ket original = make_state(n);
    const std::vector<idx> perm{3, 0, 4, 1, 2};
    const ket once =
        internal::kernels::qubit::syspermute_psi_kq(original, perm, n);
    const ket twice = internal::kernels::qubit::syspermute_psi_kq(
        once, inverse_permutation(perm), n);
    expect_same(twice, original);
    expect_same(once, reference_psi(original, perm, n));
}

TEST(qpp_qubit_syspermute_psi_kq, InputIsNotModified) {
    const idx n = 4;
    ket original = make_state(n);
    const ket before = original;
    const std::vector<idx> perm{2, 0, 3, 1};
    const ket actual =
        internal::kernels::qubit::syspermute_psi_kq(original, perm, n);
    expect_same(actual, reference_psi(before, perm, n));
    expect_same(original, before);
}

TEST(qpp_qubit_syspermute_psi_kq, ZeroQubitSystem) {
    ket scalar(1);
    scalar(0) = cplx{1.5, -0.5};
    const ket actual =
        internal::kernels::qubit::syspermute_psi_kq(scalar, {}, 0);
    expect_same(actual, scalar);
}

TEST(qpp_qubit_syspermute_rho_kq, IdentityPermutation) {
    for (idx n : {idx{1}, idx{2}, idx{4}}) {
        const cmat original = make_general_matrix(n);
        std::vector<idx> perm;
        for (idx q = 0; q < n; ++q) {
            perm.push_back(q);
        }
        const cmat actual =
            internal::kernels::qubit::syspermute_rho_kq(original, perm, n);
        expect_same(actual, original);
    }
}

TEST(qpp_qubit_syspermute_rho_kq, TwoQubitSwapKnownMatrix) {
    cmat original(4, 4);
    original << 1, 2, 3, 4,
                5, 6, 7, 8,
                9, 10, 11, 12,
                13, 14, 15, 16;
    cmat expected(4, 4);
    expected << 1, 3, 2, 4,
                9, 11, 10, 12,
                5, 7, 6, 8,
                13, 15, 14, 16;
    const cmat actual =
        internal::kernels::qubit::syspermute_rho_kq(original, {1, 0}, 2);
    expect_same(actual, expected);
}

TEST(qpp_qubit_syspermute_rho_kq, ThreeQubitAdjacentSwap) {
    cmat original = cmat::Zero(8, 8);
    original(1, 6) = cplx{0.7, 0.2};
    const cmat actual = internal::kernels::qubit::syspermute_rho_kq(
        original, {0, 2, 1}, 3);
    cmat expected = cmat::Zero(8, 8);
    // |001><110| becomes |010><101|.
    expected(2, 5) = cplx{0.7, 0.2};
    expect_same(actual, expected);
}

TEST(qpp_qubit_syspermute_rho_kq, MixedStateArbitraryPermutations) {
    const idx n = 5;
    const cmat original = make_density(n);
    const cmat before = original;
    const std::vector<std::vector<idx>> permutations{
        {0, 1, 2, 3, 4}, {0, 2, 1, 3, 4}, {4, 1, 2, 3, 0},
        {3, 0, 4, 1, 2}, {2, 4, 1, 0, 3}};
    for (const auto& perm : permutations) {
        SCOPED_TRACE(describe(perm));
        const cmat actual =
            internal::kernels::qubit::syspermute_rho_kq(original, perm, n);
        expect_same(actual, reference_rho(original, perm, n));
        EXPECT_NEAR(actual.trace().real(), 1.0, tolerance);
        expect_same(actual, actual.adjoint());
        expect_same(original, before);
    }
}

TEST(qpp_qubit_syspermute_rho_kq, GeneralNonHermitianMatrix) {
    const idx n = 4;
    const cmat original = make_general_matrix(n);
    const std::vector<std::vector<idx>> permutations{
        {0, 2, 1, 3}, {1, 0, 2, 3}, {2, 3, 0, 1}, {3, 0, 2, 1}};
    for (const auto& perm : permutations) {
        SCOPED_TRACE(describe(perm));
        const cmat actual =
            internal::kernels::qubit::syspermute_rho_kq(original, perm, n);
        expect_same(actual, reference_rho(original, perm, n));
    }
}

TEST(qpp_qubit_syspermute_rho_kq, RowMajorMatrixInput) {
    const idx n = 3;
    const cmat original = make_general_matrix(n);
    const Eigen::Matrix<cplx, Eigen::Dynamic, Eigen::Dynamic,
                        Eigen::RowMajor>
        row_major = original;
    const std::vector<idx> perm{0, 2, 1};

    const cmat actual =
        internal::kernels::qubit::syspermute_rho_kq(row_major, perm, n);
    expect_same(actual, reference_rho(original, perm, n));
}

TEST(qpp_qubit_syspermute_rho_kq, NoncontiguousBlockExpression) {
    const idx n = 2;
    const cmat original = make_general_matrix(n);
    cmat padded = cmat::Zero(6, 6);
    padded.block(1, 1, 4, 4) = original;
    const std::vector<idx> perm{1, 0};

    const cmat actual = internal::kernels::qubit::syspermute_rho_kq(
        padded.block(1, 1, 4, 4), perm, n);
    expect_same(actual, reference_rho(original, perm, n));
}

TEST(qpp_qubit_syspermute_rho_kq, PermutationAndInverseRestoreMatrix) {
    const idx n = 4;
    const cmat original = make_general_matrix(n);
    const std::vector<idx> perm{2, 0, 3, 1};
    const cmat once =
        internal::kernels::qubit::syspermute_rho_kq(original, perm, n);
    const cmat twice = internal::kernels::qubit::syspermute_rho_kq(
        once, inverse_permutation(perm), n);
    expect_same(twice, original);
    expect_same(once, reference_rho(original, perm, n));
}

TEST(qpp_qubit_syspermute_rho_kq, PureStateConsistency) {
    const idx n = 5;
    const ket psi = make_state(n);
    const cmat rho = psi * psi.adjoint();
    const std::vector<idx> perm{3, 0, 4, 1, 2};
    const ket permuted_psi =
        internal::kernels::qubit::syspermute_psi_kq(psi, perm, n);
    const cmat permuted_rho =
        internal::kernels::qubit::syspermute_rho_kq(rho, perm, n);
    const cmat expected = permuted_psi * permuted_psi.adjoint();
    expect_same(permuted_rho, expected);
    expect_same(permuted_rho, reference_rho(rho, perm, n));
}

TEST(qpp_qubit_syspermute_rho_kq, ZeroQubitSystem) {
    cmat scalar(1, 1);
    scalar(0, 0) = cplx{0.7, -0.2};
    const cmat actual =
        internal::kernels::qubit::syspermute_rho_kq(scalar, {}, 0);
    expect_same(actual, scalar);
}
