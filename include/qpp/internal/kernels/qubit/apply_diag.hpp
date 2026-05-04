/*
 * This file is part of Quantum++.
 *
 * Copyright (c) 2017 - 2026 softwareQ Inc. All rights reserved.
 *
 * MIT License
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

/**
 * @file qpp/internal/kernels/qubit/apply_diag.hpp
 * @brief Internal highly optimized critical functions for qpp::apply_diag()
 */

#ifndef QPP_INTERNAL_KERNELS_QUBIT_APPLY_DIAG_HPP_
#define QPP_INTERNAL_KERNELS_QUBIT_APPLY_DIAG_HPP_

#include <cassert>
#include <vector>
#ifndef NDEBUG
#include <set>
#endif

#include <Eigen/Dense>

#include "qpp/internal/util.hpp"
#include "qpp/types.hpp"

namespace qpp::internal::kernels::qubit {
/**
 * @brief Applies the 1-qubit diagonal gate \a A to the qubit \a i of the
 * multi-partite state vector \a state in-place
 *
 * @param state Eigen expression (modified in-place)
 * @param A Eigen expression (2 x 1 vector of diagonal elements)
 * @param i Subsystem index where the gate \a A is applied
 * @param n Number of qubits
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] void
apply_psi_1q_diag_inplace(Eigen::MatrixBase<Derived1>& state,
                          const Eigen::MatrixBase<Derived2>& A, idx i, idx n) {
    // Type and Dimension Setup
    using Scalar = typename Derived1::Scalar;
    const idx D = static_cast<idx>(std::size_t{1} << n); // Total size 2^n

    // Input Validation
    assert(i < n && "Target qubit index i must be less than n");
    assert(static_cast<idx>(state.size()) == D &&
           "State vector size must be 2^n");
    // A is a diagonal gate, typically passed as a 2x1 or 1x2 vector
    // of its diagonal entries
    assert(A.size() == 2 && "Diagonal gate A must have 2 elements");

    // Pre-calculate Bit-Shift Constants (BIG-ENDIAN)
    const idx j = n - 1 - i;
    const idx step = static_cast<idx>(std::size_t{1} << j);
    const idx jump = static_cast<idx>(std::size_t{1} << (j + 1));

    // Extract the diagonal elements
    const Scalar a00 = A.coeff(0);
    const Scalar a11 = A.coeff(1);

    // Pair-wise Amplitude Transformation
#ifdef QPP_OPENMP
// NOLINTNEXTLINE
#pragma omp parallel for
#endif // QPP_OPENMP
    for (idx L = 0; L < D; L += jump) {
        for (idx R = 0; R < step; ++R) {
            const idx k0 = L + R;
            const idx k1 = k0 + step;

            // Since A is diagonal, we simply scale the existing amplitudes
            // psi'_k0 = A[0] * psi_k0
            // psi'_k1 = A[1] * psi_k1
            state.coeffRef(k0) *= a00;
            state.coeffRef(k1) *= a11;
        }
    }
}

/**
 * @brief Applies the 1-qubit diagonal gate \a A to the qubit \a i of the
 * multi-partite state vector \a state
 *
 * @param state Eigen expression
 * @param A Eigen expression (2 x 1 vector of diagonal elements)
 * @param i Subsystem index where the gate \a A is applied
 * @param n Number of qubits
 * @return Gate \a A applied to the qubit \a i of \a state
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] expr_t<Derived1>
apply_psi_1q_diag(const Eigen::MatrixBase<Derived1>& state,
                  const Eigen::MatrixBase<Derived2>& A, idx i, idx n) {
    // A deep copy of the input state is required
    expr_t<Derived1> result = state;

    // Apply the gate in-place on the result copy
    apply_psi_1q_diag_inplace(result, A, i, n);

    return result;
}

/**
 * @brief Applies the 2-qubit diagonal gate \a A to the qubits \a i and \a j of
 * the multi-partite state vector \a state in-place
 *
 * @param state Eigen expression (modified in-place)
 * @param A Eigen expression (4 x 1 vector of diagonal elements)
 * @param i Subsystem index where the gate \a A is applied
 * @param j Subsystem index where the gate \a A is applied
 * @param n Number of qubits
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] void
apply_psi_2q_diag_inplace(Eigen::MatrixBase<Derived1>& state,
                          const Eigen::MatrixBase<Derived2>& A, idx i, idx j,
                          idx n) {
    // Type and Dimension Setup
    using Scalar = typename Derived1::Scalar;

    // Input Validation
    assert(i < n && j < n && i != j &&
           "Target qubit indices i and j must be distinct and less than n");
    // A is a diagonal gate, passed as a 4x1 or 1x4 vector
    assert(A.size() == 4 && "Diagonal gate A must have 4 elements");

#ifndef NDEBUG
    const idx D = static_cast<idx>(std::size_t{1} << n);
    assert(static_cast<idx>(state.size()) == D &&
           "State vector size must be 2^n");
#endif

    // Pre-calculate Bit-Shift Constants (BIG-ENDIAN)
    const idx p_i = n - 1 - i;
    const idx p_j = n - 1 - j;
    const idx s_i = static_cast<idx>(std::size_t{1} << p_i);
    const idx s_j = static_cast<idx>(std::size_t{1} << p_j);

    const idx D_rem = static_cast<idx>(std::size_t{1} << (n - 2));

    // Extract the 4 diagonal elements
    const Scalar a00 = A.coeff(0);
    const Scalar a01 = A.coeff(1);
    const Scalar a10 = A.coeff(2);
    const Scalar a11 = A.coeff(3);

    // Amplitude Transformation
#ifdef QPP_OPENMP
// NOLINTNEXTLINE
#pragma omp parallel for
#endif // QPP_OPENMP
    for (idx block_idx = 0; block_idx < D_rem; ++block_idx) {
        idx k00 = 0;
        idx current_block_bit = 0;

        for (idx p = 0; p < n; ++p) {
            if (p == p_i || p == p_j) {
                continue;
            }

            if ((static_cast<std::size_t>(block_idx) >> current_block_bit) &
                std::size_t{1}) {
                k00 |= static_cast<idx>(std::size_t{1} << p);
            }
            current_block_bit++;
        }

        // Calculate the indices for the 4 basis states
        const idx k01 = k00 + s_j;
        const idx k10 = k00 + s_i;
        const idx k11 = k00 + s_i + s_j;

        // Apply diagonal scaling in-place
        state.coeffRef(k00) *= a00;
        state.coeffRef(k01) *= a01;
        state.coeffRef(k10) *= a10;
        state.coeffRef(k11) *= a11;
    }
}

/**
 * @brief Applies the 2-qubit diagonal gate \a A to the qubits \a i and \a j of
 * the multi-partite state vector \a state
 *
 * @param state Eigen expression
 * @param A Eigen expression (4 x 1 vector of diagonal elements)
 * @param i Subsystem index where the gate \a A is applied
 * @param j Subsystem index where the gate \a A is applied
 * @param n Number of qubits
 * @return Gate \a A applied to the qubit \a i and \a j of \a state
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] expr_t<Derived1>
apply_psi_2q_diag(const Eigen::MatrixBase<Derived1>& state,
                  const Eigen::MatrixBase<Derived2>& A, idx i, idx j, idx n) {
    // A deep copy of the input state is required.
    expr_t<Derived1> result = state;

    // Delegate to in-place implementation
    apply_psi_2q_diag_inplace(result, A, i, j, n);

    return result;
}

/**
 * @brief Applies the 3-qubit diagonal gate \a A to the qubit \a i, \a j, and \a
 * k of the multi-partite state vector \a state in-place
 *
 * @param state Eigen expression (modified in-place)
 * @param A Eigen expression (8 x 1 vector of diagonal elements)
 * @param i Subsystem index where the gate \a A is applied
 * @param j Subsystem index where the gate \a A is applied
 * @param k Subsystem index where the gate \a A is applied
 * @param n Number of qubits
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] void
apply_psi_3q_diag_inplace(Eigen::MatrixBase<Derived1>& state,
                          const Eigen::MatrixBase<Derived2>& A, idx i, idx j,
                          idx k, idx n) {
    // Type and Dimension Setup
    using Scalar = typename Derived1::Scalar;

    // Input Validation
    assert(i < n && j < n && k < n && i != j && i != k && j != k &&
           "Target qubit indices i, j, k must be distinct and less than n");
    // A is a diagonal gate, passed as an 8x1 or 1x8 vector
    assert(A.size() == 8 && "Diagonal gate A must have 8 elements");

#ifndef NDEBUG
    const idx D = static_cast<idx>(std::size_t{1} << n);
    assert(static_cast<idx>(state.size()) == D &&
           "State vector size must be 2^n");
#endif

    // Pre-calculate Bit-Shift Constants (BIG-ENDIAN)
    const idx p_i = n - 1 - i;
    const idx p_j = n - 1 - j;
    const idx p_k = n - 1 - k;

    const idx s_i = static_cast<idx>(std::size_t{1} << p_i);
    const idx s_j = static_cast<idx>(std::size_t{1} << p_j);
    const idx s_k = static_cast<idx>(std::size_t{1} << p_k);

    const idx D_rem = static_cast<idx>(std::size_t{1} << (n - 3));

    // Extract the 8 diagonal elements
    const Scalar a000 = A.coeff(0);
    const Scalar a001 = A.coeff(1);
    const Scalar a010 = A.coeff(2);
    const Scalar a011 = A.coeff(3);
    const Scalar a100 = A.coeff(4);
    const Scalar a101 = A.coeff(5);
    const Scalar a110 = A.coeff(6);
    const Scalar a111 = A.coeff(7);

    // Amplitude Transformation
#ifdef QPP_OPENMP
// NOLINTNEXTLINE
#pragma omp parallel for
#endif // QPP_OPENMP
    for (idx block_idx = 0; block_idx < D_rem; ++block_idx) {
        idx k000 = 0;
        idx current_block_bit = 0;

        for (idx p = 0; p < n; ++p) {
            if (p == p_i || p == p_j || p == p_k) {
                continue;
            }

            if ((static_cast<std::size_t>(block_idx) >> current_block_bit) &
                std::size_t{1}) {
                k000 |= static_cast<idx>(std::size_t{1} << p);
            }
            current_block_bit++;
        }

        // Calculate the indices for the 8 basis states
        const idx k001 = k000 + s_k;
        const idx k010 = k000 + s_j;
        const idx k011 = k010 + s_k;
        const idx k100 = k000 + s_i;
        const idx k101 = k100 + s_k;
        const idx k110 = k100 + s_j;
        const idx k111 = k110 + s_k;

        // Apply diagonal scaling in-place
        state.coeffRef(k000) *= a000;
        state.coeffRef(k001) *= a001;
        state.coeffRef(k010) *= a010;
        state.coeffRef(k011) *= a011;
        state.coeffRef(k100) *= a100;
        state.coeffRef(k101) *= a101;
        state.coeffRef(k110) *= a110;
        state.coeffRef(k111) *= a111;
    }
}

/**
 * @brief Applies the 3-qubit diagonal gate \a A to the qubit \a i, \a j, and \a
 * k of the multi-partite state vector \a state
 *
 * @param state Eigen expression
 * @param A Eigen expression (8 x 1 vector of diagonal elements)
 * @param i Subsystem index where the gate \a A is applied
 * @param j Subsystem index where the gate \a A is applied
 * @param k Subsystem index where the gate \a A is applied
 * @param n Number of qubits
 * @return Gate \a A applied to the qubit \a i, \a j, and \a k of \a state
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] expr_t<Derived1>
apply_psi_3q_diag(const Eigen::MatrixBase<Derived1>& state,
                  const Eigen::MatrixBase<Derived2>& A, idx i, idx j, idx k,
                  idx n) {
    // A deep copy of the input state is required.
    expr_t<Derived1> result = state;

    // Delegate to in-place implementation
    apply_psi_3q_diag_inplace(result, A, i, j, k, n);

    return result;
}

/**
 * @brief Applies the multi-qubit diagonal gate \a A to the part \a target of
 * the multi-partite state vector \a state in-place
 *
 * @param state Eigen expression (modified in-place)
 * @param A Eigen expression (2^k x 1 vector of diagonal elements)
 * @param target Subsystem indexes where the gate \a A is applied
 * @param n Number of qubits
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] void
apply_psi_kq_diag_inplace(Eigen::MatrixBase<Derived1>& state,
                          const Eigen::MatrixBase<Derived2>& A,
                          const std::vector<idx>& target, idx n) {
    // Setup
    const idx k = static_cast<idx>(target.size());
    const idx dim = static_cast<idx>(std::size_t{1} << k);
    const idx outer_dim = static_cast<idx>(std::size_t{1} << (n - k));

    // Input Validation
    assert(static_cast<idx>(A.size()) == dim &&
           "Gate A must be a 2^k element diagonal vector");
#ifndef NDEBUG
    const idx D = static_cast<idx>(std::size_t{1} << n);
    for (idx t : target) {
        assert(t < n && "Target qubit index must be less than n");
    }
    assert(static_cast<idx>(state.size()) == D && state.cols() == 1 &&
           "State must be a 2^n x 1 column vector (ket)");
    if (k > 1) {
        std::vector<idx> sorted_target = target;
        std::sort(sorted_target.begin(), sorted_target.end());
        assert(std::adjacent_find(sorted_target.begin(), sorted_target.end()) ==
                   sorted_target.end() &&
               "Target qubit indices must be distinct");
    }
#endif

    // Pre-calculation of Inner Indices (Gate Index Mapping)
    std::vector<idx> inner_idx(dim);
    for (idx r = 0; r < dim; ++r) {
        idx index_r = 0;
        for (idx j = 0; j < k; ++j) {
            const idx target_qubit_pos = target[k - 1 - j];
            if ((static_cast<std::size_t>(r) >> j) & std::size_t{1}) {
                index_r += static_cast<idx>(std::size_t{1}
                                            << (n - 1 - target_qubit_pos));
            }
        }
        inner_idx[r] = index_r;
    }

    // Pre-calculation of Outer Indices (Spectator Block Base Index)
    std::vector<bool> is_qubit_acted_on(n, false);
    for (idx q : target) {
        is_qubit_acted_on[q] = true;
    }
    std::vector<idx> spectator_qubits;
    for (idx i = 0; i < n; ++i) {
        if (!is_qubit_acted_on[i]) {
            spectator_qubits.push_back(i);
        }
    }

    std::vector<idx> outer_idx(outer_dim);
    for (idx m = 0; m < outer_dim; ++m) {
        idx i_base = 0;
        for (idx j = 0; j < static_cast<idx>(spectator_qubits.size()); ++j) {
            if ((static_cast<std::size_t>(m) >> j) & std::size_t{1}) {
                i_base += static_cast<idx>(std::size_t{1}
                                           << (n - 1 - spectator_qubits[j]));
            }
        }
        outer_idx[m] = i_base;
    }

    // Main Parallel Loop
#ifdef QPP_OPENMP
// NOLINTNEXTLINE
#pragma omp parallel for
#endif // QPP_OPENMP
    for (idx m = 0; m < outer_dim; ++m) {
        const idx i_base = outer_idx[m];

        // Diagonal Transformation
        // Since A is diagonal, we scale each amplitude directly.
        // This avoids the O(dim^2) matrix multiplication and memory
        // allocations.
        for (idx r = 0; r < dim; ++r) {
            state.coeffRef(i_base + inner_idx[r]) *= A.coeff(r);
        }
    }
}

/**
 * @brief Applies the multi-qubit diagonal gate \a A to the part \a target of
 * the multi-partite state vector \a state
 *
 * @param state Eigen expression
 * @param A Eigen expression (2^k x 1 vector of diagonal elements)
 * @param target Subsystem indexes where the gate \a A is applied
 * @param n Number of qubits
 * @return Gate \a A applied to the part \a target of \a state
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] qpp::expr_t<Derived1>
apply_psi_kq_diag(const Eigen::MatrixBase<Derived1>& state,
                  const Eigen::MatrixBase<Derived2>& A,
                  const std::vector<idx>& target, idx n) {
    // Create the output state vector as a copy of the input.
    qpp::expr_t<Derived1> result = state;

    // Delegate to in-place implementation
    apply_psi_kq_diag_inplace(result, A, target, n);

    return result;
}

/**
 * @brief Applies the 1-qubit diagonal gate \a A to the qubit \a i of the
 * multi-partite density matrix \a state in-place
 *
 * @param state Eigen expression (modified in-place)
 * @param A Eigen expression (2 x 1 vector of diagonal elements)
 * @param i Subsystem index where the gate \a A is applied
 * @param n Number of qubits
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] void
apply_rho_1q_diag_inplace(Eigen::MatrixBase<Derived1>& state,
                          const Eigen::MatrixBase<Derived2>& A, idx i, idx n) {
    using Scalar = typename Derived1::Scalar;
    const idx D = static_cast<idx>(std::size_t{1} << n);

    // Input Validation
    assert(static_cast<idx>(state.rows()) == D &&
           static_cast<idx>(state.cols()) == D &&
           "State must be a square matrix sized 2^n x 2^n");
    assert(i < n && "Target qubit index i must be less than n");
    assert(A.size() == 2 && "Diagonal gate A must have 2 elements");

    // Bit position for qubit i (p_i) and the stride (s_i = 2^p_i)
    const idx p_i = n - 1 - i;
    const idx s_i = static_cast<idx>(std::size_t{1} << p_i);

    // Constants for index reconstruction
    const idx D_spec = D / 2; // Total number of spectator states (2^(n-1))
    const idx p_i_plus_1 = p_i + 1;
    const idx low_mask = s_i - 1; // Mask for bits below p_i

    // Extract diagonal elements and their conjugates
    const Scalar a0 = A.coeff(0);
    const Scalar a1 = A.coeff(1);
    const Scalar a0_conj = std::conj(a0);
    const Scalar a1_conj = std::conj(a1);

    // Pre-calculate products for the 2x2 block transformation
    const Scalar g00 = a0 * a0_conj; // |0><0| scaling
    const Scalar g01 = a0 * a1_conj; // |0><1| scaling
    const Scalar g10 = a1 * a0_conj; // |1><0| scaling
    const Scalar g11 = a1 * a1_conj; // |1><1| scaling

    const idx total_iterations = D_spec * D_spec;
#ifdef QPP_OPENMP
// NOLINTNEXTLINE
#pragma omp parallel for
#endif // QPP_OPENMP
    for (idx iter = 0; iter < total_iterations; ++iter) {
        const idx s = iter / D_spec;
        const idx s_prime = iter % D_spec;

        // Calculate Row Indices (r0, r1)
        const idx s_high_r =
            static_cast<idx>(static_cast<std::size_t>(s) >> p_i);
        const idx s_low_r = static_cast<idx>(
            static_cast<std::size_t>(s) & static_cast<std::size_t>(low_mask));
        const idx r0 = static_cast<idx>(
            (static_cast<std::size_t>(s_high_r) << p_i_plus_1) |
            static_cast<std::size_t>(s_low_r));
        const idx r1 = r0 + s_i;

        // Calculate Column Indices (c0, c1)
        const idx s_prime_high_c =
            static_cast<idx>(static_cast<std::size_t>(s_prime) >> p_i);
        const idx s_prime_low_c =
            static_cast<idx>(static_cast<std::size_t>(s_prime) &
                             static_cast<std::size_t>(low_mask));
        const idx c0 = static_cast<idx>(
            (static_cast<std::size_t>(s_prime_high_c) << p_i_plus_1) |
            static_cast<std::size_t>(s_prime_low_c));
        const idx c1 = c0 + s_i;

        // Apply transformation: rho'(r,c) = A(r) * rho(r,c) * conj(A(c))
        state.coeffRef(r0, c0) *= g00;
        state.coeffRef(r0, c1) *= g01;
        state.coeffRef(r1, c0) *= g10;
        state.coeffRef(r1, c1) *= g11;
    }
}

/**
 * @brief Applies the 1-qubit diagonal gate \a A to the qubit \a i of the
 * multi-partite density matrix \a state
 *
 * @param state Eigen expression
 * @param A Eigen expression (2 x 1 vector of diagonal elements)
 * @param i Subsystem index where the gate \a A is applied
 * @param n Number of qubits
 * @return Gate \a A applied to the qubit \a i of \a state
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] expr_t<Derived1>
apply_rho_1q_diag(const Eigen::MatrixBase<Derived1>& state,
                  const Eigen::MatrixBase<Derived2>& A, idx i, idx n) {
    // A deep copy of the input state is required.
    expr_t<Derived1> result = state;

    // Delegate to in-place implementation
    apply_rho_1q_diag_inplace(result, A, i, n);

    return result;
}

/**
 * @brief Applies the 2-qubit diagonal gate \a A to the qubits \a i and \a j of
 * the multi-partite density matrix \a state in-place
 *
 * @param state Eigen expression (modified in-place)
 * @param A Eigen expression (4 x 1 vector of diagonal elements)
 * @param i Subsystem index where the gate \a A is applied
 * @param j Subsystem index where the gate \a A is applied
 * @param n Number of qubits
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] void
apply_rho_2q_diag_inplace(Eigen::MatrixBase<Derived1>& state,
                          const Eigen::MatrixBase<Derived2>& A, idx i, idx j,
                          idx n) {
    using Scalar = typename Derived1::Scalar;
    const idx D [[maybe_unused]] = static_cast<idx>(state.rows());

    // Input Validation
    assert(i < n && j < n && i != j &&
           "Target qubit indices i and j must be distinct and less than n");
    assert(A.size() == 4 && "Diagonal gate A must have 4 elements");
#ifndef NDEBUG
    const idx D_expected = static_cast<idx>(std::size_t{1} << n);
    assert(D == D_expected && D == static_cast<idx>(state.cols()) &&
           "State must be a square matrix sized 2^n x 2^n");
#endif

    // Pre-calculate scaling factors G(r, c) = a_r * conj(a_c)
    Eigen::Matrix<Scalar, 4, 4> G;
    for (idx r = 0; r < 4; ++r) {
        Scalar a_r = A.coeff(r);
        for (idx c = 0; c < 4; ++c) {
            G(r, c) = a_r * std::conj(A.coeff(c));
        }
    }

    // Index Conversion (Big Endian)
    const idx p_i = n - 1 - i;
    const idx p_j = n - 1 - j;
    const idx s_i = static_cast<idx>(std::size_t{1} << p_i);
    const idx s_j = static_cast<idx>(std::size_t{1} << p_j);

    const idx D_rest = static_cast<idx>(std::size_t{1} << (n - 2));

    // Block Iteration (Parallelized over spectator rows and columns)
#ifdef QPP_OPENMP
// NOLINTNEXTLINE
#pragma omp parallel for collapse(2)
#endif // QPP_OPENMP
    for (idx r = 0; r < D_rest; ++r) {
        for (idx c = 0; c < D_rest; ++c) {
            idx r_base_row = 0;
            idx r_base_col = 0;
            idx current_r = r;
            idx current_c = c;

            // Reconstruct base indices where target bits are 0
            for (idx q = 0; q < n; ++q) {
                if (q != p_i && q != p_j) {
                    if (current_r & 1) {
                        r_base_row |= (idx{1} << q);
                    }
                    if (current_c & 1) {
                        r_base_col |= (idx{1} << q);
                    }
                    current_r >>= 1;
                    current_c >>= 1;
                }
            }

            const idx v_row[4] = {r_base_row, r_base_row + s_j,
                                  r_base_row + s_i, r_base_row + s_i + s_j};
            const idx v_col[4] = {r_base_col, r_base_col + s_j,
                                  r_base_col + s_i, r_base_col + s_i + s_j};

            // Apply transformation directly: rho_rc' = rho_rc * G_rc
            for (int row = 0; row < 4; ++row) {
                for (int col = 0; col < 4; ++col) {
                    state.coeffRef(v_row[row], v_col[col]) *= G(row, col);
                }
            }
        }
    }
}

/**
 * @brief Applies the 2-qubit diagonal gate \a A to the qubits \a i and \a j of
 * the multi-partite density matrix \a state
 *
 * @param state Eigen expression
 * @param A Eigen expression (4 x 1 vector of diagonal elements)
 * @param i Subsystem index where the gate \a A is applied
 * @param j Subsystem index where the gate \a A is applied
 * @param n Number of qubits
 * @return Gate \a A applied to the qubits \a i and \a j of \a state
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] expr_t<Derived1>
apply_rho_2q_diag(const Eigen::MatrixBase<Derived1>& state,
                  const Eigen::MatrixBase<Derived2>& A, idx i, idx j, idx n) {
    // Deep copy required for functional interface
    expr_t<Derived1> result = state;

    // Delegate to the in-place version
    apply_rho_2q_diag_inplace(result, A, i, j, n);

    return result;
}

/**
 * @brief Applies the 3-qubit diagonal gate \a A to the qubit \a i, \a j, and \a
 * k of the multi-partite density matrix \a state in-place
 *
 * @param state Eigen expression (modified in-place)
 * @param A Eigen expression (8 x 1 vector of diagonal elements)
 * @param i Subsystem index where the gate \a A is applied
 * @param j Subsystem index where the gate \a A is applied
 * @param k Subsystem index where the gate \a A is applied
 * @param n Number of qubits
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] void
apply_rho_3q_diag_inplace(Eigen::MatrixBase<Derived1>& state,
                          const Eigen::MatrixBase<Derived2>& A, idx i, idx j,
                          idx k, idx n) {
    using Scalar = typename Derived1::Scalar;
    const idx D [[maybe_unused]] = static_cast<idx>(state.rows());

    // Input Validation
    assert(i < n && j < n && k < n && i != j && i != k && j != k &&
           "Target qubit indices i, j, and k must be distinct and less than n");
    assert(A.size() == 8 && "Diagonal gate A must have 8 elements");
    assert(n >= 3 && "Need at least 3 qubits for a 3-qubit diagonal gate");
#ifndef NDEBUG
    const idx D_expected = static_cast<idx>(std::size_t{1} << n);
    assert(D == D_expected && D == static_cast<idx>(state.cols()) &&
           "State must be a square matrix of size 2^n x 2^n");
#endif

    // Pre-calculate scaling factors G(r, c) = a_r * conj(a_c)
    Eigen::Matrix<Scalar, 8, 8> G;
    for (idx r = 0; r < 8; ++r) {
        Scalar a_r = A.coeff(r);
        for (idx c = 0; c < 8; ++c) {
            G(r, c) = a_r * std::conj(A.coeff(c));
        }
    }

    // Index Conversion (Big Endian)
    const idx p_i = n - 1 - i;
    const idx p_j = n - 1 - j;
    const idx p_k = n - 1 - k;
    const idx s_i = static_cast<idx>(std::size_t{1} << p_i);
    const idx s_j = static_cast<idx>(std::size_t{1} << p_j);
    const idx s_k = static_cast<idx>(std::size_t{1} << p_k);

    const idx D_rest = static_cast<idx>(std::size_t{1} << (n - 3));

#ifdef QPP_OPENMP
// NOLINTNEXTLINE
#pragma omp parallel for collapse(2)
#endif // QPP_OPENMP
    for (idx r = 0; r < D_rest; ++r) {
        for (idx c = 0; c < D_rest; ++c) {
            idx r_base_row = 0;
            idx r_base_col = 0;
            idx current_r = r;
            idx current_c = c;

            // Reconstruct base indices where target bits are 0
            for (idx q = 0; q < n; ++q) {
                if (q != p_i && q != p_j && q != p_k) {
                    if (current_r & 1) {
                        r_base_row |= (idx{1} << q);
                    }
                    if (current_c & 1) {
                        r_base_col |= (idx{1} << q);
                    }
                    current_r >>= 1;
                    current_c >>= 1;
                }
            }

            const idx v_row[8] = {r_base_row,
                                  r_base_row + s_k,
                                  r_base_row + s_j,
                                  r_base_row + s_j + s_k,
                                  r_base_row + s_i,
                                  r_base_row + s_i + s_k,
                                  r_base_row + s_i + s_j,
                                  r_base_row + s_i + s_j + s_k};
            const idx v_col[8] = {r_base_col,
                                  r_base_col + s_k,
                                  r_base_col + s_j,
                                  r_base_col + s_j + s_k,
                                  r_base_col + s_i,
                                  r_base_col + s_i + s_k,
                                  r_base_col + s_i + s_j,
                                  r_base_col + s_i + s_j + s_k};

            // Apply transformation directly: rho_rc' = rho_rc * G_rc
            for (int row = 0; row < 8; ++row) {
                for (int col = 0; col < 8; ++col) {
                    state.coeffRef(v_row[row], v_col[col]) *= G(row, col);
                }
            }
        }
    }
}

/**
 * @brief Applies the 3-qubit diagonal gate \a A to the qubit \a i, \a j, and \a
 * k of the multi-partite density matrix \a state
 *
 * @param state Eigen expression
 * @param A Eigen expression (8 x 1 vector of diagonal elements)
 * @param i Subsystem index where the gate \a A is applied
 * @param j Subsystem index where the gate \a A is applied
 * @param k Subsystem index where the gate \a A is applied
 * @param n Number of qubits
 * @return Gate \a A applied to the qubit \a i, \a j, and \a k of \a state
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] expr_t<Derived1>
apply_rho_3q_diag(const Eigen::MatrixBase<Derived1>& state,
                  const Eigen::MatrixBase<Derived2>& A, idx i, idx j, idx k,
                  idx n) {
    // Create copy for functional behavior
    expr_t<Derived1> result = state;

    // Perform the operation in-place on the copy
    apply_rho_3q_diag_inplace(result, A, i, j, k, n);

    return result;
}

/**
 * @brief Applies the multi-qubit diagonal gate \a A to the part \a target of
 * the multi-partite density matrix \a state in-place
 *
 * @param state Eigen expression (modified in-place)
 * @param A Eigen expression (2^k x 1 vector of diagonal elements)
 * @param target Subsystem indexes where the gate \a A is applied
 * @param n Number of qubits
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] void
apply_rho_kq_diag_inplace(Eigen::MatrixBase<Derived1>& state,
                          const Eigen::MatrixBase<Derived2>& A,
                          const std::vector<idx>& target, idx n) {
    using Scalar = typename Derived1::Scalar;
    const idx D [[maybe_unused]] = static_cast<idx>(state.rows());
    const idx k = target.size();
    const idx D_k = (k == 0) ? 1 : static_cast<idx>(std::size_t{1} << k);

    // Input Validation
    assert(static_cast<idx>(A.size()) == D_k &&
           "Gate A must have 2^k elements");
    assert(static_cast<idx>(state.rows()) == D &&
           static_cast<idx>(state.cols()) == D &&
           "State must be a 2^n x 2^n matrix");
#ifndef NDEBUG
    for (idx t : target) {
        assert(t < n && "Target qubit index must be less than n");
    }
    if (k > 1) {
        std::vector<idx> sorted_target = target;
        std::sort(sorted_target.begin(), sorted_target.end());
        assert(std::unique(sorted_target.begin(), sorted_target.end()) ==
                   sorted_target.end() &&
               "Target qubit indices must be distinct");
    }
#endif

    // Pre-calculate scaling factors G(r, c) = a_r * conj(a_c)
    // Using Dynamic Eigen matrix to handle arbitrary k
    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> G(D_k, D_k);
    for (idx r = 0; r < D_k; ++r) {
        Scalar a_r = A.coeff(r);
        for (idx c = 0; c < D_k; ++c) {
            G(r, c) = a_r * std::conj(A.coeff(c));
        }
    }

    // Pre-Calculations for indices
    std::vector<idx> target_phys(k);
    std::vector<idx> P_gate_basis(k);
    for (idx l = 0; l < k; ++l) {
        target_phys[l] = n - target[l] - 1;
        P_gate_basis[l] = static_cast<idx>(std::size_t{1} << target_phys[l]);
    }

    std::vector<idx> rest_phys = fast_complement(target_phys, n);
    const idx D_rest = static_cast<idx>(std::size_t{1} << (n - k));

    // Pre-compute relative target components to avoid repeated bit-shifting
    std::vector<idx> target_components(D_k, 0);
    for (idx m = 0; m < D_k; ++m) {
        for (idx l = 0; l < k; ++l) {
            idx bit_pos = k - 1 - l;
            if ((m >> bit_pos) & 1) {
                target_components[m] += P_gate_basis[l];
            }
        }
    }

    // Block Iteration
#ifdef QPP_OPENMP
// NOLINTNEXTLINE
#pragma omp parallel for collapse(2)
#endif // QPP_OPENMP
    for (idx r = 0; r < D_rest; ++r) {
        for (idx c = 0; c < D_rest; ++c) {
            idx r_base_row = 0;
            idx r_base_col = 0;
            idx current_r = r;
            idx current_c = c;

            for (const auto& q_phys : rest_phys) {
                if (current_r & 1) {
                    r_base_row |= (idx{1} << q_phys);
                }
                if (current_c & 1) {
                    r_base_col |= (idx{1} << q_phys);
                }
                current_r >>= 1;
                current_c >>= 1;
            }

            // Apply transformation directly using pre-computed components and G
            for (idx row = 0; row < D_k; ++row) {
                const idx actual_row = r_base_row + target_components[row];
                for (idx col = 0; col < D_k; ++col) {
                    state.coeffRef(actual_row,
                                   r_base_col + target_components[col]) *=
                        G(row, col);
                }
            }
        }
    }
}

/**
 * @brief Applies the multi-qubit diagonal gate \a A to the part \a target of
 * the multi-partite density matrix \a state
 *
 * @param state Eigen expression
 * @param A Eigen expression (2^k x 1 vector of diagonal elements)
 * @param target Subsystem indexes where the gate \a A is applied
 * @param n Number of qubits
 * @return Gate \a A applied to the part \a target of \a state
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] expr_t<Derived1>
apply_rho_kq_diag(const Eigen::MatrixBase<Derived1>& state,
                  const Eigen::MatrixBase<Derived2>& A,
                  const std::vector<idx>& target, idx n) {
    // Functional deep copy
    expr_t<Derived1> result = state;

    // Apply transformation in-place on the copy
    apply_rho_kq_diag_inplace(result, A, target, n);

    return result;
}
} // namespace qpp::internal::kernels::qubit

#endif /* QPP_INTERNAL_KERNELS_QUBIT_APPLY_DIAG_HPP_ */
