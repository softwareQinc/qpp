/*
 * This file is part of Quantum++.
 *
 * Copyright (c) 2017 - 2025 softwareQ Inc. All rights reserved.
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
 * \file qpp/internal/critical.hpp
 * \brief Internal highly optimized critical functions
 */

#ifndef QPP_INTERNAL_CRITICAL_HPP_
#define QPP_INTERNAL_CRITICAL_HPP_

#include <cassert>
#include <vector>

#ifndef NDEBUG
#include <set>
#endif

#include <Eigen/Dense>

#include "qpp/types.hpp"

namespace qpp {
namespace internal {
/**
 * \brief Applies the 1-qubit gate \a A to the qubit \a i of the
 * multi-partite state vector \a state
 *
 * \param state Eigen expression
 * \param A Eigen expression (2x2 matrix)
 * \param i Subsystem index where the gate \a A is applied
 * \param n Number of qubits
 * \return Gate \a A applied to the qubit \a i of \a state
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] expr_t<Derived1>
apply_psi_1q(const Eigen::MatrixBase<Derived1>& state,
             const Eigen::MatrixBase<Derived2>& A, idx i, idx n) {
    // Type and Dimension Setup
    using Scalar = typename Derived1::Scalar;
    const idx D = 1ULL << n; // Total size of the state vector (2^n)

    // Input Validation
#ifndef NDEBUG
    assert(i < n && "Target qubit index i must be less than n");
    assert(static_cast<idx>(state.size()) == D &&
           "State vector size must be 2^n");
    assert(A.rows() == 2 && A.cols() == 2 && "Gate A must be a 2x2 matrix");
#endif

    // A deep copy of the input state is required. The new state must be
    // calculated using the *old* values of psi_k0 and psi_k1 simultaneously.
    expr_t<Derived1> result = state;

    // Pre-calculate Bit-Shift Constants (BIG-ENDIAN)
    // The effective bit position 'j' corresponding to qubit 'i' in big-endian
    // is (n - 1 - i).
    const idx j = n - 1 - i;

    // 'step' is 2^j. This is the difference in index between |...0...> and
    // |...1...> at bit j.
    const idx step = 1ULL << j;

    // 'jump' is 2^(j+1). This is the size of the block that repeats.
    const idx jump = 1ULL << (j + 1);

    // Extract the 2x2 gate elements
    const Scalar a00 = A.coeff(0, 0);
    const Scalar a01 = A.coeff(0, 1);
    const Scalar a10 = A.coeff(1, 0);
    const Scalar a11 = A.coeff(1, 1);

    // Pair-wise Amplitude Transformation
    // The outer loop (L) iterates over all blocks of size 'jump'.
    // This loop is perfectly independent and is the primary target for
    // parallelization.
#ifdef QPP_OPENMP
// NOLINTNEXTLINE
#pragma omp parallel for
#endif // QPP_OPENMP
    for (idx L = 0; L < D; L += jump) {
        // The inner loop (R) iterates over the lower part of the block, from 0
        // to (step - 1).
        for (idx R = 0; R < step; ++R) {
            // k0 is the index where the j-th bit is 0: L + R
            const idx k0 = L + R;
            // k1 is the index where the j-th bit is 1: k0 + step
            const idx k1 = k0 + step;

            // Fetch the ORIGINAL amplitudes from the INPUT 'state'
            // Using const reference for safety and potential optimization
            const Scalar& psi_k0 = state.coeff(k0);
            const Scalar& psi_k1 = state.coeff(k1);

            // Apply the 2x2 matrix-vector multiplication and write to 'result'
            // psi'_k0 = A[0,0] * psi_k0 + A[0,1] * psi_k1
            result.coeffRef(k0) = (a00 * psi_k0) + (a01 * psi_k1);

            // psi'_k1 = A[1,0] * psi_k0 + A[1,1] * psi_k1
            result.coeffRef(k1) = (a10 * psi_k0) + (a11 * psi_k1);
        }
    }

    return result;
}

/**
 * \brief Applies the 2-qubit gate \a A to the qubits \a i and \a j of the
 * multi-partite state vector \a state
 *
 * \param state Eigen expression
 * \param A Eigen expression (4x4 matrix)
 * \param i Subsystem index where the gate \a A is applied
 * \param j Subsystem index where the gate \a A is applied
 * \param n Number of qubits
 * \return Gate \a A applied to the qubit \a i and \a j of \a state
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] expr_t<Derived1>
apply_psi_2q(const Eigen::MatrixBase<Derived1>& state,
             const Eigen::MatrixBase<Derived2>& A, idx i, idx j, idx n) {
    // Type and Dimension Setup
    using Scalar = typename Derived1::Scalar;

    // Input Validation
#ifndef NDEBUG
    const idx D = 1ULL << n; // Total size of the state vector (2^n)
    assert(i < n && j < n && i != j &&
           "Target qubit indices i and j must be distinct and less than n");
    assert(static_cast<idx>(state.size()) == D &&
           "State vector size must be 2^n");
    assert(A.rows() == 4 && A.cols() == 4 && "Gate A must be a 4x4 matrix");
#endif

    // A deep copy of the input state is required.
    expr_t<Derived1> result = state;

    // Pre-calculate Bit-Shift Constants (BIG-ENDIAN)
    // The effective bit position 'p' for a physical qubit 'q' in big-endian is
    // (n - 1 - q).
    const idx p_i = n - 1 - i;
    const idx p_j = n - 1 - j;

    // 's_i' and 's_j' are the power-of-2 values that flip the bit at p_i and
    // p_j.
    const idx s_i = 1ULL << p_i;
    const idx s_j = 1ULL << p_j;

    // The number of remaining bits that define the independent blocks (2^(n-2)
    // total blocks)
    const idx D_rem = 1ULL << (n - 2);

    // Extract the 4x4 gate elements
    // Using .coeff() for fast, direct access. A is assumed to be 4x4.
    const Scalar a00 = A.coeff(0, 0);
    const Scalar a01 = A.coeff(0, 1);
    const Scalar a02 = A.coeff(0, 2);
    const Scalar a03 = A.coeff(0, 3);

    const Scalar a10 = A.coeff(1, 0);
    const Scalar a11 = A.coeff(1, 1);
    const Scalar a12 = A.coeff(1, 2);
    const Scalar a13 = A.coeff(1, 3);

    const Scalar a20 = A.coeff(2, 0);
    const Scalar a21 = A.coeff(2, 1);
    const Scalar a22 = A.coeff(2, 2);
    const Scalar a23 = A.coeff(2, 3);

    const Scalar a30 = A.coeff(3, 0);
    const Scalar a31 = A.coeff(3, 1);
    const Scalar a32 = A.coeff(3, 2);
    const Scalar a33 = A.coeff(3, 3);

    // Pair-wise Amplitude Transformation (Robust Loop for all i, j)
#ifdef QPP_OPENMP
// NOLINTNEXTLINE
#pragma omp parallel for
#endif // QPP_OPENMP
    // block_idx iterates over all 2^(n-2) base indices where bits p_i and p_j
    // are zero
    for (idx block_idx = 0; block_idx < D_rem; ++block_idx) {

        // Reconstruct the base index k00
        idx k00 = 0;
        idx current_block_bit =
            0; // The current bit position in 'block_idx' (0 to n-3)

        // Iterate over all n physical bit positions (p=0 to n-1, p=0 is
        // MSB/q_n-1)
        for (idx p = 0; p < n; ++p) {
            // Skip the two target bit positions p_i and p_j. These bits are
            // always 0 in k00.
            if (p == p_i || p == p_j) {
                continue;
            }

            // If the 'current_block_bit' of 'block_idx' is 1, set the p-th bit
            // in k00
            if ((block_idx >> current_block_bit) & 1) {
                k00 |= (1ULL << p);
            }
            current_block_bit++;
        }

        // Calculate the other three indices from k00, s_i, and s_j

        // k01: |...0_i 1_j...> (A index 1)
        const idx k01 = k00 + s_j;
        // k10: |...1_i 0_j...> (A index 2)
        const idx k10 = k00 + s_i;
        // k11: |...1_i 1_j...> (A index 3)
        const idx k11 = k00 + s_i + s_j;

        // Fetch the ORIGINAL amplitudes
        const Scalar psi00 = state.coeff(k00);
        const Scalar psi01 = state.coeff(k01);
        const Scalar psi10 = state.coeff(k10);
        const Scalar psi11 = state.coeff(k11);

        // Apply the 4x4 matrix-vector multiplication A * [psi00, psi01,
        // psi10, psi11]^T

        // Row 0 -> result_k00 (|0_i 0_j>)
        result.coeffRef(k00) =
            (a00 * psi00) + (a01 * psi01) + (a02 * psi10) + (a03 * psi11);

        // Row 1 -> result_k01 (|0_i 1_j>)
        result.coeffRef(k01) =
            (a10 * psi00) + (a11 * psi01) + (a12 * psi10) + (a13 * psi11);

        // Row 2 -> result_k10 (|1_i 0_j>)
        result.coeffRef(k10) =
            (a20 * psi00) + (a21 * psi01) + (a22 * psi10) + (a23 * psi11);

        // Row 3 -> result_k11 (|1_i 1_j>)
        result.coeffRef(k11) =
            (a30 * psi00) + (a31 * psi01) + (a32 * psi10) + (a33 * psi11);
    }

    return result;
}

/**
 * \brief Applies the 3-qubit gate \a A to the qubit \a i, \a j, and \a k of the
 * multi-partite state vector \a state
 *
 * \param state Eigen expression
 * \param A Eigen expression (8x8 matrix)
 * \param i Subsystem index where the gate \a A is applied
 * \param j Subsystem index where the gate \a A is applied
 * \param k Subsystem index where the gate \a A is applied
 * \param n Number of qubits
 * \return Gate \a A applied to the qubit \a i, \a j, and \a k of \a state
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] expr_t<Derived1>
apply_psi_3q(const Eigen::MatrixBase<Derived1>& state,
             const Eigen::MatrixBase<Derived2>& A, idx i, idx j, idx k, idx n) {
    // Type and Dimension Setup
    using Scalar = typename Derived1::Scalar;

    // Input Validation
#ifndef NDEBUG
    const idx D = 1ULL << n; // Total size of the state vector (2^n)
    assert(i < n && j < n && k < n && i != j && i != k && j != k &&
           "Target qubit indices i, j, k must be distinct and less than n");
    assert(static_cast<idx>(state.size()) == D &&
           "State vector size must be 2^n");
    assert(A.rows() == 8 && A.cols() == 8 && "Gate A must be an 8x8 matrix");
#endif

    // A deep copy of the input state is required.
    expr_t<Derived1> result = state;

    // Pre-calculate Bit-Shift Constants (BIG-ENDIAN)
    // The effective bit position 'p' for a physical qubit 'q' in big-endian is
    // (n - 1 - q).
    const idx p_i = n - 1 - i;
    const idx p_j = n - 1 - j;
    const idx p_k = n - 1 - k;

    // 's_i', 's_j', 's_k' are the power-of-2 values that flip the bit at p_i,
    // p_j, p_k.
    const idx s_i = 1ULL << p_i;
    const idx s_j = 1ULL << p_j;
    const idx s_k = 1ULL << p_k;

    // The number of remaining bits that define the independent blocks (2^(n-3)
    // total blocks)
    const idx D_rem = 1ULL << (n - 3);

    // Extract the 8x8 gate elements (A is assumed to be 8x8)
    // A[r, c] maps to basis states |c_i c_j c_k> -> |r_i r_j r_k>
    const Scalar A00 = A.coeff(0, 0), A01 = A.coeff(0, 1), A02 = A.coeff(0, 2),
                 A03 = A.coeff(0, 3);
    const Scalar A04 = A.coeff(0, 4), A05 = A.coeff(0, 5), A06 = A.coeff(0, 6),
                 A07 = A.coeff(0, 7);

    const Scalar A10 = A.coeff(1, 0), A11 = A.coeff(1, 1), A12 = A.coeff(1, 2),
                 A13 = A.coeff(1, 3);
    const Scalar A14 = A.coeff(1, 4), A15 = A.coeff(1, 5), A16 = A.coeff(1, 6),
                 A17 = A.coeff(1, 7);

    const Scalar A20 = A.coeff(2, 0), A21 = A.coeff(2, 1), A22 = A.coeff(2, 2),
                 A23 = A.coeff(2, 3);
    const Scalar A24 = A.coeff(2, 4), A25 = A.coeff(2, 5), A26 = A.coeff(2, 6),
                 A27 = A.coeff(2, 7);

    const Scalar A30 = A.coeff(3, 0), A31 = A.coeff(3, 1), A32 = A.coeff(3, 2),
                 A33 = A.coeff(3, 3);
    const Scalar A34 = A.coeff(3, 4), A35 = A.coeff(3, 5), A36 = A.coeff(3, 6),
                 A37 = A.coeff(3, 7);

    const Scalar A40 = A.coeff(4, 0), A41 = A.coeff(4, 1), A42 = A.coeff(4, 2),
                 A43 = A.coeff(4, 3);
    const Scalar A44 = A.coeff(4, 4), A45 = A.coeff(4, 5), A46 = A.coeff(4, 6),
                 A47 = A.coeff(4, 7);

    const Scalar A50 = A.coeff(5, 0), A51 = A.coeff(5, 1), A52 = A.coeff(5, 2),
                 A53 = A.coeff(5, 3);
    const Scalar A54 = A.coeff(5, 4), A55 = A.coeff(5, 5), A56 = A.coeff(5, 6),
                 A57 = A.coeff(5, 7);

    const Scalar A60 = A.coeff(6, 0), A61 = A.coeff(6, 1), A62 = A.coeff(6, 2),
                 A63 = A.coeff(6, 3);
    const Scalar A64 = A.coeff(6, 4), A65 = A.coeff(6, 5), A66 = A.coeff(6, 6),
                 A67 = A.coeff(6, 7);

    const Scalar A70 = A.coeff(7, 0), A71 = A.coeff(7, 1), A72 = A.coeff(7, 2),
                 A73 = A.coeff(7, 3);
    const Scalar A74 = A.coeff(7, 4), A75 = A.coeff(7, 5), A76 = A.coeff(7, 6),
                 A77 = A.coeff(7, 7);

    // Robust Amplitude Transformation
#ifdef QPP_OPENMP
// NOLINTNEXTLINE
#pragma omp parallel for
#endif // QPP_OPENMP
    // block_idx iterates over all 2^(n-3) base indices where bits p_i, p_j, p_k
    // are zero
    for (idx block_idx = 0; block_idx < D_rem; ++block_idx) {

        // Reconstruct the base index k000
        idx k000 = 0;
        idx current_block_bit =
            0; // The current bit position in 'block_idx' (0 to n-4)

        // Iterate over all n physical bit positions (p=0 to n-1, p=0 is
        // MSB/q_n-1)
        for (idx p = 0; p < n; ++p) {
            // Skip the three target bit positions p_i, p_j, p_k. These bits are
            // always 0 in k000.
            if (p == p_i || p == p_j || p == p_k) {
                continue;
            }

            // If the 'current_block_bit' of 'block_idx' is 1, set the p-th bit
            // in k000
            if ((block_idx >> current_block_bit) & 1) {
                k000 |= (1ULL << p);
            }
            current_block_bit++;
        }

        // Calculate the other seven indices from k000
        // k0xy is k000 + (x * s_j) + (y * s_k)
        // k1xy is k000 + s_i + (x * s_j) + (y * s_k)

        // k000 is defined
        const idx k001 = k000 + s_k; // |001>
        const idx k010 = k000 + s_j; // |010>
        const idx k011 = k010 + s_k; // |011>

        const idx k100 = k000 + s_i; // |100>
        const idx k101 = k100 + s_k; // |101>
        const idx k110 = k100 + s_j; // |110>
        const idx k111 = k110 + s_k; // |111>

        // Fetch the ORIGINAL amplitudes from the INPUT 'state' (8
        // inputs)
        const Scalar psi000 = state.coeff(k000);
        const Scalar psi001 = state.coeff(k001);
        const Scalar psi010 = state.coeff(k010);
        const Scalar psi011 = state.coeff(k011);
        const Scalar psi100 = state.coeff(k100);
        const Scalar psi101 = state.coeff(k101);
        const Scalar psi110 = state.coeff(k110);
        const Scalar psi111 = state.coeff(k111);

        // Apply the 8x8 matrix-vector multiplication A * [psi000, ...,
        // psi111]^T; grouping multiplications in parentheses for
        // consistency.

        // Row 0 -> result_k000 (|000>)
        result.coeffRef(k000) =
            (A00 * psi000) + (A01 * psi001) + (A02 * psi010) + (A03 * psi011) +
            (A04 * psi100) + (A05 * psi101) + (A06 * psi110) + (A07 * psi111);

        // Row 1 -> result_k001 (|001>)
        result.coeffRef(k001) =
            (A10 * psi000) + (A11 * psi001) + (A12 * psi010) + (A13 * psi011) +
            (A14 * psi100) + (A15 * psi101) + (A16 * psi110) + (A17 * psi111);

        // Row 2 -> result_k010 (|010>)
        result.coeffRef(k010) =
            (A20 * psi000) + (A21 * psi001) + (A22 * psi010) + (A23 * psi011) +
            (A24 * psi100) + (A25 * psi101) + (A26 * psi110) + (A27 * psi111);

        // Row 3 -> result_k011 (|011>)
        result.coeffRef(k011) =
            (A30 * psi000) + (A31 * psi001) + (A32 * psi010) + (A33 * psi011) +
            (A34 * psi100) + (A35 * psi101) + (A36 * psi110) + (A37 * psi111);

        // Row 4 -> result_k100 (|100>)
        result.coeffRef(k100) =
            (A40 * psi000) + (A41 * psi001) + (A42 * psi010) + (A43 * psi011) +
            (A44 * psi100) + (A45 * psi101) + (A46 * psi110) + (A47 * psi111);

        // Row 5 -> result_k101 (|101>)
        result.coeffRef(k101) =
            (A50 * psi000) + (A51 * psi001) + (A52 * psi010) + (A53 * psi011) +
            (A54 * psi100) + (A55 * psi101) + (A56 * psi110) + (A57 * psi111);

        // Row 6 -> result_k110 (|110>)
        result.coeffRef(k110) =
            (A60 * psi000) + (A61 * psi001) + (A62 * psi010) + (A63 * psi011) +
            (A64 * psi100) + (A65 * psi101) + (A66 * psi110) + (A67 * psi111);

        // Row 7 -> result_k111 (|111>)
        result.coeffRef(k111) =
            (A70 * psi000) + (A71 * psi001) + (A72 * psi010) + (A73 * psi011) +
            (A74 * psi100) + (A75 * psi101) + (A76 * psi110) + (A77 * psi111);
    }

    return result;
}

/**
 * \brief Applies the multi-qubit gate \a A to the part \a target of the
 * multi-partite state vector \a state
 *
 * \param state Eigen expression
 * \param A Eigen expression
 * \param target Subsystem indexes where the gate \a A is applied
 * \param n Number of qubits
 * \return Gate \a A applied to the part \a target of \a state
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] qpp::expr_t<Derived1>
apply_psi_kq(const Eigen::MatrixBase<Derived1>& state,
             const Eigen::MatrixBase<Derived2>& A,
             const std::vector<idx>& target, idx n) {
    // Type Aliases
    using Scalar = typename Derived1::Scalar;
    using EigenVector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

    // Setup
    const idx k = static_cast<idx>(target.size());
    const idx dim = 1ULL << k;
    const idx outer_dim = 1ULL << (n - k);

    // Input Validation
#ifndef NDEBUG
    // D is the dimension of the state (2^n)
    const idx D = 1ULL << n;

    // Check Target Qubit Indices are valid
    for (idx t : target) {
        assert(t < n && "Target qubit index must be less than n");
    }

    // Check Gate Dimension: A must be a 2^k x 2^k matrix
    assert(static_cast<idx>(A.rows()) == dim &&
           static_cast<idx>(A.cols()) == dim &&
           "Gate A must be a 2^k x 2^k matrix, where k is the number of target "
           "qubits");

    // Check State Dimension: state must be 2^n x 1 vector
    assert(static_cast<idx>(state.size()) == D && state.cols() == 1 &&
           "State must be a 2^n x 1 column vector (ket)");

    // Check Target Qubit Indices are distinct
    if (k > 1) {
        std::vector<idx> sorted_target = target;
        std::sort(sorted_target.begin(), sorted_target.end());
        assert(std::adjacent_find(sorted_target.begin(), sorted_target.end()) ==
                   sorted_target.end() &&
               "Target qubit indices must be distinct");
    }
#endif

    // Create the output state vector as a copy of the input.
    // This allows in-place updates via parallel access to non-overlapping
    // blocks.
    qpp::expr_t<Derived1> result = state;

    // Pre-calculation of Inner Indices (Gate Index Mapping)
    // Maps local index r to its global index contribution.
    std::vector<idx> inner_idx(dim);
    for (idx r = 0; r < dim; ++r) {
        idx index_r = 0;
        for (idx j = 0; j < k; ++j) {
            const idx target_qubit_pos = target[k - 1 - j]; // Permutation fix
            if ((r >> j) & 1) {
                // BIG-ENDIAN contribution
                index_r += (1ULL << (n - 1 - target_qubit_pos));
            }
        }
        inner_idx[r] = index_r;
    }

    // Pre-calculation of Outer Indices (Spectator Block Base Index)
    // Maps spectator state index m to its global base index (i_base).
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
            if ((m >> j) & 1) {
                i_base += (1ULL << (n - 1 - spectator_qubits[j]));
            }
        }
        outer_idx[m] = i_base;
    }

    // Main Parallel Loop (Iterate over spectator blocks m)
#ifdef QPP_OPENMP
// NOLINTNEXTLINE
#pragma omp parallel for
#endif // QPP_OPENMP
    for (idx m = 0; m < outer_dim; ++m) {
        const idx i_base = outer_idx[m];

        // Input Block Caching (The Optimization)
        // Read the scattered input block (dim elements) into a contiguous
        // temporary vector. This ensures the GEMV in step 2 benefits from
        // perfect cache locality.
        EigenVector input_block(dim);
        for (idx c = 0; c < dim; ++c) {
            input_block(c) = state(i_base + inner_idx[c]);
        }

        // Core Transformation (Block GEMV)
        // The output of the matrix multiplication is stored in new_block
        EigenVector new_block = A * input_block;

        // Write Output Block (Scattered Write)
        // Write the contiguous result block (new_block) back to the scattered
        // global indices.
        for (idx r = 0; r < dim; ++r) {
            const idx i_out = i_base + inner_idx[r];
            result(i_out) = new_block(r);
        }
    }

    return result;
}

/**
 * \brief Applies the 1-qubit gate \a A to the qubit \a i of the
 * multi-partite density matrix \a state
 *
 * \param state Eigen expression
 * \param A Eigen expression (2x2 matrix)
 * \param i Subsystem index where the gate \a A is applied
 * \param n Number of qubits
 * \return Gate \a A applied to the qubit \a i of \a state
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] expr_t<Derived1>
apply_rho_1q(const Eigen::MatrixBase<Derived1>& state,
             const Eigen::MatrixBase<Derived2>& A, idx i, idx n) {
    using Scalar = typename Derived1::Scalar;
    using Matrix2 = Eigen::Matrix2<Scalar>;
    const idx D = 1ULL << n; // total Hilbert space dimension

    // Input Validation
#ifndef NDEBUG
    assert(i < n && "Target qubit index i must be less than n");
    assert(static_cast<idx>(state.rows()) == D &&
           static_cast<idx>(state.cols()) == D &&
           "State must be a square matrix sized 2^n x 2^n");
    assert(A.rows() == 2 && A.cols() == 2 && "Gate A must be a 2x2 matrix");
#endif

    expr_t<Derived1> result = state; // deep copy

    // Bit position for qubit i (p_i) and the stride (s_i = 2^p_i)
    const idx p_i = n - 1 - i;
    const idx s_i = 1ULL << p_i;

    // Constants for index reconstruction
    const idx D_spec = D / 2; // Total number of spectator states (2^(n-1))
    const idx p_i_plus_1 = p_i + 1;
    const idx low_mask = s_i - 1; // Mask for bits below p_i

    const Matrix2 A_block = A;
    const Matrix2 A_dagger = A_block.adjoint();

    // The primary optimization: Iterate over D_spec * D_spec blocks directly,
    // eliminating conditional checks inside the loops.

    const idx total_iterations = D_spec * D_spec;
#ifdef QPP_OPENMP
// NOLINTNEXTLINE
#pragma omp parallel for
#endif // QPP_OPENMP
    for (idx i = 0; i < total_iterations; ++i) {
        const idx s = i / D_spec;
        const idx s_prime = i % D_spec;
        // Calculate Row Indices (r0, r1) from spectator state 's'
        // r0 = index for |s, 0> (inserting '0' at bit p_i)
        const idx s_high_r = s >> p_i;
        const idx s_low_r = s & low_mask;
        const idx r0 = (s_high_r << p_i_plus_1) | s_low_r;
        const idx r1 = r0 + s_i; // r1 = index for |s, 1>

        // Calculate Column Indices (c0, c1) from spectator state
        // 's_prime' c0 = index for |s', 0> (inserting '0' at bit p_i)
        const idx s_prime_high_c = s_prime >> p_i;
        const idx s_prime_low_c = s_prime & low_mask;
        const idx c0 = (s_prime_high_c << p_i_plus_1) | s_prime_low_c;
        const idx c1 = c0 + s_i; // c1 = index for |s', 1>

        // Extract the current 2x2 block rho[r0:r1, c0:c1]
        // We use Matrix2cd (fixed size) for maximum optimization by Eigen.
        Matrix2 rho_block;
        rho_block(0, 0) = state.coeff(r0, c0);
        rho_block(0, 1) = state.coeff(r0, c1);
        rho_block(1, 0) = state.coeff(r1, c0);
        rho_block(1, 1) = state.coeff(r1, c1);

        // Apply the transformation: A * rho_block * A^\dagger
        // Eigen will fully unroll and vectorize this small matrix
        // multiplication.
        const Matrix2 result_block = A_block * rho_block * A_dagger;

        // 5. Write the transformed 2x2 block back
        result.coeffRef(r0, c0) = result_block(0, 0);
        result.coeffRef(r0, c1) = result_block(0, 1);
        result.coeffRef(r1, c0) = result_block(1, 0);
        result.coeffRef(r1, c1) = result_block(1, 1);
    }

    return result;
}

/**
 * \brief Applies the 2-qubit gate \a A to the qubits \a i and \a j of the
 * multi-partite density matrix \a state
 *
 * \param state Eigen expression
 * \param A Eigen expression (4x4 matrix)
 * \param i Subsystem index where the gate \a A is applied
 * \param j Subsystem index where the gate \a A is applied
 * \param n Number of qubits
 * \return Gate \a A applied to the qubits \a i and \a j of \a state
 */
template <typename Derived1, typename Derived2>
expr_t<Derived1> apply_rho_2q(const Eigen::MatrixBase<Derived1>& state,
                              const Eigen::MatrixBase<Derived2>& A, idx i,
                              idx j, idx n) {

    const idx D = static_cast<idx>(state.rows());

    // Input Validation
#ifndef NDEBUG
    const idx D_expected = (1ULL << n);
    assert(i < n && j < n && i != j &&
           "Target qubit indices i and j must be distinct and less than n");
    assert(D == D_expected && D == static_cast<idx>(state.cols()) &&
           "State must be a square matrix sized 2^n x 2^n");
    assert(A.rows() == 4 && A.cols() == 4 && "Gate A must be a 4x4 matrix");
#endif

    // Index Conversion (Big Endian to Little Endian)
    // Convert user's MSB-first indices (0 is MSB) to physical LSB-first indices
    // (0 is LSB)
    const idx i_phys = n - i - 1;
    const idx j_phys = n - j - 1;

    // Type and Dimension Setup
    using Scalar = typename Derived1::Scalar;
    using ComputeMatrixType = Eigen::Matrix<Scalar, -1, -1>;
    using ComputeBlockType = Eigen::Matrix<Scalar, 4, 4>;

    const ComputeMatrixType rho_cd = state.template cast<Scalar>().eval();
    const ComputeBlockType U = A.template cast<Scalar>().eval();
    const ComputeBlockType U_adj = U.adjoint();

    // Size of the N-2 qubit subsystem
    const idx D_rest = (n >= 2) ? (1ULL << (n - 2)) : 1;

    // Powers of 2 are calculated using the physical (LSB-first) indices.
    const idx P_i = 1ULL << i_phys;
    const idx P_j = 1ULL << j_phys;

    ComputeMatrixType rho_prime_cd = ComputeMatrixType::Zero(D, D);

    // Block Iteration (Parallelized)
    // Parallelize the outermost loop (row-blocks). This is safe as each thread
    // writes to a non-overlapping region of rho_prime_cd.
#ifdef QPP_OPENMP
// NOLINTNEXTLINE
#pragma omp parallel for default(none) shared(                                 \
        D_rest, n, P_i, P_j, rho_cd, U, U_adj, rho_prime_cd, i_phys, j_phys)
#endif // QPP_OPENMP
    // Loop over row-blocks
    for (idx r = 0; r < D_rest; ++r) {

        // Variables private to this thread
        idx r_base_row = 0;
        idx current_r = r;

        // Calculate the 'base' index for the row block (rest qubits)
        for (idx q = 0; q < n;
             ++q) { // q is the physical index (0=LSB, n-1=MSB)
            // Skip the physical qubits we are acting on (i_phys and j_phys).
            if (q != i_phys && q != j_phys) {
                if (current_r & 1) {
                    r_base_row |= (1ULL << q);
                }
                current_r >>= 1;
            }
        }

        // Define the 4 *row* indices, aligned with the gate A's basis |q_i
        // q_j> P_i corresponds to the MSB of the 4x4 block, P_j to the LSB.
        const idx vec_k_row[4] = {
            r_base_row,            // |00>
            r_base_row + P_j,      // |01> (LSB set)
            r_base_row + P_i,      // |10> (MSB set)
            r_base_row + P_i + P_j // |11> (Both set)
        };

        // Inner loop (col-blocks)
        for (idx c = 0; c < D_rest; ++c) {
            // Variables private to this inner iteration
            idx r_base_col = 0;
            idx current_c = c;
            ComputeBlockType M;
            ComputeBlockType M_prime;
            idx vec_k_col[4];

            // 1b. Calculate the 'base' index for the column block (rest qubits)
            for (idx q = 0; q < n; ++q) {
                if (q != i_phys && q != j_phys) {
                    if (current_c & 1) {
                        r_base_col |= (1ULL << q);
                    }
                    current_c >>= 1;
                }
            }

            // Define the 4 *column* indices
            vec_k_col[0] = r_base_col;
            vec_k_col[1] = r_base_col + P_j;
            vec_k_col[2] = r_base_col + P_i;
            vec_k_col[3] = r_base_col + P_i + P_j;

            // Extract the 4x4 block M from the *input* state
            for (int row = 0; row < 4; ++row) {
                for (int col = 0; col < 4; ++col) {
                    M(row, col) = rho_cd(vec_k_row[row], vec_k_col[col]);
                }
            }

            // Apply the gate: M' = U * M * U_adj (Optimized 4x4
            // multiplication)
            M_prime = (U * M * U_adj).eval();

            // 5. Insert the resulting 4x4 block M_prime back into the result
            for (int row = 0; row < 4; ++row) {
                for (int col = 0; col < 4; ++col) {
                    rho_prime_cd(vec_k_row[row], vec_k_col[col]) =
                        M_prime(row, col);
                }
            }
        } // end col loop
    } // end row loop (parallelized)

    // Cast back to the original's scalar type
    return rho_prime_cd.template cast<typename Derived1::Scalar>();
}

/**
 * \brief Applies the 3-qubit gate \a A to the qubit \a i, \a j, and \a k of the
 * multi-partite density matrix \a state
 *
 * \param state Eigen expression
 * \param A Eigen expression (8x8 matrix)
 * \param i Subsystem index where the gate \a A is applied
 * \param j Subsystem index where the gate \a A is applied
 * \param k Subsystem index where the gate \a A is applied
 * \param n Number of qubits
 * \return Gate \a A applied to the qubit \a target of \a state
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] expr_t<Derived1>
apply_rho_3q(const Eigen::MatrixBase<Derived1>& state,
             const Eigen::MatrixBase<Derived2>& A, idx i, idx j, idx k, idx n) {

    const idx D = static_cast<idx>(state.rows());

    // Input Validation
#ifndef NDEBUG
    const idx D_expected = (1ULL << n);
    assert(i < n && j < n && k < n && i != j && i != k && j != k &&
           "Target qubit indices i, j, and k must be distinct and less than n");
    assert(D == D_expected && D == static_cast<idx>(state.cols()) &&
           "State must be a square matrix of size 2^n x 2^n");
    assert(A.rows() == 8 && A.cols() == 8 && "Gate A must be an 8x8 matrix");
    assert(n >= 3 && "Need at least 3 qubits for a 3-qubit gate");
#endif

    // Index Conversion (Big Endian to Little Endian
    // Convert user's MSB-first indices (0 is MSB) to physical LSB-first indices
    // (0 is LSB)
    const idx i_phys = n - i - 1;
    const idx j_phys = n - j - 1;
    const idx k_phys = n - k - 1;

    // Type and Dimension Setup
    using Scalar = typename Derived1::Scalar;
    using ComputeMatrixType = Eigen::Matrix<Scalar, -1, -1>;
    // The block size is 8x8 for a 3-qubit gate
    using ComputeBlockType = Eigen::Matrix<Scalar, 8, 8>;

    const ComputeMatrixType rho_cd = state.template cast<Scalar>().eval();
    const ComputeBlockType U = A.template cast<Scalar>().eval();
    const ComputeBlockType U_adj = U.adjoint();

    // Size of the N-3 qubit subsystem (the "rest" of the system)
    const idx D_rest = (n >= 3) ? (1ULL << (n - 3)) : 1;

    // Powers of 2 corresponding to the physical qubit indices
    const idx P_i = 1ULL << i_phys;
    const idx P_j = 1ULL << j_phys;
    const idx P_k = 1ULL << k_phys;

    ComputeMatrixType rho_prime_cd = ComputeMatrixType::Zero(D, D);

    // Block Iteration (Parallelized)
    // Loop over row-blocks (r) and column-blocks (c) corresponding to the rest
    // of the system (D_rest x D_rest blocks). Each block is 8x8.
#ifdef QPP_OPENMP
// NOLINTNEXTLINE
#pragma omp parallel for default(none)                                         \
    shared(D_rest, n, P_i, P_j, P_k, rho_cd, U, U_adj, rho_prime_cd, i_phys,   \
               j_phys, k_phys)
#endif // QPP_OPENMP
    // Loop over row-blocks (rest index)
    for (idx r = 0; r < D_rest; ++r) {

        // Variables private to this thread
        idx r_base_row = 0;
        idx current_r = r;

        // Calculate the 'base' index for the row block (rest qubits)
        // This index represents the |rest> part of the basis state |i j k rest>
        for (idx q = 0; q < n;
             ++q) { // q is the physical index (0=LSB, n-1=MSB)
            // Skip the physical qubits we are acting on (i_phys, j_phys,
            // k_phys)
            if (q != i_phys && q != j_phys && q != k_phys) {
                if (current_r & 1) {
                    r_base_row |= (1ULL << q);
                }
                current_r >>= 1;
            }
        }

        // Define the 8 *row* indices, aligned with the gate A's basis |q_i
        // q_j q_k> Assuming A's basis order is |000> to |111> with i=MSB,
        // j=Mid, k=LSB.
        const idx vec_k_row[8] = {
            r_base_row,                  // |000>
            r_base_row + P_k,            // |001>
            r_base_row + P_j,            // |010>
            r_base_row + P_j + P_k,      // |011>
            r_base_row + P_i,            // |100>
            r_base_row + P_i + P_k,      // |101>
            r_base_row + P_i + P_j,      // |110>
            r_base_row + P_i + P_j + P_k // |111>
        };

        // Inner loop (col-blocks)
        for (idx c = 0; c < D_rest; ++c) {
            // Variables private to this inner iteration
            idx r_base_col = 0;
            idx current_c = c;
            ComputeBlockType M;
            ComputeBlockType M_prime;
            idx vec_k_col[8];

            // Calculate the 'base' index for the column block (rest qubits)
            for (idx q = 0; q < n; ++q) {
                if (q != i_phys && q != j_phys && q != k_phys) {
                    if (current_c & 1) {
                        r_base_col |= (1ULL << q);
                    }
                    current_c >>= 1;
                }
            }

            // Define the 8 *column* indices (using the same order)
            vec_k_col[0] = r_base_col;
            vec_k_col[1] = r_base_col + P_k;
            vec_k_col[2] = r_base_col + P_j;
            vec_k_col[3] = r_base_col + P_j + P_k;
            vec_k_col[4] = r_base_col + P_i;
            vec_k_col[5] = r_base_col + P_i + P_k;
            vec_k_col[6] = r_base_col + P_i + P_j;
            vec_k_col[7] = r_base_col + P_i + P_j + P_k;

            // Extract the 8x8 block M from the *input* state
            for (int row = 0; row < 8; ++row) {
                for (int col = 0; col < 8; ++col) {
                    M(row, col) = rho_cd(vec_k_row[row], vec_k_col[col]);
                }
            }

            // Apply the gate: M' = U * M * U_adj (Optimized 8x8
            // multiplication)
            M_prime = (U * M * U_adj).eval();

            // Insert the resulting 8x8 block M_prime back into the result
            for (int row = 0; row < 8; ++row) {
                for (int col = 0; col < 8; ++col) {
                    rho_prime_cd(vec_k_row[row], vec_k_col[col]) =
                        M_prime(row, col);
                }
            }
        } // end col loop
    } // end row loop (parallelized)

    // Return
    // Cast back to the original's scalar type
    return rho_prime_cd.template cast<typename Derived1::Scalar>();
}

/**
 * \brief Applies the multi-qubit gate \a A to the part \a target of the
 * multi-partite density matrix \a state
 *
 * \param state Eigen expression
 * \param A Eigen expression
 * \param target Subsystem indexes where the gate \a A is applied
 * \return Gate \a A applied to the part \a target of \a state
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] expr_t<Derived1>
apply_rho_kq(const Eigen::MatrixBase<Derived1>& state,
             const Eigen::MatrixBase<Derived2>& A,
             const std::vector<idx>& target, idx n) {
    using Scalar = typename Derived1::Scalar;
    using ComputeMatrixType = Eigen::Matrix<Scalar, -1, -1>;

    const idx D = static_cast<idx>(state.rows());
    const idx k = target.size(); // Gate qubits

    // Calculate block dimension D_k = 2^k
    const idx D_k = (k == 0) ? 1 : (1ULL << k);

    // Input Validation
#ifndef NDEBUG
    // Check Target Qubit Indices are valid
    for (idx t : target) {
        assert(t < n && "Target qubit index must be less than n");
    }

    // Check Gate Dimension: A must be a 2^k x 2^k matrix
    assert(static_cast<idx>(A.rows()) == D_k &&
           static_cast<idx>(A.cols()) == D_k &&
           "Gate A must be a 2^k x 2^k matrix, where k is the number of target "
           "qubits");

    // Check State Dimension: state must be 2^n x 2^n matrix
    assert(static_cast<idx>(state.rows()) == D &&
           static_cast<idx>(state.cols()) == D &&
           "State must be a 2^n x 2^n matrix");

    // Validate target indices are distinct and in range
    if (k > 1) {
        std::vector<idx> sorted_target = target;
        std::sort(sorted_target.begin(), sorted_target.end());
        assert(std::unique(sorted_target.begin(), sorted_target.end()) ==
                   sorted_target.end() &&
               "Target qubit indices must be distinct");
    }
#endif

    // Pre-Calculations for Index Mapping

    // Convert user's logical (MSB-first) target indices to physical (LSB-first)
    // indices and determine their powers of 2.
    // P_gate_basis[l] is the power of 2 for the physical index of target[l].
    std::vector<idx> target_phys(k);
    std::vector<idx> P_gate_basis(k);
    for (idx l = 0; l < k; ++l) {
        // Physical index: 0 is LSB, n-1 is MSB
        target_phys[l] = n - target[l] - 1;
        P_gate_basis[l] = 1ULL << target_phys[l];
    }

    // Determine the physical indices that are *not* in the target set (the
    // "rest")
    std::vector<idx> rest_phys;
    // Iterate over all physical indices (0 to n-1)
    for (idx q_phys = 0; q_phys < n; ++q_phys) {
        bool is_target = false;
        for (const auto& t_phys : target_phys) {
            if (q_phys == t_phys) {
                is_target = true;
                break;
            }
        }
        if (!is_target) {
            rest_phys.push_back(q_phys);
        }
    }

    // Size of the N-k qubit subsystem
    const idx D_rest = (n >= k) ? (1ULL << (n - k)) : 1;

    // Matrix Setup
    const ComputeMatrixType rho_cd = state.template cast<Scalar>().eval();
    const ComputeMatrixType U = A.template cast<Scalar>().eval();
    const ComputeMatrixType U_adj = U.adjoint();

    ComputeMatrixType rho_prime_cd = ComputeMatrixType::Zero(D, D);

    // Block Iteration (Parallelized)
    // Loop over row-blocks (r) and column-blocks (c) corresponding to the
    // rest of the system (D_rest x D_rest blocks). Each block is D_k x D_k.
#ifdef QPP_OPENMP
// NOLINTNEXTLINE
#pragma omp parallel for default(none)                                         \
    shared(D_rest, n, k, rest_phys, target, target_phys, P_gate_basis, rho_cd, \
               U, U_adj, rho_prime_cd, D_k)
#endif // QPP_OPENMP
    // Loop over row-blocks (rest index)
    for (idx r = 0; r < D_rest; ++r) {

        // Variables private to this thread
        idx r_base_row = 0;
        idx current_r = r;

        // Calculate the 'base' index for the row block (rest qubits)
        // This index represents the |rest> part of the basis state
        for (const auto& q_phys : rest_phys) {
            if (current_r & 1) {
                r_base_row |= (1ULL << q_phys);
            }
            current_r >>= 1;
        }

        // Pre-calculate the D_k *row* indices for the block.
        // The index m is the linear index in the gate's basis |q_t0 ...
        // q_{t(k-1)}>, where bit l of m determines the state of target[k-1-l]
        // (MSB/LSB flip)
        std::vector<idx> vec_k_row(D_k);
        for (idx m = 0; m < D_k; ++m) { // m iterates over the 2^k basis states
            idx target_component = 0;
            for (idx l = 0; l < k; ++l) {
                // l=0 is MSB of the gate (target[0]), l=k-1 is LSB
                // (target[k-1]) The bit corresponding to target[l] is bit
                // (k-1-l) of m.
                idx bit_pos = k - 1 - l;
                if ((m >> bit_pos) & 1) {
                    target_component += P_gate_basis[l];
                }
            }
            vec_k_row[m] = r_base_row + target_component;
        }

        // Inner loop (col-blocks)
        for (idx c = 0; c < D_rest; ++c) {
            // Variables private to this inner iteration
            idx r_base_col = 0;
            idx current_c = c;
            ComputeMatrixType M(D_k, D_k);
            ComputeMatrixType M_prime;
            std::vector<idx> vec_k_col(D_k);

            // Calculate the 'base' index for the column block (rest qubits)
            for (const auto& q_phys : rest_phys) {
                if (current_c & 1) {
                    r_base_col |= (1ULL << q_phys);
                }
                current_c >>= 1;
            }

            // Pre-calculate the D_k *column* indices (using the same order)
            for (idx m = 0; m < D_k; ++m) {
                idx target_component = 0;
                for (idx l = 0; l < k; ++l) {
                    idx bit_pos = k - 1 - l;
                    if ((m >> bit_pos) & 1) {
                        target_component += P_gate_basis[l];
                    }
                }
                vec_k_col[m] = r_base_col + target_component;
            }

            // Extract the D_k x D_k block M from the *input* state
            for (idx row = 0; row < D_k; ++row) {
                for (idx col = 0; col < D_k; ++col) {
                    M(row, col) = rho_cd(vec_k_row[row], vec_k_col[col]);
                }
            }

            // Apply the gate: M' = U * M * U_adj
            M_prime = (U * M * U_adj).eval();

            // Insert the resulting block M_prime back into the result
            for (idx row = 0; row < D_k; ++row) {
                for (idx col = 0; col < D_k; ++col) {
                    rho_prime_cd(vec_k_row[row], vec_k_col[col]) =
                        M_prime(row, col);
                }
            }
        } // end col loop
    } // end row loop (parallelized)

    // Return
    return rho_prime_cd.template cast<typename Derived1::Scalar>();
}

/**
 * \brief Applies the controlled 1-qubit gate \a A to the qubit \a i of the
 * multi-partite state vector \a state
 *
 * The gate is applied only if the state of the control qubits \a ctrl matches
 * the pattern defined by \a shift (0 for |1> control, 1 for |0> control)
 *
 * \param state Eigen expression
 * \param A Eigen expression (2x2 matrix)
 * \param i Target subsystem index
 * \param ctrl Vector of control qubit indices
 * \param shift Vector of control values (0: positive/|1> control, 1:
 * negative/|0> control)
 * \param n Number of qubits
 * \return Controlled gate \a A applied to the qubit \a i of \a state
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] expr_t<Derived1>
apply_ctrl_psi_1q(const Eigen::MatrixBase<Derived1>& state,
                  const Eigen::MatrixBase<Derived2>& A,
                  const std::vector<idx>& ctrl, idx i,
                  const std::vector<idx>& shift, idx n) {
    // Type and Dimension Setup
    using Scalar = typename Derived1::Scalar;
    const idx D = 1ULL << n; // Total size of the state vector (2^n)

    // Input Validation
#ifndef NDEBUG
    assert(i < n && "Target qubit index i must be less than n");
    assert(static_cast<idx>(state.size()) == D &&
           "State vector size must be 2^n");
    assert(A.rows() == 2 && A.cols() == 2 && "Gate A must be a 2x2 matrix");
    assert(ctrl.size() == shift.size() &&
           "ctrl and shift vectors must have the same size");

    // Validate control qubits and shift values
    for (idx c_idx = 0; c_idx < ctrl.size(); ++c_idx) {
        const idx c = ctrl[c_idx];
        assert(c < n && "Control qubit index c must be less than n");
        assert(c != i && "Target qubit i cannot also be a control qubit");
        assert(
            (shift[c_idx] == 0 || shift[c_idx] == 1) &&
            "Shift value must be 0 (positive control) or 1 (negative control)");
    }
#endif

    // A deep copy of the input state is required. If the control condition
    // is not met, the state remains unchanged for that amplitude pair.
    expr_t<Derived1> result = state;

    // --- Bit-Shift Constants (BIG-ENDIAN) ---
    // The effective bit position 'j' corresponding to qubit 'i' in big-endian
    // is (n - 1 - i).
    const idx j = n - 1 - i;

    // 'step' is 2^j. Index difference between |...0...> and |...1...> at bit j.
    const idx step = 1ULL << j;

    // 'jump' is 2^(j+1). Size of the block that repeats.
    const idx jump = 1ULL << (j + 1);

    // --- Control Mask and Expected Pattern Calculation ---
    // Mask for controls that MUST be |1> (shift == 0, Positive Control)
    idx expected_pattern_for_ones = 0;
    // Mask for controls that MUST be |0> (shift == 1, Negative Control)
    idx expected_zero_mask = 0;

    // We iterate over the list of control qubits (ctrl)
    for (idx c_idx = 0; c_idx < ctrl.size(); ++c_idx) {
        const idx c = ctrl[c_idx];
        const idx j_c =
            n - 1 - c; // Big-endian bit position for control qubit c
        const idx bit = (1ULL << j_c); // Pre-calculate the bit to be set

        // If shift[c_idx] == 0: Positive Control (requires |1>)
        if (shift[c_idx] == 0) {
            expected_pattern_for_ones |= bit;
        }
        // If shift[c_idx] == 1: Negative Control (requires |0>)
        else {
            // Since validation ensures shift is 0 or 1, this is the shift == 1
            // case
            expected_zero_mask |= bit;
        }
    }
    // --- Gate Matrix Elements ---
    const Scalar a00 = A.coeff(0, 0);
    const Scalar a01 = A.coeff(0, 1);
    const Scalar a10 = A.coeff(1, 0);
    const Scalar a11 = A.coeff(1, 1);

    // Pair-wise Amplitude Transformation
#ifdef QPP_OPENMP
// NOLINTNEXTLINE
#pragma omp parallel for
#endif // QPP_OPENMP
    for (idx L = 0; L < D; L += jump) {
        // The inner loop (R) iterates over the lower part of the block.
        for (idx R = 0; R < step; ++R) {
            // k0 is the index where the j-th bit (target qubit) is 0: L + R
            const idx k0 = L + R;
            // k1 is the index where the j-th bit (target qubit) is 1: k0 + step
            [[maybe_unused]] const idx k1 = k0 + step;

            // CONTROL CHECK: Apply gate A only if ALL conditions are met:
            // 1. All positive control bits are 1: (k0 &
            // expected_pattern_for_ones) == expected_pattern_for_ones
            // 2. All negative control bits are 0: (k0 & expected_zero_mask) ==
            // 0
            if (((k0 & expected_pattern_for_ones) ==
                 expected_pattern_for_ones) &&
                ((k0 & expected_zero_mask) == 0)) {

                // Fetch the ORIGINAL amplitudes from the INPUT 'state'
                const Scalar& psi_k0 = state.coeff(k0);
                const Scalar& psi_k1 = state.coeff(k1);

                // Apply the 2x2 matrix-vector multiplication and write to
                // 'result'
                result.coeffRef(k0) = (a00 * psi_k0) + (a01 * psi_k1);
                result.coeffRef(k1) = (a10 * psi_k0) + (a11 * psi_k1);
            }
        }
    }

    return result;
}

/**
 * \brief Applies the controlled 2-qubit gate \a A to the qubits \a i and \a j
 * of the multi-partite state vector \a state
 *
 * The gate is applied only if the state of the control qubits \a ctrl matches
 * the pattern defined by \a shift (0 for |1> control, 1 for |0> control)
 *
 * \param state Eigen expression
 * \param A Eigen expression (4x4 matrix)
 * \param ctrl Vector of control qubit indices
 * \param i Target subsystem index
 * \param j Target subsystem index
 * \param shift Vector of control values (0: positive/|1> control, 1:
 * negative/|0> control)
 * \param n Number of qubits
 * \return Controlled gate \a A applied to the qubits \a i and \a j of \a state
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] expr_t<Derived1>
apply_ctrl_psi_2q(const Eigen::MatrixBase<Derived1>& state,
                  const Eigen::MatrixBase<Derived2>& A,
                  const std::vector<idx>& ctrl, idx i, idx j,
                  const std::vector<idx>& shift, idx n) {
    using Scalar = typename Derived1::Scalar;
    const idx D = 1ULL << n;

#ifndef NDEBUG
    assert(i < n && j < n && i != j &&
           "Target qubit indices i and j must be distinct and less than n");
    assert(static_cast<idx>(state.size()) == D &&
           "State vector size must be 2^n");
    assert(A.rows() == 4 && A.cols() == 4 && "Gate A must be a 4x4 matrix");
    assert(ctrl.size() == shift.size() &&
           "ctrl and shift vectors must have the same size");
    for (idx c_idx = 0; c_idx < ctrl.size(); ++c_idx) {
        const idx c = ctrl[c_idx];
        assert(c < n && "Control qubit index c must be less than n");
        assert(c != i && c != j &&
               "Target qubits i and j cannot also be control qubits");
        assert(
            (shift[c_idx] == 0 || shift[c_idx] == 1) &&
            "Shift value must be 0 (positive control) or 1 (negative control)");
    }
#endif

    expr_t<Derived1> result =
        state; // deep copy; we'll overwrite only where control holds

    // Big-endian bit positions
    const idx p_i = n - 1 - i;
    const idx p_j = n - 1 - j;
    const idx s_i = (1ULL << p_i);
    const idx s_j = (1ULL << p_j);

    // build control mask/value using big-endian positions
    idx control_mask = 0;
    idx control_value = 0;
    for (idx c_idx = 0; c_idx < ctrl.size(); ++c_idx) {
        const idx c = ctrl[c_idx];
        const idx p_c = n - 1 - c; // big-endian position
        const idx bit = (1ULL << p_c);
        control_mask |= bit;
        if (shift[c_idx] == 0) { // require |1>
            control_value |= bit;
        }
        // else require |0> -> leave bit 0 in control_value
    }

    // extract A entries
    const Scalar a00 = A.coeff(0, 0), a01 = A.coeff(0, 1), a02 = A.coeff(0, 2),
                 a03 = A.coeff(0, 3);
    const Scalar a10 = A.coeff(1, 0), a11 = A.coeff(1, 1), a12 = A.coeff(1, 2),
                 a13 = A.coeff(1, 3);
    const Scalar a20 = A.coeff(2, 0), a21 = A.coeff(2, 1), a22 = A.coeff(2, 2),
                 a23 = A.coeff(2, 3);
    const Scalar a30 = A.coeff(3, 0), a31 = A.coeff(3, 1), a32 = A.coeff(3, 2),
                 a33 = A.coeff(3, 3);

    // Iterate over every base index where target bits are 0.
    // That visits each 4-amplitude block exactly once.
#ifdef QPP_OPENMP
// NOLINTNEXTLINE
#pragma omp parallel for
#endif // QPP_OPENMP
    for (idx k00 = 0; k00 < D; ++k00) {
        // skip indices that do not have both target bits == 0
        if ((k00 & s_i) || (k00 & s_j)) {
            continue;
        }

        // control check (control bits of k00 must match control_value under
        // control_mask)
        if ((k00 & control_mask) != control_value) {
            continue;
        }

        // indices for the 2-qubit subspace: |0_i 0_j>, |0_i 1_j>, |1_i 0_j>,
        // |1_i 1_j>
        const idx k01 = k00 | s_j;       // flip j
        const idx k10 = k00 | s_i;       // flip i
        const idx k11 = k00 | s_i | s_j; // flip both

        // read original amplitudes from 'state'
        const Scalar psi00 = state.coeff(k00);
        const Scalar psi01 = state.coeff(k01);
        const Scalar psi10 = state.coeff(k10);
        const Scalar psi11 = state.coeff(k11);

        // write transformed amplitudes into 'result' (based on A *
        // [psi00,psi01,psi10,psi11]^T)
        result.coeffRef(k00) =
            (a00 * psi00) + (a01 * psi01) + (a02 * psi10) + (a03 * psi11);
        result.coeffRef(k01) =
            (a10 * psi00) + (a11 * psi01) + (a12 * psi10) + (a13 * psi11);
        result.coeffRef(k10) =
            (a20 * psi00) + (a21 * psi01) + (a22 * psi10) + (a23 * psi11);
        result.coeffRef(k11) =
            (a30 * psi00) + (a31 * psi01) + (a32 * psi10) + (a33 * psi11);
    }

    return result;
}

/**
 * \brief Applies the multi-controlled qubit gate \a A to the part \a target of
 * the multi-partite state vector \a state
 *
 * The gate is applied only if the state of the control qubits \a ctrl matches
 * the pattern defined by \a shift (0 for |1> control, 1 for |0> control)
 *
 * \param state Eigen expression
 * \param A Eigen expression (4x4 matrix)
 * \param ctrl Vector of control qubit indices
 * \param target Subsystem indexes where the gate \a A is applied
 * \param shift Vector of control values (0: positive/|1> control, 1:
 * negative/|0> control)
 * \param n Number of qubits
 * \return Controlled gate \a A applied to the part \a target of \a state
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] qpp::expr_t<Derived1>
apply_ctrl_psi_kq(const Eigen::MatrixBase<Derived1>& state,
                  const Eigen::MatrixBase<Derived2>& A,
                  const std::vector<idx>& ctrl, const std::vector<idx>& target,
                  const std::vector<idx>& shift, idx n) {
    // Type Aliases
    using Scalar = typename Derived1::Scalar;
    using EigenVector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

    // Setup
    const idx k = static_cast<idx>(target.size()); // Number of target qubits
    const idx dim = 1ULL << k;             // Target subspace dimension (2^k)
    const idx outer_dim = 1ULL << (n - k); // Number of spectator blocks

    // Input Validation
#ifndef NDEBUG
    // D is the dimension of the state (2^n)
    const idx D = 1ULL << n;

    // Check Target Qubit Indices are valid
    for (idx t : target) {
        assert(t < n && "Target qubit index must be less than n");
    }

    // Check Gate Dimension: A must be a 2^k x 2^k matrix
    assert(static_cast<idx>(A.rows()) == dim &&
           static_cast<idx>(A.cols()) == dim &&
           "Gate A must be a 2^k x 2^k matrix, where k is the number of target "
           "qubits");

    // Check State Dimension: state must be 2^n x 1 vector
    assert(static_cast<idx>(state.size()) == D && state.cols() == 1 &&
           "State must be a 2^n x 1 column vector (ket)");

    // Check Target Qubit Indices are distinct
    if (k > 1) {
        std::vector<idx> sorted_target = target;
        std::sort(sorted_target.begin(), sorted_target.end());
        assert(std::adjacent_find(sorted_target.begin(), sorted_target.end()) ==
                   sorted_target.end() &&
               "Target qubit indices must be distinct");
    }

    // Check Control Qubit Indices and overlap
    assert(ctrl.size() == shift.size() &&
           "ctrl and shift vectors must have the same size");

    // Create a set of target indices for quick lookup
    std::set<idx> target_set(target.begin(), target.end());

    // Validate control qubits and shift values
    for (idx c_idx = 0; c_idx < ctrl.size(); ++c_idx) {
        const idx c = ctrl[c_idx];
        assert(c < n && "Control qubit index c must be less than n");
        assert(target_set.find(c) == target_set.end() &&
               "Control qubit cannot also be a target qubit");
        assert(
            (shift[c_idx] == 0 || shift[c_idx] == 1) &&
            "Shift value must be 0 (positive control) or 1 (negative control)");
    }
#endif

    // Create the output state vector as a copy of the input.
    // This allows in-place updates via parallel access to non-overlapping
    // blocks and ensures blocks not meeting the control condition remain
    // unchanged.
    qpp::expr_t<Derived1> result = state;

    // --- Control Mask and Expected Pattern Calculation ---
    // Mask for controls that MUST be |1> (shift == 0, Positive Control)
    idx expected_pattern_for_ones = 0;
    // Mask for controls that MUST be |0> (shift == 1, Negative Control)
    idx expected_zero_mask = 0;

    // We iterate over the list of control qubits (ctrl)
    for (idx c_idx = 0; c_idx < ctrl.size(); ++c_idx) {
        const idx c = ctrl[c_idx];
        const idx j_c =
            n - 1 - c; // Big-endian bit position for control qubit c
        const idx bit = (1ULL << j_c); // Pre-calculate the bit to be set

        // If shift[c_idx] == 0: Positive Control (requires |1>)
        if (shift[c_idx] == 0) {
            expected_pattern_for_ones |= bit;
        }
        // If shift[c_idx] == 1: Negative Control (requires |0>)
        else {
            expected_zero_mask |= bit;
        }
    }

    // Pre-calculation of Inner Indices (Gate Index Mapping)
    // Maps local index r to its global index contribution.
    std::vector<idx> inner_idx(dim);
    for (idx r = 0; r < dim; ++r) {
        idx index_r = 0;
        for (idx j = 0; j < k; ++j) {
            const idx target_qubit_pos = target[k - 1 - j]; // Permutation fix
            if ((r >> j) & 1) {
                // BIG-ENDIAN contribution
                index_r += (1ULL << (n - 1 - target_qubit_pos));
            }
        }
        inner_idx[r] = index_r;
    }

    // Pre-calculation of Outer Indices (Spectator Block Base Index)
    // Maps spectator state index m to its global base index (i_base).
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
            if ((m >> j) & 1) {
                i_base += (1ULL << (n - 1 - spectator_qubits[j]));
            }
        }
        outer_idx[m] = i_base;
    }

    // Main Parallel Loop (Iterate over spectator blocks m)
#ifdef QPP_OPENMP
// NOLINTNEXTLINE
#pragma omp parallel for
#endif // QPP_OPENMP
    for (idx m = 0; m < outer_dim; ++m) {
        const idx i_base = outer_idx[m];

        // CONTROL CHECK: Apply gate A only if ALL control conditions are met
        // by the spectator state (represented by i_base)
        // 1. All positive control bits are 1: (i_base &
        // expected_pattern_for_ones) == expected_pattern_for_ones
        // 2. All negative control bits are 0: (i_base & expected_zero_mask) ==
        // 0
        if (((i_base & expected_pattern_for_ones) ==
             expected_pattern_for_ones) &&
            ((i_base & expected_zero_mask) == 0)) {

            // Input Block Caching (The Optimization)
            // Read the scattered input block (dim elements) into a contiguous
            // temporary vector.
            EigenVector input_block(dim);
            for (idx c = 0; c < dim; ++c) {
                input_block(c) = state(i_base + inner_idx[c]);
            }

            // Core Transformation (Block GEMV)
            // The output of the matrix multiplication is stored in new_block
            EigenVector new_block = A * input_block;

            // Write Output Block (Scattered Write)
            // Write the contiguous result block (new_block) back to the
            // scattered global indices.
            for (idx r = 0; r < dim; ++r) {
                const idx i_out = i_base + inner_idx[r];
                result(i_out) = new_block(r);
            }
        }
        // If control is not met, the original values in result (from the deep
        // copy) are kept.
    }

    return result;
}

/**
 * \brief Applies the controlled 1-qubit gate \a A to the qubit \a i of the
 * multi-partite density matrix \a state
 *
 * The gate is applied only if the state of the control qubits \a ctrl matches
 * the pattern defined by \a shift (0 for |1> control, 1 for |0> control)
 *
 * \param state Eigen expression
 * \param A Eigen expression (2x2 matrix)
 * \param i Target subsystem index
 * \param ctrl Vector of control qubit indices
 * \param shift Vector of control values (0: positive/|1> control, 1:
 * negative/|0> control)
 * \param n Number of qubits
 * \return Controlled gate \a A applied to the qubits \a i of \a state
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] expr_t<Derived1>
apply_ctrl_rho_1q(const Eigen::MatrixBase<Derived1>& state,
                  const Eigen::MatrixBase<Derived2>& A,
                  const std::vector<idx>& ctrl, idx i,
                  const std::vector<idx>& shift, idx n) {
    // Type and Dimension Setup
    using Scalar = typename Derived1::Scalar;
    using Matrix2 = Eigen::Matrix2<Scalar>;
    const idx D = 1ULL << n; // total Hilbert space dimension

    // Input Validation
#ifndef NDEBUG
    assert(i < n && "Target qubit index i must be less than n");
    assert(static_cast<idx>(state.rows()) == D &&
           static_cast<idx>(state.cols()) == D &&
           "State must be a square matrix sized 2^n x 2^n");
    assert(A.rows() == 2 && A.cols() == 2 && "Gate A must be a 2x2 matrix");
    assert(ctrl.size() == shift.size() &&
           "ctrl and shift vectors must have the same size");

    // Validate control qubits and overlap
    std::set<idx> target_set = {i};
    for (idx c_idx = 0; c_idx < ctrl.size(); ++c_idx) {
        const idx c = ctrl[c_idx];
        assert(c < n && "Control qubit index c must be less than n");
        assert(target_set.find(c) == target_set.end() &&
               "Control qubit cannot also be a target qubit");
        assert(
            (shift[c_idx] == 0 || shift[c_idx] == 1) &&
            "Shift value must be 0 (positive control) or 1 (negative control)");
    }
#endif

    // A deep copy is essential, as blocks not meeting the control condition
    // must remain unchanged (Case 4: I * rho * I).
    expr_t<Derived1> result = state;

    // --- Control Mask and Expected Pattern Calculation ---
    // The mask check will be applied to the row/column base index (which
    // represents the spectator state).
    idx expected_pattern_for_ones = 0; // Requires |1>
    idx expected_zero_mask = 0;        // Requires |0>

    for (idx c_idx = 0; c_idx < ctrl.size(); ++c_idx) {
        const idx c = ctrl[c_idx];
        const idx j_c =
            n - 1 - c; // Big-endian bit position for control qubit c
        const idx bit = (1ULL << j_c);

        if (shift[c_idx] == 0) {
            expected_pattern_for_ones |= bit;
        } else {
            expected_zero_mask |= bit;
        }
    }

    // Helper to check if a base index satisfies the control conditions
    auto control_is_met = [&](idx k_base) {
        return ((k_base & expected_pattern_for_ones) ==
                expected_pattern_for_ones) &&
               ((k_base & expected_zero_mask) == 0);
    };

    // --- Indexing Constants (same as apply_rho_1q) ---
    const idx p_i = n - 1 - i;
    const idx s_i = 1ULL << p_i;
    const idx D_spec = D / 2;
    const idx p_i_plus_1 = p_i + 1;
    const idx low_mask = s_i - 1;

    const Matrix2 A_block = A;
    const Matrix2 A_dagger = A_block.adjoint();

    const idx total_iterations = D_spec * D_spec;

    // Iterate over D_spec * D_spec blocks
#ifdef QPP_OPENMP
// NOLINTNEXTLINE
#pragma omp parallel for
#endif // QPP_OPENMP
    for (idx i = 0; i < total_iterations; ++i) {
        const idx s = i / D_spec;
        const idx s_prime = i % D_spec;

        // Calculate Row Indices (r0, r1) from spectator state 's'
        const idx s_high_r = s >> p_i;
        const idx s_low_r = s & low_mask;
        const idx r0 = (s_high_r << p_i_plus_1) | s_low_r;
        const idx r1 = r0 + s_i;

        const bool row_ctrl = control_is_met(r0);

        // Calculate Column Indices (c0, c1) from spectator state 's_prime'
        const idx s_prime_high_c = s_prime >> p_i;
        const idx s_prime_low_c = s_prime & low_mask;
        const idx c0 = (s_prime_high_c << p_i_plus_1) | s_prime_low_c;
        const idx c1 = c0 + s_i;

        const bool col_ctrl = control_is_met(c0);

        // Case 4: Skip if control not met for either row or column
        if (!row_ctrl && !col_ctrl) {
            continue;
        }

        // Extract the current 2x2 block rho[r0:r1, c0:c1]
        Matrix2 rho_block;
        rho_block(0, 0) = state.coeff(r0, c0);
        rho_block(0, 1) = state.coeff(r0, c1);
        rho_block(1, 0) = state.coeff(r1, c0);
        rho_block(1, 1) = state.coeff(r1, c1);

        Matrix2 result_block;
        if (row_ctrl && col_ctrl) {
            // Case 1: A * rho_block * A^\dagger
            result_block = A_block * rho_block * A_dagger;
        } else if (row_ctrl) {
            // Case 2: A * rho_block * I
            result_block = A_block * rho_block;
        } else { // col_ctrl must be true here
            // Case 3: I * rho_block * A^\dagger
            result_block = rho_block * A_dagger;
        }

        // Write the transformed 2x2 block back
        result.coeffRef(r0, c0) = result_block(0, 0);
        result.coeffRef(r0, c1) = result_block(0, 1);
        result.coeffRef(r1, c0) = result_block(1, 0);
        result.coeffRef(r1, c1) = result_block(1, 1);
    }

    return result;
}

/**
 * \brief Applies the controlled 2-qubit gate \a A to the qubits \a i and \a j
 * of the multi-partite density matrix \a state
 *
 * The gate is applied only if the state of the control qubits \a ctrl matches
 * the pattern defined by \a shift (0 for |1> control, 1 for |0> control)
 *
 * \param state Eigen expression
 * \param A Eigen expression (4x4 matrix)
 * \param i Target subsystem index
 * \param j Target subsystem index
 * \param ctrl Vector of control qubit indices
 * \param shift Vector of control values (0: positive/|1> control, 1:
 * negative/|0> control)
 * \param n Number of qubits
 * \return Controlled gate \a A applied to the qubits \a i and \a j of \a state
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] expr_t<Derived1>
apply_ctrl_rho_2q(const Eigen::MatrixBase<Derived1>& state,
                  const Eigen::MatrixBase<Derived2>& A,
                  const std::vector<idx>& ctrl, idx i, idx j,
                  const std::vector<idx>& shift, idx n) {
    // Input Validation
#ifndef NDEBUG
    const idx D = static_cast<idx>(state.rows());
    const idx D_expected = (1ULL << n);
    assert(i < n && j < n && i != j &&
           "Target qubit indices i and j must be distinct and less than n");
    assert(D == D_expected && D == static_cast<idx>(state.cols()) &&
           "State must be a square matrix sized 2^n x 2^n");
    assert(A.rows() == 4 && A.cols() == 4 && "Gate A must be a 4x4 matrix");
    assert(ctrl.size() == shift.size() &&
           "ctrl and shift vectors must have the same size");

    // Validate control qubits and overlap
    std::set<idx> target_set = {i, j};
    for (idx c_idx = 0; c_idx < ctrl.size(); ++c_idx) {
        const idx c = ctrl[c_idx];
        assert(c < n && "Control qubit index c must be less than n");
        assert(target_set.find(c) == target_set.end() &&
               "Control qubit cannot also be a target qubit");
        assert(
            (shift[c_idx] == 0 || shift[c_idx] == 1) &&
            "Shift value must be 0 (positive control) or 1 (negative control)");
    }
#endif

    // Index Conversion (MSB-first index 'i' -> LSB-first physical index
    // 'i_phys')
    const idx i_phys = n - i - 1;
    const idx j_phys = n - j - 1;

    // Type and Dimension Setup
    using Scalar = typename Derived1::Scalar;
    using ComputeMatrixType = Eigen::Matrix<Scalar, -1, -1>;
    using ComputeBlockType = Eigen::Matrix<Scalar, 4, 4>;

    // Use a deep copy for the result (rho_prime_cd) since blocks not meeting
    // the control condition are not modified.
    ComputeMatrixType rho_prime_cd = state.template cast<Scalar>().eval();
    // Use the original (copied) matrix for reading input blocks (rho_cd)
    const ComputeMatrixType& rho_cd = rho_prime_cd;

    const ComputeBlockType U = A.template cast<Scalar>().eval();
    const ComputeBlockType U_adj = U.adjoint();

    const idx D_rest = (n >= 2) ? (1ULL << (n - 2)) : 1;

    // Powers of 2 are calculated using the physical (LSB-first) indices.
    const idx P_i = 1ULL << i_phys;
    const idx P_j = 1ULL << j_phys;

    // --- Control Mask and Expected Pattern Calculation ---
    // The mask check must be done against the global index (Big Endian).
    idx expected_pattern_for_ones = 0; // Requires |1>
    idx expected_zero_mask = 0;        // Requires |0>

    for (idx c_idx = 0; c_idx < ctrl.size(); ++c_idx) {
        const idx c = ctrl[c_idx];
        const idx j_c =
            n - 1 - c; // Big-endian bit position for control qubit c
        const idx bit = (1ULL << j_c);

        if (shift[c_idx] == 0) {
            expected_pattern_for_ones |= bit;
        } else {
            expected_zero_mask |= bit;
        }
    }

    // Helper to check if a base index satisfies the control conditions
    auto control_is_met = [&](idx k_base) {
        return ((k_base & expected_pattern_for_ones) ==
                expected_pattern_for_ones) &&
               ((k_base & expected_zero_mask) == 0);
    };

    // Block Iteration (Parallelized)
#ifdef QPP_OPENMP
// NOLINTNEXTLINE
#pragma omp parallel for default(none)                                         \
    shared(D_rest, n, P_i, P_j, rho_cd, U, U_adj, rho_prime_cd, i_phys,        \
               j_phys, control_is_met)
#endif // QPP_OPENMP
    // Loop over row-blocks
    for (idx r = 0; r < D_rest; ++r) {

        // Variables private to this thread
        idx r_base_row = 0;
        idx current_r = r;

        // Calculate the 'base' index for the row block (rest qubits)
        for (idx q = 0; q < n;
             ++q) { // q is the physical index (0=LSB, n-1=MSB)
            if (q != i_phys && q != j_phys) {
                if (current_r & 1) {
                    r_base_row |= (1ULL << q);
                }
                current_r >>= 1;
            }
        }

        const bool row_ctrl = control_is_met(r_base_row);

        // Define the 4 *row* indices, aligned with the gate A's basis |q_i q_j>
        const idx vec_k_row[4] = {
            r_base_row,            // |00>
            r_base_row + P_j,      // |01> (LSB set)
            r_base_row + P_i,      // |10> (MSB set)
            r_base_row + P_i + P_j // |11> (Both set)
        };

        // Inner loop (col-blocks)
        for (idx c = 0; c < D_rest; ++c) {
            // Variables private to this inner iteration
            idx r_base_col = 0;
            idx current_c = c;
            ComputeBlockType M;
            ComputeBlockType M_prime;
            idx vec_k_col[4];

            // Calculate the 'base' index for the column block (rest qubits)
            for (idx q = 0; q < n; ++q) {
                if (q != i_phys && q != j_phys) {
                    if (current_c & 1) {
                        r_base_col |= (1ULL << q);
                    }
                    current_c >>= 1;
                }
            }

            const bool col_ctrl = control_is_met(r_base_col);

            // Case 4: Skip if control not met for both
            if (!row_ctrl && !col_ctrl) {
                continue;
            }

            // Define the 4 *column* indices
            vec_k_col[0] = r_base_col;
            vec_k_col[1] = r_base_col + P_j;
            vec_k_col[2] = r_base_col + P_i;
            vec_k_col[3] = r_base_col + P_i + P_j;

            // Extract the 4x4 block M from the *input* state
            for (int row = 0; row < 4; ++row) {
                for (int col = 0; col < 4; ++col) {
                    M(row, col) = rho_cd(vec_k_row[row], vec_k_col[col]);
                }
            }

            // Apply the gate based on control condition
            if (row_ctrl && col_ctrl) {
                // Case 1: U * M * U_adj
                M_prime = (U * M * U_adj).eval();
            } else if (row_ctrl) {
                // Case 2: U * M
                M_prime = (U * M).eval();
            } else { // col_ctrl must be true here
                // Case 3: M * U_adj
                M_prime = (M * U_adj).eval();
            }

            // Insert the resulting 4x4 block M_prime back into the result
            for (int row = 0; row < 4; ++row) {
                for (int col = 0; col < 4; ++col) {
                    rho_prime_cd(vec_k_row[row], vec_k_col[col]) =
                        M_prime(row, col);
                }
            }
        } // end col loop
    } // end row loop (parallelized)

    // Cast back to the original's scalar type
    return rho_prime_cd.template cast<typename Derived1::Scalar>();
}

/**
 * \brief Applies the controlled multi-qubit gate \a A to the part \a target of
 * the multi-partite density matrix \a state
 *
 * The gate is applied only if the state of the control qubits \a ctrl matches
 * the pattern defined by \a shift (0 for |1> control, 1 for |0> control)
 *
 * \param state Eigen expression
 * \param A Eigen expression
 * \param target Subsystem indexes where the gate \a A is applied
 * \param ctrl Vector of control qubit indices
 * \param shift Vector of control values (0: positive/|1> control, 1:
 * negative/|0> control)
 * \param n Number of qubits
 * \return Controlled gate \a A applied to the part \a target of \a state
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] expr_t<Derived1>
apply_ctrl_rho_kq(const Eigen::MatrixBase<Derived1>& state,
                  const Eigen::MatrixBase<Derived2>& A,
                  const std::vector<idx>& ctrl, const std::vector<idx>& target,
                  const std::vector<idx>& shift, idx n) {
    using Scalar = typename Derived1::Scalar;
    using ComputeMatrixType = Eigen::Matrix<Scalar, -1, -1>;

    const idx k = target.size(); // Gate qubits

    // Calculate block dimension D_k = 2^k
    const idx D_k = (k == 0) ? 1 : (1ULL << k);

    // --- Input Validation ---
#ifndef NDEBUG
    const idx D = static_cast<idx>(state.rows());
    // Standard dimension checks
    assert(static_cast<idx>(A.rows()) == D_k &&
           static_cast<idx>(A.cols()) == D_k &&
           "Gate A must be a 2^k x 2^k matrix");
    assert(static_cast<idx>(state.rows()) == D &&
           static_cast<idx>(state.cols()) == D &&
           "State must be a 2^n x 2^n matrix");
    assert(ctrl.size() == shift.size() &&
           "ctrl and shift vectors must have the same size");

    // Validate target and control qubit indices and check for overlap
    // Use a temporary vector to check distinctness and overlap efficiently
    std::vector<idx> check_set;
    check_set.reserve(target.size() + ctrl.size());
    check_set.insert(check_set.end(), target.begin(), target.end());
    check_set.insert(check_set.end(), ctrl.begin(), ctrl.end());

    // Check target indices distinctness
    std::vector<idx> sorted_target = target;
    std::sort(sorted_target.begin(), sorted_target.end());
    assert(std::unique(sorted_target.begin(), sorted_target.end()) ==
               sorted_target.end() &&
           "Target qubit indices must be distinct");

    // Check overlap between target and controls
    std::vector<idx> sorted_check_set = check_set;
    std::sort(sorted_check_set.begin(), sorted_check_set.end());
    assert(std::unique(sorted_check_set.begin(), sorted_check_set.end()) ==
               sorted_check_set.end() &&
           "Target and control qubit indices must not overlap");

    for (idx i = 0; i < check_set.size(); ++i) {
        assert(check_set[i] < n && "Qubit index must be less than n");
    }

    for (idx s : shift) {
        assert(
            (s == 0 || s == 1) &&
            "Shift value must be 0 (positive control) or 1 (negative control)");
    }
#endif

    // --- Pre-Calculations for Index Mapping ---

    // 1. Convert user's logical (MSB-first) target indices to physical
    // (LSB-first) indices and determine their powers of 2 (P_gate_basis).
    std::vector<idx> target_phys(k);
    std::vector<idx> P_gate_basis(k);

    // Marker to quickly determine which physical indices are target qubits
    std::vector<bool> is_target_phys(n, false);

    for (idx l = 0; l < k; ++l) {
        const idx phys_idx = n - target[l] - 1; // Physical index (0=LSB)
        target_phys[l] = phys_idx;
        P_gate_basis[l] = 1ULL << phys_idx;
        is_target_phys[phys_idx] = true;
    }

    // 2. Determine the physical indices that are *not* in the target set (the
    // "rest").
    std::vector<idx> rest_phys;
    rest_phys.reserve(n);

    // This efficiently calculates the complement of target_phys within [0, n-1]
    for (idx q_phys = 0; q_phys < n; ++q_phys) {
        if (!is_target_phys[q_phys]) {
            rest_phys.push_back(q_phys);
        }
    }
    // Size of the N-k qubit subsystem
    const idx D_rest = (n >= k) ? (1ULL << (n - k)) : 1;

    // --- Control Mask and Expected Pattern Calculation ---
    // Masks are based on global index positions (MSB-first logic).
    idx expected_pattern_for_ones =
        0;                      // Bits that must be 1 for control to meet
    idx expected_zero_mask = 0; // Bits that must be 0 for control to meet

    for (idx c_idx = 0; c_idx < ctrl.size(); ++c_idx) {
        const idx c = ctrl[c_idx];
        const idx bit = (1ULL << (n - 1 - c)); // Global index position marker

        if (shift[c_idx] == 0) { // Positive control (|1>)
            expected_pattern_for_ones |= bit;
        } else { // Negative control (|0>)
            expected_zero_mask |= bit;
        }
    }

    // Helper to check if a global base index satisfies the control conditions
    auto control_is_met = [&](idx k_base) {
        // All required '1' bits are set AND all required '0' bits are unset
        return ((k_base & expected_pattern_for_ones) ==
                expected_pattern_for_ones) &&
               ((k_base & expected_zero_mask) == 0);
    };

    // --- Matrix Setup ---
    const ComputeMatrixType U = A.template cast<Scalar>().eval();
    const ComputeMatrixType U_adj = U.adjoint();

    // Initialize result as a deep copy. Blocks that are not transformed (Case
    // 4) will remain unchanged (I * M * I = M).
    ComputeMatrixType rho_prime_cd = state.template cast<Scalar>().eval();
    // Read from the copy, write back to it
    const ComputeMatrixType& rho_cd = rho_prime_cd;

    // Block Iteration (Parallelized)
#ifdef QPP_OPENMP
// NOLINTNEXTLINE
#pragma omp parallel for default(none)                                         \
    shared(D_rest, n, k, rest_phys, P_gate_basis, rho_cd, U, U_adj,            \
               rho_prime_cd, D_k, control_is_met)
#endif // QPP_OPENMP
    // Loop over row-blocks (rest index)
    for (idx r = 0; r < D_rest; ++r) {

        // Calculate the 'base' index for the row block (rest qubits)
        idx r_base_row = 0;
        idx current_r = r;
        for (const auto& q_phys : rest_phys) {
            if (current_r & 1) {
                r_base_row |= (1ULL << q_phys);
            }
            current_r >>= 1;
        }

        const bool row_ctrl = control_is_met(r_base_row);

        // Pre-calculate the D_k *row* indices for the block.
        std::vector<idx> vec_k_row(D_k);
        for (idx m = 0; m < D_k; ++m) {
            idx target_component = 0;
            for (idx l = 0; l < k; ++l) {
                // l=0 is target[0], which corresponds to MSB in the gate
                idx bit_pos = k - 1 - l;
                if ((m >> bit_pos) & 1) {
                    target_component += P_gate_basis[l];
                }
            }
            vec_k_row[m] = r_base_row + target_component;
        }

        // Inner loop (col-blocks)
        for (idx c = 0; c < D_rest; ++c) {
            // Calculate the 'base' index for the column block (rest
            // qubits)
            idx r_base_col = 0;
            idx current_c = c;
            for (const auto& q_phys : rest_phys) {
                if (current_c & 1) {
                    r_base_col |= (1ULL << q_phys);
                }
                current_c >>= 1;
            }

            const bool col_ctrl = control_is_met(r_base_col);

            // Case 4: Skip if control not met for both (block is I*M*I =
            // M, already present)
            if (!row_ctrl && !col_ctrl) {
                continue;
            }

            ComputeMatrixType M(D_k, D_k);
            ComputeMatrixType M_prime;
            std::vector<idx> vec_k_col(D_k);

            // Pre-calculate the D_k *column* indices
            for (idx m = 0; m < D_k; ++m) {
                idx target_component = 0;
                for (idx l = 0; l < k; ++l) {
                    idx bit_pos = k - 1 - l;
                    if ((m >> bit_pos) & 1) {
                        target_component += P_gate_basis[l];
                    }
                }
                vec_k_col[m] = r_base_col + target_component;
            }

            // Extract the D_k x D_k block M from the *input* state
            // (rho_cd)
            for (idx row = 0; row < D_k; ++row) {
                for (idx col = 0; col < D_k; ++col) {
                    M(row, col) = rho_cd(vec_k_row[row], vec_k_col[col]);
                }
            }

            // Apply the transformation: M' = U_ctrl * M * U_ctrl^\dagger
            if (row_ctrl && col_ctrl) {
                // Case 1: U * M * U_adj
                M_prime = (U * M * U_adj).eval();
            } else if (row_ctrl) {
                // Case 2: U * M (U * M * I)
                M_prime = (U * M).eval();
            } else { // col_ctrl must be true here
                // Case 3: M * U_adj (I * M * U_adj)
                M_prime = (M * U_adj).eval();
            }

            // Insert the resulting block M_prime back into the result
            // (rho_prime_cd)
            for (idx row = 0; row < D_k; ++row) {
                for (idx col = 0; col < D_k; ++col) {
                    rho_prime_cd(vec_k_row[row], vec_k_col[col]) =
                        M_prime(row, col);
                }
            }
        } // end col loop
    } // end row loop (parallelized)

    // Return
    return rho_prime_cd.template cast<typename Derived1::Scalar>();
}

} /* namespace internal */
} /* namespace qpp */

#endif /* QPP_INTERNAL_CRITICAL_HPP_ */
