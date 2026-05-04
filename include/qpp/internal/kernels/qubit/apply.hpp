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
 * @file qpp/internal/kernels/qubit/apply.hpp
 * @brief Internal highly optimized critical functions for qpp::apply()
 */

#ifndef QPP_INTERNAL_KERNELS_QUBIT_APPLY_HPP_
#define QPP_INTERNAL_KERNELS_QUBIT_APPLY_HPP_

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
 * @brief Applies the 1-qubit gate \a A to the qubit \a i of the
 * multi-partite state vector \a state in-place
 *
 * @param state Eigen expression (modified in-place)
 * @param A Eigen expression (2x2 matrix)
 * @param i Subsystem index where the gate \a A is applied
 * @param n Number of qubits
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] void
apply_psi_1q_inplace(Eigen::MatrixBase<Derived1>& state,
                     const Eigen::MatrixBase<Derived2>& A, idx i, idx n) {
    // Type and Dimension Setup
    using Scalar = typename Derived1::Scalar;
    const idx D = static_cast<idx>(
        std::size_t{1} << n); // Total size of the state vector (2^n)

    // Input Validation
    assert(i < n && "Target qubit index i must be less than n");
    assert(static_cast<idx>(state.size()) == D &&
           "State vector size must be 2^n");
    assert(A.rows() == 2 && A.cols() == 2 && "Gate A must be a 2x2 matrix");

    // Pre-calculate Bit-Shift Constants (BIG-ENDIAN)
    // The effective bit position 'j' corresponding to qubit 'i' in big-endian
    // is (n - 1 - i).
    const idx j = n - 1 - i;

    // 'step' is 2^j. This is the difference in index between |...0...> and
    // |...1...> at bit j.
    const idx step = static_cast<idx>(std::size_t{1} << j);

    // 'jump' is 2^(j+1). This is the size of the block that repeats.
    const idx jump = static_cast<idx>(std::size_t{1} << (j + 1));

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
            // We must store these in local variables to allow in-place update
            const Scalar psi_k0 = state.coeff(k0);
            const Scalar psi_k1 = state.coeff(k1);

            // Apply the 2x2 matrix-vector multiplication and write to 'state'
            // psi'_k0 = A[0,0] * psi_k0 + A[0,1] * psi_k1
            state.coeffRef(k0) = (a00 * psi_k0) + (a01 * psi_k1);

            // psi'_k1 = A[1,0] * psi_k0 + A[1,1] * psi_k1
            state.coeffRef(k1) = (a10 * psi_k0) + (a11 * psi_k1);
        }
    }
}

/**
 * @brief Applies the 1-qubit gate \a A to the qubit \a i of the
 * multi-partite state vector \a state
 *
 * @param state Eigen expression
 * @param A Eigen expression (2x2 matrix)
 * @param i Subsystem index where the gate \a A is applied
 * @param n Number of qubits
 * @return Gate \a A applied to the qubit \a i of \a state
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] expr_t<Derived1>
apply_psi_1q(const Eigen::MatrixBase<Derived1>& state,
             const Eigen::MatrixBase<Derived2>& A, idx i, idx n) {
    // A deep copy of the input state is required
    expr_t<Derived1> result = state;

    // Apply the gate in-place on the result copy
    apply_psi_1q_inplace(result, A, i, n);

    return result;
}

/**
 * @brief Applies the 2-qubit gate \a A to the qubits \a i and \a j of the
 * multi-partite state vector \a state in-place
 *
 * @param state Eigen expression (modified in-place)
 * @param A Eigen expression (4x4 matrix)
 * @param i Subsystem index where the gate \a A is applied
 * @param j Subsystem index where the gate \a A is applied
 * @param n Number of qubits
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] void
apply_psi_2q_inplace(Eigen::MatrixBase<Derived1>& state,
                     const Eigen::MatrixBase<Derived2>& A, idx i, idx j,
                     idx n) {
    // Type and Dimension Setup
    using Scalar = typename Derived1::Scalar;

    // Input Validation
    assert(i < n && j < n && i != j &&
           "Target qubit indices i and j must be distinct and less than n");
    assert(A.rows() == 4 && A.cols() == 4 && "Gate A must be a 4x4 matrix");
#ifndef NDEBUG
    const idx D = static_cast<idx>(
        std::size_t{1} << n); // Total size of the state vector (2^n)
    assert(static_cast<idx>(state.size()) == D &&
           "State vector size must be 2^n");
#endif

    // Pre-calculate Bit-Shift Constants (BIG-ENDIAN)
    // The effective bit position 'p' for a physical qubit 'q' in big-endian is
    // (n - 1 - q).
    const idx p_i = n - 1 - i;
    const idx p_j = n - 1 - j;

    // 's_i' and 's_j' are the power-of-2 values that flip the bit at p_i and
    // p_j.
    const idx s_i = static_cast<idx>(std::size_t{1} << p_i);
    const idx s_j = static_cast<idx>(std::size_t{1} << p_j);

    // The number of remaining bits that define the independent blocks (2^(n-2)
    // total blocks)
    const idx D_rem = static_cast<idx>(std::size_t{1} << (n - 2));

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
        // The current bit position in 'block_idx' (0 to n-3)
        idx current_block_bit = 0;

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
            if ((static_cast<std::size_t>(block_idx) >> current_block_bit) &
                std::size_t{1}) {
                k00 |= static_cast<idx>(std::size_t{1} << p);
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
        // Stored in local variables to allow in-place modification of 'state'
        const Scalar psi00 = state.coeff(k00);
        const Scalar psi01 = state.coeff(k01);
        const Scalar psi10 = state.coeff(k10);
        const Scalar psi11 = state.coeff(k11);

        // Apply the 4x4 matrix-vector multiplication A * [psi00, psi01,
        // psi10, psi11]^T

        // Row 0 -> state_k00 (|0_i 0_j>)
        state.coeffRef(k00) =
            (a00 * psi00) + (a01 * psi01) + (a02 * psi10) + (a03 * psi11);

        // Row 1 -> state_k01 (|0_i 1_j>)
        state.coeffRef(k01) =
            (a10 * psi00) + (a11 * psi01) + (a12 * psi10) + (a13 * psi11);

        // Row 2 -> state_k10 (|1_i 0_j>)
        state.coeffRef(k10) =
            (a20 * psi00) + (a21 * psi01) + (a22 * psi10) + (a23 * psi11);

        // Row 3 -> state_k11 (|1_i 1_j>)
        state.coeffRef(k11) =
            (a30 * psi00) + (a31 * psi01) + (a32 * psi10) + (a33 * psi11);
    }
}

/**
 * @brief Applies the 2-qubit gate \a A to the qubits \a i and \a j of the
 * multi-partite state vector \a state
 *
 * @param state Eigen expression
 * @param A Eigen expression (4x4 matrix)
 * @param i Subsystem index where the gate \a A is applied
 * @param j Subsystem index where the gate \a A is applied
 * @param n Number of qubits
 * @return Gate \a A applied to the qubit \a i and \a j of \a state
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] expr_t<Derived1>
apply_psi_2q(const Eigen::MatrixBase<Derived1>& state,
             const Eigen::MatrixBase<Derived2>& A, idx i, idx j, idx n) {
    // A deep copy of the input state is required.
    expr_t<Derived1> result = state;

    // Delegate to in-place implementation
    apply_psi_2q_inplace(result, A, i, j, n);

    return result;
}

/**
 * @brief Applies the 3-qubit gate \a A to the qubit \a i, \a j, and \a k of the
 * multi-partite state vector \a state in-place
 *
 * @param state Eigen expression (modified in-place)
 * @param A Eigen expression (8x8 matrix)
 * @param i Subsystem index where the gate \a A is applied
 * @param j Subsystem index where the gate \a A is applied
 * @param k Subsystem index where the gate \a A is applied
 * @param n Number of qubits
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] void
apply_psi_3q_inplace(Eigen::MatrixBase<Derived1>& state,
                     const Eigen::MatrixBase<Derived2>& A, idx i, idx j, idx k,
                     idx n) {
    // Type and Dimension Setup
    using Scalar = typename Derived1::Scalar;

    // Input Validation
    assert(i < n && j < n && k < n && i != j && i != k && j != k &&
           "Target qubit indices i, j, k must be distinct and less than n");
    assert(A.rows() == 8 && A.cols() == 8 && "Gate A must be an 8x8 matrix");
#ifndef NDEBUG
    const idx D = static_cast<idx>(
        std::size_t{1} << n); // Total size of the state vector (2^n)
    assert(static_cast<idx>(state.size()) == D &&
           "State vector size must be 2^n");
#endif

    // Pre-calculate Bit-Shift Constants (BIG-ENDIAN)
    // The effective bit position 'p' for a physical qubit 'q' in big-endian is
    // (n - 1 - q).
    const idx p_i = n - 1 - i;
    const idx p_j = n - 1 - j;
    const idx p_k = n - 1 - k;

    // 's_i', 's_j', 's_k' are the power-of-2 values that flip the bit at p_i,
    // p_j, p_k.
    const idx s_i = static_cast<idx>(std::size_t{1} << p_i);
    const idx s_j = static_cast<idx>(std::size_t{1} << p_j);
    const idx s_k = static_cast<idx>(std::size_t{1} << p_k);

    // The number of remaining bits that define the independent blocks (2^(n-3)
    // total blocks)
    const idx D_rem = static_cast<idx>(std::size_t{1} << (n - 3));

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
            if ((static_cast<std::size_t>(block_idx) >> current_block_bit) &
                std::size_t{1}) {
                k000 |= static_cast<idx>(std::size_t{1} << p);
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
        // We must store these in local variables to allow in-place update
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

        // Row 0 -> state_k000 (|000>)
        state.coeffRef(k000) =
            (A00 * psi000) + (A01 * psi001) + (A02 * psi010) + (A03 * psi011) +
            (A04 * psi100) + (A05 * psi101) + (A06 * psi110) + (A07 * psi111);

        // Row 1 -> state_k001 (|001>)
        state.coeffRef(k001) =
            (A10 * psi000) + (A11 * psi001) + (A12 * psi010) + (A13 * psi011) +
            (A14 * psi100) + (A15 * psi101) + (A16 * psi110) + (A17 * psi111);

        // Row 2 -> state_k010 (|010>)
        state.coeffRef(k010) =
            (A20 * psi000) + (A21 * psi001) + (A22 * psi010) + (A23 * psi011) +
            (A24 * psi100) + (A25 * psi101) + (A26 * psi110) + (A27 * psi111);

        // Row 3 -> state_k011 (|011>)
        state.coeffRef(k011) =
            (A30 * psi000) + (A31 * psi001) + (A32 * psi010) + (A33 * psi011) +
            (A34 * psi100) + (A35 * psi101) + (A36 * psi110) + (A37 * psi111);

        // Row 4 -> state_k100 (|100>)
        state.coeffRef(k100) =
            (A40 * psi000) + (A41 * psi001) + (A42 * psi010) + (A43 * psi011) +
            (A44 * psi100) + (A45 * psi101) + (A46 * psi110) + (A47 * psi111);

        // Row 5 -> state_k101 (|101>)
        state.coeffRef(k101) =
            (A50 * psi000) + (A51 * psi001) + (A52 * psi010) + (A53 * psi011) +
            (A54 * psi100) + (A55 * psi101) + (A56 * psi110) + (A57 * psi111);

        // Row 6 -> state_k110 (|110>)
        state.coeffRef(k110) =
            (A60 * psi000) + (A61 * psi001) + (A62 * psi010) + (A63 * psi011) +
            (A64 * psi100) + (A65 * psi101) + (A66 * psi110) + (A67 * psi111);

        // Row 7 -> state_k111 (|111>)
        state.coeffRef(k111) =
            (A70 * psi000) + (A71 * psi001) + (A72 * psi010) + (A73 * psi011) +
            (A74 * psi100) + (A75 * psi101) + (A76 * psi110) + (A77 * psi111);
    }
}

/**
 * @brief Applies the 3-qubit gate \a A to the qubit \a i, \a j, and \a k of the
 * multi-partite state vector \a state
 *
 * @param state Eigen expression
 * @param A Eigen expression (8x8 matrix)
 * @param i Subsystem index where the gate \a A is applied
 * @param j Subsystem index where the gate \a A is applied
 * @param k Subsystem index where the gate \a A is applied
 * @param n Number of qubits
 * @return Gate \a A applied to the qubit \a i, \a j, and \a k of \a state
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] expr_t<Derived1>
apply_psi_3q(const Eigen::MatrixBase<Derived1>& state,
             const Eigen::MatrixBase<Derived2>& A, idx i, idx j, idx k, idx n) {
    // A deep copy of the input state is required.
    expr_t<Derived1> result = state;

    // Delegate to in-place implementation
    apply_psi_3q_inplace(result, A, i, j, k, n);

    return result;
}

/**
 * @brief Applies the multi-qubit gate \a A to the part \a target of the
 * multi-partite state vector \a state in-place
 *
 * @param state Eigen expression (modified in-place)
 * @param A Eigen expression
 * @param target Subsystem indexes where the gate \a A is applied
 * @param n Number of qubits
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] void
apply_psi_kq_inplace(Eigen::MatrixBase<Derived1>& state,
                     const Eigen::MatrixBase<Derived2>& A,
                     const std::vector<idx>& target, idx n) {
    // Type Aliases
    using Scalar = typename Derived1::Scalar;
    using EigenVector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

    // Setup
    const idx k = static_cast<idx>(target.size());
    const idx dim = static_cast<idx>(std::size_t{1} << k);
    const idx outer_dim = static_cast<idx>(std::size_t{1} << (n - k));

    // Input Validation
    // Check Gate Dimension: A must be a 2^k x 2^k matrix
    assert(static_cast<idx>(A.rows()) == dim &&
           static_cast<idx>(A.cols()) == dim &&
           "Gate A must be a 2^k x 2^k matrix, where k is the number of target "
           "qubits");
#ifndef NDEBUG
    // D is the dimension of the state (2^n)
    const idx D = static_cast<idx>(std::size_t{1} << n);

    // Check Target Qubit Indices are valid
    for (idx t : target) {
        assert(t < n && "Target qubit index must be less than n");
    }

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

    // Pre-calculation of Inner Indices (Gate Index Mapping)
    // Maps local index r to its global index contribution.
    std::vector<idx> inner_idx(dim);
    for (idx r = 0; r < dim; ++r) {
        idx index_r = 0;
        for (idx j = 0; j < k; ++j) {
            const idx target_qubit_pos = target[k - 1 - j]; // Permutation fix
            if ((static_cast<std::size_t>(r) >> j) & std::size_t{1}) {
                // BIG-ENDIAN contribution
                index_r += static_cast<idx>(std::size_t{1}
                                            << (n - 1 - target_qubit_pos));
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
            if ((static_cast<std::size_t>(m) >> j) & std::size_t{1}) {
                i_base += static_cast<idx>(std::size_t{1}
                                           << (n - 1 - spectator_qubits[j]));
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
        // perfect cache locality and allows the state to be modified in-place.
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
            state(i_out) = new_block(r);
        }
    }
}

/**
 * @brief Applies the multi-qubit gate \a A to the part \a target of the
 * multi-partite state vector \a state
 *
 * @param state Eigen expression
 * @param A Eigen expression
 * @param target Subsystem indexes where the gate \a A is applied
 * @param n Number of qubits
 * @return Gate \a A applied to the part \a target of \a state
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] qpp::expr_t<Derived1>
apply_psi_kq(const Eigen::MatrixBase<Derived1>& state,
             const Eigen::MatrixBase<Derived2>& A,
             const std::vector<idx>& target, idx n) {
    // Create the output state vector as a copy of the input.
    qpp::expr_t<Derived1> result = state;

    // Delegate to in-place implementation
    apply_psi_kq_inplace(result, A, target, n);

    return result;
}

/**
 * @brief Applies the 1-qubit gate \a A to the qubit \a i of the
 * multi-partite density matrix \a state in-place
 *
 * @param state Eigen expression (modified in-place)
 * @param A Eigen expression (2x2 matrix)
 * @param i Subsystem index where the gate \a A is applied
 * @param n Number of qubits
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] void
apply_rho_1q_inplace(Eigen::MatrixBase<Derived1>& state,
                     const Eigen::MatrixBase<Derived2>& A, idx i, idx n) {
    using Scalar = typename Derived1::Scalar;
    using Matrix2 = Eigen::Matrix<Scalar, 2, 2>;
    const idx D = static_cast<idx>(std::size_t{1} << n);

    // Input Validation
    assert(static_cast<idx>(state.rows()) == D &&
           static_cast<idx>(state.cols()) == D &&
           "State must be a square matrix sized 2^n x 2^n");
    assert(i < n && "Target qubit index i must be less than n");
    assert(A.rows() == 2 && A.cols() == 2 && "Gate A must be a 2x2 matrix");

    // Bit position for qubit i (p_i) and the stride (s_i = 2^p_i)
    const idx p_i = n - 1 - i;
    const idx s_i = static_cast<idx>(std::size_t{1} << p_i);

    // Constants for index reconstruction
    const idx D_spec = D / 2; // Total number of spectator states (2^(n-1))
    const idx p_i_plus_1 = p_i + 1;
    const idx low_mask = s_i - 1; // Mask for bits below p_i

    const Matrix2 A_block = A;
    const Matrix2 A_dagger = A_block.adjoint();

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

        // Extract the current 2x2 block from 'state'
        Matrix2 rho_block;
        rho_block(0, 0) = state.coeff(r0, c0);
        rho_block(0, 1) = state.coeff(r0, c1);
        rho_block(1, 0) = state.coeff(r1, c0);
        rho_block(1, 1) = state.coeff(r1, c1);

        // Apply transformation: A * rho_block * A^\dagger
        const Matrix2 result_block = A_block * rho_block * A_dagger;

        // Write transformed block back to 'state' in-place
        state.coeffRef(r0, c0) = result_block(0, 0);
        state.coeffRef(r0, c1) = result_block(0, 1);
        state.coeffRef(r1, c0) = result_block(1, 0);
        state.coeffRef(r1, c1) = result_block(1, 1);
    }
}

/**
 * @brief Applies the 1-qubit gate \a A to the qubit \a i of the
 * multi-partite density matrix \a state
 *
 * @param state Eigen expression
 * @param A Eigen expression (2x2 matrix)
 * @param i Subsystem index where the gate \a A is applied
 * @param n Number of qubits
 * @return Gate \a A applied to the qubit \a i of \a state
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] expr_t<Derived1>
apply_rho_1q(const Eigen::MatrixBase<Derived1>& state,
             const Eigen::MatrixBase<Derived2>& A, idx i, idx n) {
    // A deep copy of the input state is required.
    expr_t<Derived1> result = state;

    // Delegate to in-place implementation
    apply_rho_1q_inplace(result, A, i, n);

    return result;
}

/**
 * @brief Applies the 2-qubit gate \a A to the qubits \a i and \a j of the
 * multi-partite density matrix \a state in-place
 *
 * @param state Eigen expression (modified in-place)
 * @param A Eigen expression (4x4 matrix)
 * @param i Subsystem index where the gate \a A is applied
 * @param j Subsystem index where the gate \a A is applied
 * @param n Number of qubits
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] void
apply_rho_2q_inplace(Eigen::MatrixBase<Derived1>& state,
                     const Eigen::MatrixBase<Derived2>& A, idx i, idx j,
                     idx n) {
    const idx D [[maybe_unused]] = static_cast<idx>(state.rows());

    // Input Validation
    assert(i < n && j < n && i != j &&
           "Target qubit indices i and j must be distinct and less than n");
    assert(A.rows() == 4 && A.cols() == 4 && "Gate A must be a 4x4 matrix");
#ifndef NDEBUG
    const idx D_expected = static_cast<idx>(std::size_t{1} << n);
    assert(D == D_expected && D == static_cast<idx>(state.cols()) &&
           "State must be a square matrix sized 2^n x 2^n");
#endif

    // Index Conversion (Big Endian to Little Endian)
    const idx i_phys = n - i - 1;
    const idx j_phys = n - j - 1;

    // Type and Dimension Setup
    using Scalar = typename Derived1::Scalar;
    using ComputeBlockType = Eigen::Matrix<Scalar, 4, 4>;

    const ComputeBlockType U = A.template cast<Scalar>();
    const ComputeBlockType U_adj = U.adjoint();

    // Size of the N-2 qubit subsystem
    const idx D_rest =
        (n >= 2) ? static_cast<idx>(std::size_t{1} << (n - 2)) : 1;

    const idx P_i = static_cast<idx>(std::size_t{1} << i_phys);
    const idx P_j = static_cast<idx>(std::size_t{1} << j_phys);

    // Block Iteration (Parallelized)
#ifdef QPP_OPENMP
// NOLINTNEXTLINE
#pragma omp parallel for default(none)                                         \
    shared(D_rest, n, P_i, P_j, state, U, U_adj, i_phys, j_phys)
#endif // QPP_OPENMP
    for (idx r = 0; r < D_rest; ++r) {
        idx r_base_row = 0;
        idx current_r = r;

        for (idx q = 0; q < n; ++q) {
            if (q != i_phys && q != j_phys) {
                if (static_cast<std::size_t>(current_r) & std::size_t{1}) {
                    r_base_row |= static_cast<idx>(std::size_t{1} << q);
                }
                current_r >>= 1;
            }
        }

        const idx vec_k_row[4] = {r_base_row, r_base_row + P_j,
                                  r_base_row + P_i, r_base_row + P_i + P_j};

        for (idx c = 0; c < D_rest; ++c) {
            idx r_base_col = 0;
            idx current_c = c;
            ComputeBlockType M;
            ComputeBlockType M_prime;
            idx vec_k_col[4];

            for (idx q = 0; q < n; ++q) {
                if (q != i_phys && q != j_phys) {
                    if (static_cast<std::size_t>(current_c) & std::size_t{1}) {
                        r_base_col |= static_cast<idx>(std::size_t{1} << q);
                    }
                    current_c >>= 1;
                }
            }

            vec_k_col[0] = r_base_col;
            vec_k_col[1] = r_base_col + P_j;
            vec_k_col[2] = r_base_col + P_i;
            vec_k_col[3] = r_base_col + P_i + P_j;

            // Extract the 4x4 block M from the state
            for (int row = 0; row < 4; ++row) {
                for (int col = 0; col < 4; ++col) {
                    M(row, col) = state.coeff(vec_k_row[row], vec_k_col[col]);
                }
            }

            // Apply the gate: M' = U * M * U_adj
            M_prime.noalias() = U * M * U_adj;

            // Insert back in-place
            for (int row = 0; row < 4; ++row) {
                for (int col = 0; col < 4; ++col) {
                    state.coeffRef(vec_k_row[row], vec_k_col[col]) =
                        M_prime(row, col);
                }
            }
        }
    }
}

/**
 * @brief Applies the 2-qubit gate \a A to the qubits \a i and \a j of the
 * multi-partite density matrix \a state
 *
 * @param state Eigen expression
 * @param A Eigen expression (4x4 matrix)
 * @param i Subsystem index where the gate \a A is applied
 * @param j Subsystem index where the gate \a A is applied
 * @param n Number of qubits
 * @return Gate \a A applied to the qubits \a i and \a j of \a state
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] expr_t<Derived1>
apply_rho_2q(const Eigen::MatrixBase<Derived1>& state,
             const Eigen::MatrixBase<Derived2>& A, idx i, idx j, idx n) {
    // Deep copy required for functional interface
    expr_t<Derived1> result = state;

    // Delegate to the in-place version
    apply_rho_2q_inplace(result, A, i, j, n);

    return result;
}

/**
 * @brief Applies the 3-qubit gate \a A to the qubit \a i, \a j, and \a k of the
 * multi-partite density matrix \a state in-place
 *
 * @param state Eigen expression (modified in-place)
 * @param A Eigen expression (8x8 matrix)
 * @param i Subsystem index where the gate \a A is applied
 * @param j Subsystem index where the gate \a A is applied
 * @param k Subsystem index where the gate \a A is applied
 * @param n Number of qubits
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] void
apply_rho_3q_inplace(Eigen::MatrixBase<Derived1>& state,
                     const Eigen::MatrixBase<Derived2>& A, idx i, idx j, idx k,
                     idx n) {
    const idx D [[maybe_unused]] = static_cast<idx>(state.rows());

    // Input Validation
    assert(i < n && j < n && k < n && i != j && i != k && j != k &&
           "Target qubit indices i, j, and k must be distinct and less than n");
    assert(A.rows() == 8 && A.cols() == 8 && "Gate A must be an 8x8 matrix");
    assert(n >= 3 && "Need at least 3 qubits for a 3-qubit gate");
#ifndef NDEBUG
    const idx D_expected = static_cast<idx>(std::size_t{1} << n);
    assert(D == D_expected && D == static_cast<idx>(state.cols()) &&
           "State must be a square matrix of size 2^n x 2^n");
#endif

    // Index Conversion (Big Endian to Little Endian)
    const idx i_phys = n - i - 1;
    const idx j_phys = n - j - 1;
    const idx k_phys = n - k - 1;

    // Type and Dimension Setup
    using Scalar = typename Derived1::Scalar;
    using ComputeBlockType = Eigen::Matrix<Scalar, 8, 8>;

    const ComputeBlockType U = A.template cast<Scalar>();
    const ComputeBlockType U_adj = U.adjoint();

    // Size of the N-3 qubit subsystem
    const idx D_rest =
        (n >= 3) ? static_cast<idx>(std::size_t{1} << (n - 3)) : 1;

    const idx P_i = static_cast<idx>(std::size_t{1} << i_phys);
    const idx P_j = static_cast<idx>(std::size_t{1} << j_phys);
    const idx P_k = static_cast<idx>(std::size_t{1} << k_phys);

    // Block Iteration (Parallelized)
#ifdef QPP_OPENMP
// NOLINTNEXTLINE
#pragma omp parallel for default(none)                                         \
    shared(D_rest, n, P_i, P_j, P_k, state, U, U_adj, i_phys, j_phys, k_phys)
#endif // QPP_OPENMP
    for (idx r = 0; r < D_rest; ++r) {
        idx r_base_row = 0;
        idx current_r = r;

        for (idx q = 0; q < n; ++q) {
            if (q != i_phys && q != j_phys && q != k_phys) {
                if (static_cast<std::size_t>(current_r) & std::size_t{1}) {
                    r_base_row |= static_cast<idx>(std::size_t{1} << q);
                }
                current_r >>= 1;
            }
        }

        const idx vec_k_row[8] = {r_base_row,
                                  r_base_row + P_k,
                                  r_base_row + P_j,
                                  r_base_row + P_j + P_k,
                                  r_base_row + P_i,
                                  r_base_row + P_i + P_k,
                                  r_base_row + P_i + P_j,
                                  r_base_row + P_i + P_j + P_k};

        for (idx c = 0; c < D_rest; ++c) {
            idx r_base_col = 0;
            idx current_c = c;
            ComputeBlockType M;
            ComputeBlockType M_prime;
            idx vec_k_col[8];

            for (idx q = 0; q < n; ++q) {
                if (q != i_phys && q != j_phys && q != k_phys) {
                    if (static_cast<std::size_t>(current_c) & std::size_t{1}) {
                        r_base_col |= static_cast<idx>(std::size_t{1} << q);
                    }
                    current_c >>= 1;
                }
            }

            vec_k_col[0] = r_base_col;
            vec_k_col[1] = r_base_col + P_k;
            vec_k_col[2] = r_base_col + P_j;
            vec_k_col[3] = r_base_col + P_j + P_k;
            vec_k_col[4] = r_base_col + P_i;
            vec_k_col[5] = r_base_col + P_i + P_k;
            vec_k_col[6] = r_base_col + P_i + P_j;
            vec_k_col[7] = r_base_col + P_i + P_j + P_k;

            // Extract the 8x8 block M from the state
            for (int row = 0; row < 8; ++row) {
                for (int col = 0; col < 8; ++col) {
                    M(row, col) = state.coeff(vec_k_row[row], vec_k_col[col]);
                }
            }

            // Apply the gate: M' = U * M * U_adj
            M_prime.noalias() = U * M * U_adj;

            // Insert back in-place
            for (int row = 0; row < 8; ++row) {
                for (int col = 0; col < 8; ++col) {
                    state.coeffRef(vec_k_row[row], vec_k_col[col]) =
                        M_prime(row, col);
                }
            }
        }
    }
}

/**
 * @brief Applies the 3-qubit gate \a A to the qubit \a i, \a j, and \a k of the
 * multi-partite density matrix \a state
 *
 * @param state Eigen expression
 * @param A Eigen expression (8x8 matrix)
 * @param i Subsystem index where the gate \a A is applied
 * @param j Subsystem index where the gate \a A is applied
 * @param k Subsystem index where the gate \a A is applied
 * @param n Number of qubits
 * @return Gate \a A applied to the qubit \a i, \a j, and \a k of \a state
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] expr_t<Derived1>
apply_rho_3q(const Eigen::MatrixBase<Derived1>& state,
             const Eigen::MatrixBase<Derived2>& A, idx i, idx j, idx k, idx n) {
    // Create copy for functional behavior
    expr_t<Derived1> result = state;

    // Perform the operation in-place on the copy
    apply_rho_3q_inplace(result, A, i, j, k, n);

    return result;
}

/**
 * @brief Applies the multi-qubit gate \a A to the part \a target of the
 * multi-partite density matrix \a state in-place
 *
 * @param state Eigen expression (modified in-place)
 * @param A Eigen expression
 * @param target Subsystem indexes where the gate \a A is applied
 * @param n Number of qubits
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] void
apply_rho_kq_inplace(Eigen::MatrixBase<Derived1>& state,
                     const Eigen::MatrixBase<Derived2>& A,
                     const std::vector<idx>& target, idx n) {
    using Scalar = typename Derived1::Scalar;
    using ComputeMatrixType = Eigen::Matrix<Scalar, -1, -1>;

    const idx D [[maybe_unused]] = static_cast<idx>(state.rows());
    const idx k = target.size();

    const idx D_k = (k == 0) ? 1 : static_cast<idx>(std::size_t{1} << k);

    // Input Validation
    assert(static_cast<idx>(A.rows()) == D_k &&
           static_cast<idx>(A.cols()) == D_k &&
           "Gate A must be a 2^k x 2^k matrix");
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

    // Pre-Calculations
    std::vector<idx> target_phys(k);
    std::vector<idx> P_gate_basis(k);
    for (idx l = 0; l < k; ++l) {
        target_phys[l] = n - target[l] - 1;
        P_gate_basis[l] = static_cast<idx>(std::size_t{1} << target_phys[l]);
    }

    std::vector<idx> rest_phys = fast_complement(target_phys, n);
    const idx D_rest =
        (n >= k) ? static_cast<idx>(std::size_t{1} << (n - k)) : 1;

    const ComputeMatrixType U = A.template cast<Scalar>();
    const ComputeMatrixType U_adj = U.adjoint();

    // Block Iteration (Parallelized)
#ifdef QPP_OPENMP
// NOLINTNEXTLINE
#pragma omp parallel for default(none)                                         \
    shared(D_rest, n, k, rest_phys, target, target_phys, P_gate_basis, state,  \
               U, U_adj, D_k)
#endif // QPP_OPENMP
    for (idx r = 0; r < D_rest; ++r) {
        idx r_base_row = 0;
        idx current_r = r;

        for (const auto& q_phys : rest_phys) {
            if (static_cast<std::size_t>(current_r) & std::size_t{1}) {
                r_base_row |= static_cast<idx>(std::size_t{1} << q_phys);
            }
            current_r >>= 1;
        }

        std::vector<idx> vec_k_row(D_k);
        for (idx m = 0; m < D_k; ++m) {
            idx target_component = 0;
            for (idx l = 0; l < k; ++l) {
                idx bit_pos = k - 1 - l;
                if ((static_cast<std::size_t>(m) >> bit_pos) & std::size_t{1}) {
                    target_component += P_gate_basis[l];
                }
            }
            vec_k_row[m] = r_base_row + target_component;
        }

        for (idx c = 0; c < D_rest; ++c) {
            idx r_base_col = 0;
            idx current_c = c;
            ComputeMatrixType M(D_k, D_k);
            ComputeMatrixType M_prime;
            std::vector<idx> vec_k_col(D_k);

            for (const auto& q_phys : rest_phys) {
                if (static_cast<std::size_t>(current_c) & std::size_t{1}) {
                    r_base_col |= static_cast<idx>(std::size_t{1} << q_phys);
                }
                current_c >>= 1;
            }

            for (idx m = 0; m < D_k; ++m) {
                idx target_component = 0;
                for (idx l = 0; l < k; ++l) {
                    idx bit_pos = k - 1 - l;
                    if ((static_cast<std::size_t>(m) >> bit_pos) &
                        std::size_t{1}) {
                        target_component += P_gate_basis[l];
                    }
                }
                vec_k_col[m] = r_base_col + target_component;
            }

            // Extract the D_k x D_k block M
            for (idx row = 0; row < D_k; ++row) {
                for (idx col = 0; col < D_k; ++col) {
                    M(row, col) = state.coeff(vec_k_row[row], vec_k_col[col]);
                }
            }

            // Apply the gate: M' = U * M * U_adj
            M_prime.noalias() = U * M * U_adj;

            // Insert back in-place
            for (idx row = 0; row < D_k; ++row) {
                for (idx col = 0; col < D_k; ++col) {
                    state.coeffRef(vec_k_row[row], vec_k_col[col]) =
                        M_prime(row, col);
                }
            }
        }
    }
}

/**
 * @brief Applies the multi-qubit gate \a A to the part \a target of the
 * multi-partite density matrix \a state
 *
 * @param state Eigen expression
 * @param A Eigen expression
 * @param target Subsystem indexes where the gate \a A is applied
 * @param n Number of qubits
 * @return Gate \a A applied to the part \a target of \a state
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] expr_t<Derived1>
apply_rho_kq(const Eigen::MatrixBase<Derived1>& state,
             const Eigen::MatrixBase<Derived2>& A,
             const std::vector<idx>& target, idx n) {
    // Functional deep copy
    expr_t<Derived1> result = state;

    // Apply transformation in-place on the copy
    apply_rho_kq_inplace(result, A, target, n);

    return result;
}
} // namespace qpp::internal::kernels::qubit

#endif /* QPP_INTERNAL_KERNELS_QUBIT_APPLY_HPP_ */
