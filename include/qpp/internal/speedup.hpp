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
 * \file qpp/internal/speedup.hpp
 * \brief Internal highly optimized functions
 */

#ifndef QPP_INTERNAL_SPEEDUP_HPP_
#define QPP_INTERNAL_SPEEDUP_HPP_

#include <vector>

#include <Eigen/Dense>

#include "qpp/types.hpp"

namespace qpp {
namespace internal {
/**
 * \brief Applies the 1-qubit gate \a A to the qubit \a i of the
 * multi-partite state vector \a state
 *
 * \param state Eigen expression
 * \param A Eigen expression
 * \param i Subsystem index where the gate \a A is applied
 * \paran n Number of qubits
 * \return Gate \a A applied to the qubit \a target of \a state
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] expr_t<Derived1>
apply_psi_1q(const Eigen::MatrixBase<Derived1>& state,
             const Eigen::MatrixBase<Derived2>& A, idx i, idx n) {
    // 1. Type and Dimension Setup
    using Scalar = typename Derived1::Scalar;
    const idx D = 1ULL << n; // Total size of the state vector (2^n)

    // A deep copy of the input state is required. The new state must be
    // calculated using the *old* values of psi_k0 and psi_k1 simultaneously.
    expr_t<Derived1> result = state;

    // 2. Pre-calculate Bit-Shift Constants (BIG-ENDIAN)
    // The effective bit position 'j' corresponding to qubit 'i' in big-endian
    // is (n - 1 - i).
    const idx j = n - 1 - i;

    // 'step' is 2^j. This is the difference in index between |...0...> and
    // |...1...> at bit j.
    const idx step = 1ULL << j;

    // 'jump' is 2^(j+1). This is the size of the block that repeats.
    const idx jump = 1ULL << (j + 1);

    // 3. Extract the 2x2 gate elements
    const Scalar a00 = A.coeff(0, 0);
    const Scalar a01 = A.coeff(0, 1);
    const Scalar a10 = A.coeff(1, 0);
    const Scalar a11 = A.coeff(1, 1);

    // 4. Pair-wise Amplitude Transformation
    // The outer loop (L) iterates over all blocks of size 'jump'.
    // This loop is perfectly independent and is the primary target for
    // parallelization.
#ifdef QPP_OPENMP
// NOLINTNEXTLINE
#pragma omp parallel for collapse(1)
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
            result.coeffRef(k0) = a00 * psi_k0 + a01 * psi_k1;

            // psi'_k1 = A[1,0] * psi_k0 + A[1,1] * psi_k1
            result.coeffRef(k1) = a10 * psi_k0 + a11 * psi_k1;
        }
    }

    return result;
}

/**
 * \brief Applies the 2-qubit gate \a A to the qubits \a i and \a j of the
 * multi-partite state vector \a state
 *
 * \param state Eigen expression
 * \param A Eigen expression
 * \param i Subsystem index where the gate \a A is applied
 * \param j Subsystem index where the gate \a A is applied
 * \paran n Number of qubits
 * \return Gate \a A applied to the qubit \a target of \a state
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] expr_t<Derived1>
apply_psi_2q(const Eigen::MatrixBase<Derived1>& state,
             const Eigen::MatrixBase<Derived2>& A, idx i, idx j, idx n) {

    // 1. Type and Dimension Setup
    using Scalar = typename Derived1::Scalar;
    const idx D = 1ULL << n; // Total size of the state vector (2^n)

    // A deep copy of the input state is required.
    expr_t<Derived1> result = state;

    // 2. Pre-calculate Bit-Shift Constants (BIG-ENDIAN)
    // The effective bit position 'p' for a physical qubit 'q' in big-endian is
    // (n - 1 - q).
    const idx p_i = n - 1 - i;
    const idx p_j = n - 1 - j;

    // 's_i' and 's_j' are the power-of-2 values that flip the bit at p_i and
    // p_j.
    const idx s_i = 1ULL << p_i;
    const idx s_j = 1ULL << p_j;

    // Determine the overall 'jump' and 'step' for the loop structure.
    // The loop jump is determined by the highest (most significant) bit
    // position involved, plus one.
    const idx p_max = std::max(p_i, p_j);
    const idx jump = 1ULL << (p_max + 1);

    // The inner loop iterates up to the smallest 'step' value.
    const idx step_min = std::min(s_i, s_j);

    // 3. Extract the 4x4 gate elements
    // Using .coeff() for fast, direct access. A is assumed to be 4x4.
    const Scalar u00 = A.coeff(0, 0);
    const Scalar u01 = A.coeff(0, 1);
    const Scalar u02 = A.coeff(0, 2);
    const Scalar u03 = A.coeff(0, 3);

    const Scalar u10 = A.coeff(1, 0);
    const Scalar u11 = A.coeff(1, 1);
    const Scalar u12 = A.coeff(1, 2);
    const Scalar u13 = A.coeff(1, 3);

    const Scalar u20 = A.coeff(2, 0);
    const Scalar u21 = A.coeff(2, 1);
    const Scalar u22 = A.coeff(2, 2);
    const Scalar u23 = A.coeff(2, 3);

    const Scalar u30 = A.coeff(3, 0);
    const Scalar u31 = A.coeff(3, 1);
    const Scalar u32 = A.coeff(3, 2);
    const Scalar u33 = A.coeff(3, 3);

    // 4. Pair-wise Amplitude Transformation
    // The outer loop (L) iterates over the higher bits not involved in the
    // gate.
#ifdef QPP_OPENMP
// NOLINTNEXTLINE
#pragma omp parallel for collapse(1)
#endif // QPP_OPENMP
    for (idx L = 0; L < D; L += jump) {
        // The inner loop (R) iterates over the lower bits not involved in the
        // gate. It runs up to the step of the lower bit position, effectively
        // iterating over the bits below the target qubits.
        for (idx R = 0; R < step_min; ++R) {

            // Calculate the four coupled indices based on L, R, s_i, and s_j
            // k00: |...0_i 0_j...>
            const idx k00 = L + R;
            // k01: |...0_i 1_j...>
            const idx k01 = k00 + s_j;
            // k10: |...1_i 0_j...>
            const idx k10 = k00 + s_i;
            // k11: |...1_i 1_j...>
            const idx k11 = k00 + s_i + s_j;

            // Fetch the ORIGINAL amplitudes from the INPUT 'state'
            const Scalar psi00 = state.coeff(k00);
            const Scalar psi01 = state.coeff(k01);
            const Scalar psi10 = state.coeff(k10);
            const Scalar psi11 = state.coeff(k11);

            // Apply the 4x4 matrix-vector multiplication U * [psi00, psi01,
            // psi10, psi11]^T

            // result_k00 (Row 0 of U)
            result.coeffRef(k00) =
                (u00 * psi00) + (u01 * psi01) + (u02 * psi10) + (u03 * psi11);

            // result_k01 (Row 1 of U)
            result.coeffRef(k01) =
                (u10 * psi00) + (u11 * psi01) + (u12 * psi10) + (u13 * psi11);

            // result_k10 (Row 2 of U)
            result.coeffRef(k10) =
                (u20 * psi00) + (u21 * psi01) + (u22 * psi10) + (u23 * psi11);

            // result_k11 (Row 3 of U)
            result.coeffRef(k11) =
                (u30 * psi00) + (u31 * psi01) + (u32 * psi10) + (u33 * psi11);
        }
    }

    return result;
}

/**
 * \brief Applies the 3-qubit gate \a A to the qubit \a i, \a j, and \a k of the
 * multi-partite state vector \a state
 *
 * \param state Eigen expression
 * \param A Eigen expression
 * \param i Subsystem index where the gate \a A is applied
 * \param j Subsystem index where the gate \a A is applied
 * \param k Subsystem index where the gate \a A is applied
 * \paran n Number of qubits
 * \return Gate \a A applied to the qubit \a target of \a state
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] expr_t<Derived1>
apply_psi_3q(const Eigen::MatrixBase<Derived1>& state,
             const Eigen::MatrixBase<Derived2>& A, idx i, idx j, idx k, idx n) {
}

/**
 * \brief Applies the multi-qubit gate \a A to the part \a target of the
 * multi-partite state vector \a state
 *
 * \param state Eigen expression
 * \param A Eigen expression
 * \param target Subsystem indexes where the gate \a A is applied
 * \return Gate \a A applied to the part \a target of \a state
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] expr_t<Derived1>
apply_psi_multi(const Eigen::MatrixBase<Derived1>& state,
                const Eigen::MatrixBase<Derived2>& A,
                const std::vector<idx>& target) {}

/**
 * \brief Applies the 1-qubit gate \a A to the qubit \a i of the
 * multi-partite density matrix \a state
 *
 * \param state Eigen expression
 * \param A Eigen expression
 * \param i Subsystem index where the gate \a A is applied
 * \paran n Number of qubits
 * \return Gate \a A applied to the qubit \a target of \a state
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] expr_t<Derived1>
apply_rho_1q(const Eigen::MatrixBase<Derived1>& state,
             const Eigen::MatrixBase<Derived2>& A, idx i, idx n) {}

/**
 * \brief Applies the 2-qubit gate \a A to the qubits \a i and \a j of the
 * multi-partite density matrix \a state
 *
 * \param state Eigen expression
 * \param A Eigen expression
 * \param i Subsystem index where the gate \a A is applied
 * \param j Subsystem index where the gate \a A is applied
 * \paran n Number of qubits
 * \return Gate \a A applied to the qubit \a target of \a state
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] expr_t<Derived1>
apply_rho_2q(const Eigen::MatrixBase<Derived1>& state,
             const Eigen::MatrixBase<Derived2>& A, idx i, idx j, idx n) {}

/**
 * \brief Applies the 3-qubit gate \a A to the qubit \a i, \a j, and \a k of the
 * multi-partite density matrix \a state
 *
 * \param state Eigen expression
 * \param A Eigen expression
 * \param i Subsystem index where the gate \a A is applied
 * \param j Subsystem index where the gate \a A is applied
 * \param k Subsystem index where the gate \a A is applied
 * \paran n Number of qubits
 * \return Gate \a A applied to the qubit \a target of \a state
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] expr_t<Derived1>
apply_rho_3q(const Eigen::MatrixBase<Derived1>& state,
             const Eigen::MatrixBase<Derived2>& A, idx i, idx j, idx k, idx n) {
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
apply_rho_multi(const Eigen::MatrixBase<Derived1>& state,
                const Eigen::MatrixBase<Derived2>& A,
                const std::vector<idx>& target) {}

} /* namespace internal */
} /* namespace qpp */

#endif /* QPP_INTERNAL_SPEEDUP_HPP_ */
