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
 * @file qpp/internal/kernels/qubit/apply_ctrl.hpp
 * @brief Internal highly optimized critical functions for qpp::applyCTRL()
 */

#ifndef QPP_INTERNAL_KERNELS_QUBIT_APPLY_CTRL_HPP_
#define QPP_INTERNAL_KERNELS_QUBIT_APPLY_CTRL_HPP_

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
 * @brief Applies the controlled 1-qubit gate \a A to the qubit \a i of the
 * multi-partite state vector \a state
 *
 * The gate is applied only if the state of the control qubits \a ctrl matches
 * the pattern defined by \a shift (0 for |1> control, 1 for |0> control)
 *
 * @param state Eigen expression
 * @param A Eigen expression (2x2 matrix)
 * @param i Target subsystem index
 * @param ctrl Vector of control qubit indices
 * @param shift Vector of control values (0: positive/|1> control, 1:
 * negative/|0> control)
 * @param n Number of qubits
 * @return Controlled gate \a A applied to the qubit \a i of \a state
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
    assert(i < n && "Target qubit index i must be less than n");
    assert(static_cast<idx>(state.size()) == D &&
           "State vector size must be 2^n");
    assert(A.rows() == 2 && A.cols() == 2 && "Gate A must be a 2x2 matrix");
    assert(ctrl.size() == shift.size() &&
           "ctrl and shift vectors must have the same size");
    const idx ctrl_size = static_cast<idx>(ctrl.size());
#ifndef NDEBUG
    // Validate control qubits and shift values
    for (idx c_idx = 0; c_idx < ctrl_size; ++c_idx) {
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
    for (idx c_idx = 0; c_idx < ctrl_size; ++c_idx) {
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
 * @brief Applies the controlled 2-qubit gate \a A to the qubits \a i and \a j
 * of the multi-partite state vector \a state
 *
 * The gate is applied only if the state of the control qubits \a ctrl matches
 * the pattern defined by \a shift (0 for |1> control, 1 for |0> control)
 *
 * @param state Eigen expression
 * @param A Eigen expression (4x4 matrix)
 * @param ctrl Vector of control qubit indices
 * @param i Target subsystem index
 * @param j Target subsystem index
 * @param shift Vector of control values (0: positive/|1> control, 1:
 * negative/|0> control)
 * @param n Number of qubits
 * @return Controlled gate \a A applied to the qubits \a i and \a j of \a state
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] expr_t<Derived1>
apply_ctrl_psi_2q(const Eigen::MatrixBase<Derived1>& state,
                  const Eigen::MatrixBase<Derived2>& A,
                  const std::vector<idx>& ctrl, idx i, idx j,
                  const std::vector<idx>& shift, idx n) {
    using Scalar = typename Derived1::Scalar;
    const idx D = 1ULL << n;

    // Input Validation
    assert(i < n && j < n && i != j &&
           "Target qubit indices i and j must be distinct and less than n");
    assert(static_cast<idx>(state.size()) == D &&
           "State vector size must be 2^n");
    assert(A.rows() == 4 && A.cols() == 4 && "Gate A must be a 4x4 matrix");
    assert(ctrl.size() == shift.size() &&
           "ctrl and shift vectors must have the same size");
    const idx ctrl_size = static_cast<idx>(ctrl.size());
#ifndef NDEBUG
    for (idx c_idx = 0; c_idx < ctrl_size; ++c_idx) {
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
    for (idx c_idx = 0; c_idx < ctrl_size; ++c_idx) {
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
 * @brief Applies the multi-controlled qubit gate \a A to the part \a target of
 * the multi-partite state vector \a state
 *
 * The gate is applied only if the state of the control qubits \a ctrl matches
 * the pattern defined by \a shift (0 for |1> control, 1 for |0> control)
 *
 * @param state Eigen expression
 * @param A Eigen expression (4x4 matrix)
 * @param ctrl Vector of control qubit indices
 * @param target Subsystem indexes where the gate \a A is applied
 * @param shift Vector of control values (0: positive/|1> control, 1:
 * negative/|0> control)
 * @param n Number of qubits
 * @return Controlled gate \a A applied to the part \a target of \a state
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
    // Check Gate Dimension: A must be a 2^k x 2^k matrix
    assert(static_cast<idx>(A.rows()) == dim &&
           static_cast<idx>(A.cols()) == dim &&
           "Gate A must be a 2^k x 2^k matrix, where k is the number of target "
           "qubits");
    // Check Control Qubit Indices and overlap
    assert(ctrl.size() == shift.size() &&
           "ctrl and shift vectors must have the same size");
    const idx ctrl_size = static_cast<idx>(ctrl.size());
#ifndef NDEBUG
    // D is the dimension of the state (2^n)
    const idx D = 1ULL << n;

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

    // Create a set of target indices for quick lookup
    std::set<idx> target_set(target.begin(), target.end());

    // Validate control qubits and shift values
    for (idx c_idx = 0; c_idx < ctrl_size; ++c_idx) {
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
    for (idx c_idx = 0; c_idx < ctrl_size; ++c_idx) {
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
 * @brief Applies the controlled 1-qubit gate \a A to the qubit \a i of the
 * multi-partite density matrix \a state
 *
 * The gate is applied only if the state of the control qubits \a ctrl matches
 * the pattern defined by \a shift (0 for |1> control, 1 for |0> control)
 *
 * @param state Eigen expression
 * @param A Eigen expression (2x2 matrix)
 * @param i Target subsystem index
 * @param ctrl Vector of control qubit indices
 * @param shift Vector of control values (0: positive/|1> control, 1:
 * negative/|0> control)
 * @param n Number of qubits
 * @return Controlled gate \a A applied to the qubits \a i of \a state
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
    assert(i < n && "Target qubit index i must be less than n");
    assert(static_cast<idx>(state.rows()) == D &&
           static_cast<idx>(state.cols()) == D &&
           "State must be a square matrix sized 2^n x 2^n");
    assert(A.rows() == 2 && A.cols() == 2 && "Gate A must be a 2x2 matrix");
    assert(ctrl.size() == shift.size() &&
           "ctrl and shift vectors must have the same size");
    const idx ctrl_size = static_cast<idx>(ctrl.size());
#ifndef NDEBUG
    // Validate control qubits and overlap
    std::set<idx> target_set = {i};
    for (idx c_idx = 0; c_idx < ctrl_size; ++c_idx) {
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

    for (idx c_idx = 0; c_idx < ctrl_size; ++c_idx) {
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
 * @brief Applies the controlled 2-qubit gate \a A to the qubits \a i and \a j
 * of the multi-partite density matrix \a state
 *
 * The gate is applied only if the state of the control qubits \a ctrl matches
 * the pattern defined by \a shift (0 for |1> control, 1 for |0> control)
 *
 * @param state Eigen expression
 * @param A Eigen expression (4x4 matrix)
 * @param i Target subsystem index
 * @param j Target subsystem index
 * @param ctrl Vector of control qubit indices
 * @param shift Vector of control values (0: positive/|1> control, 1:
 * negative/|0> control)
 * @param n Number of qubits
 * @return Controlled gate \a A applied to the qubits \a i and \a j of \a state
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] expr_t<Derived1>
apply_ctrl_rho_2q(const Eigen::MatrixBase<Derived1>& state,
                  const Eigen::MatrixBase<Derived2>& A,
                  const std::vector<idx>& ctrl, idx i, idx j,
                  const std::vector<idx>& shift, idx n) {
    // Input Validation
    assert(i < n && j < n && i != j &&
           "Target qubit indices i and j must be distinct and less than n");
    assert(A.rows() == 4 && A.cols() == 4 && "Gate A must be a 4x4 matrix");
    assert(ctrl.size() == shift.size() &&
           "ctrl and shift vectors must have the same size");
    const idx ctrl_size = static_cast<idx>(ctrl.size());
#ifndef NDEBUG
    const idx D = static_cast<idx>(state.rows());
    const idx D_expected = (1ULL << n);
    assert(D == D_expected && D == static_cast<idx>(state.cols()) &&
           "State must be a square matrix sized 2^n x 2^n");

    // Validate control qubits and overlap
    std::set<idx> target_set = {i, j};
    for (idx c_idx = 0; c_idx < ctrl_size; ++c_idx) {
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

    for (idx c_idx = 0; c_idx < ctrl_size; ++c_idx) {
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
 * @brief Applies the controlled multi-qubit gate \a A to the part \a target of
 * the multi-partite density matrix \a state
 *
 * The gate is applied only if the state of the control qubits \a ctrl matches
 * the pattern defined by \a shift (0 for |1> control, 1 for |0> control)
 *
 * @param state Eigen expression
 * @param A Eigen expression
 * @param target Subsystem indexes where the gate \a A is applied
 * @param ctrl Vector of control qubit indices
 * @param shift Vector of control values (0: positive/|1> control, 1:
 * negative/|0> control)
 * @param n Number of qubits
 * @return Controlled gate \a A applied to the part \a target of \a state
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

    // Input Validation
    // Standard dimension checks
    assert(static_cast<idx>(A.rows()) == D_k &&
           static_cast<idx>(A.cols()) == D_k &&
           "Gate A must be a 2^k x 2^k matrix");
    assert(static_cast<idx>(A.rows()) == D_k &&
           static_cast<idx>(A.cols()) == D_k &&
           "Gate A must be a 2^k x 2^k matrix");
    assert(ctrl.size() == shift.size() &&
           "ctrl and shift vectors must have the same size");
    const idx ctrl_size = static_cast<idx>(ctrl.size());
#ifndef NDEBUG
    const idx D = static_cast<idx>(state.rows());
    assert(static_cast<idx>(state.rows()) == D &&
           static_cast<idx>(state.cols()) == D &&
           "State must be a 2^n x 2^n matrix");

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

    for (idx i = 0; i < static_cast<idx>(check_set.size()); ++i) {
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

    for (idx l = 0; l < k; ++l) {
        const idx phys_idx = n - target[l] - 1; // Physical index (0=LSB)
        target_phys[l] = phys_idx;
        P_gate_basis[l] = 1ULL << phys_idx;
    }

    // 2. Determine the physical indices that are *not* in the target set (the
    // "rest").
    std::vector<idx> rest_phys = fast_complement(target_phys, n);

    // Size of the N-k qubit subsystem
    const idx D_rest = (n >= k) ? (1ULL << (n - k)) : 1;

    // --- Control Mask and Expected Pattern Calculation ---
    // Masks are based on global index positions (MSB-first logic).
    idx expected_pattern_for_ones =
        0;                      // Bits that must be 1 for control to meet
    idx expected_zero_mask = 0; // Bits that must be 0 for control to meet

    for (idx c_idx = 0; c_idx < ctrl_size; ++c_idx) {
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

    return rho_prime_cd.template cast<typename Derived1::Scalar>();
}
} // namespace qpp::internal::kernels::qubit

#endif /* QPP_INTERNAL_KERNELS_QUBIT_APPLY_CTRL_HPP_ */
