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
 * @file qpp/internal/kernels/qubit/apply_ctrl_diag.hpp
 * @brief Internal highly optimized critical functions for qpp::applyCTRL_diag()
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
 * @brief Applies the controlled 1-qubit diagonal gate \a A to the qubit \a i of
 * the multi-partite state vector \a state in-place
 *
 * @param state Eigen expression (modified in-place)
 * @param A Eigen expression (2 x 1 vector of diagonal elements)
 * @param i Target subsystem index
 * @param ctrl Vector of control qubit indices
 * @param shift Vector of control values (0: positive/|1> control, 1:
 * negative/|0> control)
 * @param n Number of qubits
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] void apply_ctrl_psi_1q_diag_inplace(
    Eigen::MatrixBase<Derived1>& state, const Eigen::MatrixBase<Derived2>& A,
    const std::vector<idx>& ctrl, idx i, const std::vector<idx>& shift, idx n) {
    using Scalar = typename Derived1::Scalar;
    const idx D = static_cast<idx>(std::size_t{1} << n);

    // Input Validation
    assert(i < n && "Target qubit index i must be less than n");
    assert(static_cast<idx>(state.size()) == D &&
           "State vector size must be 2^n");
    assert(A.size() == 2 && "Diagonal gate A must have 2 elements");
    assert(ctrl.size() == shift.size() &&
           "ctrl and shift vectors must have the same size");
    const idx ctrl_size = static_cast<idx>(ctrl.size());

#ifndef NDEBUG
    for (idx c_idx = 0; c_idx < ctrl_size; ++c_idx) {
        const idx c = ctrl[c_idx];
        assert(c < n && "Control qubit index c must be less than n");
        assert(c != i && "Target qubit i cannot also be a control qubit");
        assert((shift[c_idx] == 0 || shift[c_idx] == 1) &&
               "Shift value must be 0 or 1");
    }
#endif

    const idx j = n - 1 - i;
    const idx step = static_cast<idx>(std::size_t{1} << j);
    const idx jump = static_cast<idx>(std::size_t{1} << (j + 1));

    idx expected_pattern_for_ones = 0;
    idx expected_zero_mask = 0;

    for (idx c_idx = 0; c_idx < ctrl_size; ++c_idx) {
        const idx c = ctrl[c_idx];
        const idx bit = static_cast<idx>(std::size_t{1} << (n - 1 - c));
        if (shift[c_idx] == 0) {
            expected_pattern_for_ones |= bit;
        } else {
            expected_zero_mask |= bit;
        }
    }

    // Extract diagonal elements
    const Scalar a0 = A.coeff(0);
    const Scalar a1 = A.coeff(1);

#ifdef QPP_OPENMP
// NOLINTNEXTLINE
#pragma omp parallel for
#endif // QPP_OPENMP
    for (idx L = 0; L < D; L += jump) {
        for (idx R = 0; R < step; ++R) {
            const idx k0 = L + R;
            const idx k1 = k0 + step;

            // Check if control conditions are met
            if (((static_cast<std::size_t>(k0) &
                  static_cast<std::size_t>(expected_pattern_for_ones)) ==
                 static_cast<std::size_t>(expected_pattern_for_ones)) &&
                ((static_cast<std::size_t>(k0) &
                  static_cast<std::size_t>(expected_zero_mask)) == 0)) {

                // Apply diagonal elements directly
                state.coeffRef(k0) *= a0;
                state.coeffRef(k1) *= a1;
            }
        }
    }
}

/**
 * @brief Applies the controlled 1-qubit diagonal gate \a A to the qubit \a i of
 * the multi-partite state vector \a state
 *
 * @param state Eigen expression
 * @param A Eigen expression (2 x 1 vector of diagonal elements)
 * @param i Target subsystem index
 * @param ctrl Vector of control qubit indices
 * @param shift Vector of control values
 * @param n Number of qubits
 * @return Controlled gate \a A applied to the qubit \a i of \a state
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] expr_t<Derived1>
apply_ctrl_psi_1q_diag(const Eigen::MatrixBase<Derived1>& state,
                       const Eigen::MatrixBase<Derived2>& A,
                       const std::vector<idx>& ctrl, idx i,
                       const std::vector<idx>& shift, idx n) {
    // Deep copy for functional return
    expr_t<Derived1> result = state;

    // Apply transformation in-place
    apply_ctrl_psi_1q_diag_inplace(result, A, ctrl, i, shift, n);

    return result;
}

/**
 * @brief Applies the controlled 2-qubit diagonal gate \a A to the qubits \a i
 * and \a j of the multi-partite state vector \a state in-place
 *
 * @param state Eigen expression (modified in-place)
 * @param A Eigen expression (4 x 1 vector of diagonal elements)
 * @param ctrl Vector of control qubit indices
 * @param i Target subsystem index
 * @param j Target subsystem index
 * @param shift Vector of control values (0: positive/|1> control, 1:
 * negative/|0> control)
 * @param n Number of qubits
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] void
apply_ctrl_psi_2q_diag_inplace(Eigen::MatrixBase<Derived1>& state,
                               const Eigen::MatrixBase<Derived2>& A,
                               const std::vector<idx>& ctrl, idx i, idx j,
                               const std::vector<idx>& shift, idx n) {
    using Scalar = typename Derived1::Scalar;
    const idx D = static_cast<idx>(std::size_t{1} << n);

    // Input Validation
    assert(i < n && j < n && i != j &&
           "Target qubit indices i and j must be distinct and less than n");
    assert(static_cast<idx>(state.size()) == D &&
           "State vector size must be 2^n");
    assert(A.size() == 4 && "Diagonal gate A must have 4 elements");
    assert(ctrl.size() == shift.size() &&
           "ctrl and shift vectors must have the same size");
    const idx ctrl_size = static_cast<idx>(ctrl.size());

#ifndef NDEBUG
    for (idx c_idx = 0; c_idx < ctrl_size; ++c_idx) {
        const idx c = ctrl[c_idx];
        assert(c < n && "Control qubit index c must be less than n");
        assert(c != i && c != j &&
               "Target qubits cannot also be control qubits");
        assert((shift[c_idx] == 0 || shift[c_idx] == 1) &&
               "Shift value must be 0 or 1");
    }
#endif

    // Big-endian bit positions
    const idx p_i = n - 1 - i;
    const idx p_j = n - 1 - j;
    const idx s_i = static_cast<idx>(std::size_t{1} << p_i);
    const idx s_j = static_cast<idx>(std::size_t{1} << p_j);

    idx control_mask = 0;
    idx control_value = 0;
    for (idx c_idx = 0; c_idx < ctrl_size; ++c_idx) {
        const idx c = ctrl[c_idx];
        const idx bit = static_cast<idx>(std::size_t{1} << (n - 1 - c));
        control_mask |= bit;
        if (shift[c_idx] == 0) {
            control_value |= bit;
        }
    }

    // Extract diagonal elements
    const Scalar a0 = A.coeff(0);
    const Scalar a1 = A.coeff(1);
    const Scalar a2 = A.coeff(2);
    const Scalar a3 = A.coeff(3);

#ifdef QPP_OPENMP
// NOLINTNEXTLINE
#pragma omp parallel for
#endif // QPP_OPENMP
    for (idx k00 = 0; k00 < D; ++k00) {
        // Skip indices that do not have both target bits == 0
        if ((static_cast<std::size_t>(k00) & static_cast<std::size_t>(s_i)) ||
            (static_cast<std::size_t>(k00) & static_cast<std::size_t>(s_j))) {
            continue;
        }

        // Control check
        if ((static_cast<std::size_t>(k00) &
             static_cast<std::size_t>(control_mask)) !=
            static_cast<std::size_t>(control_value)) {
            continue;
        }

        const idx k01 = k00 | s_j;
        const idx k10 = k00 | s_i;
        const idx k11 = k00 | s_i | s_j;

        // Apply diagonal scaling directly
        state.coeffRef(k00) *= a0;
        state.coeffRef(k01) *= a1;
        state.coeffRef(k10) *= a2;
        state.coeffRef(k11) *= a3;
    }
}

/**
 * @brief Applies the controlled 2-qubit diagonal gate \a A to the qubits \a i
 * and \a j of the multi-partite state vector \a state
 *
 * @param state Eigen expression
 * @param A Eigen expression (4 x 1 vector of diagonal elements)
 * @param ctrl Vector of control qubit indices
 * @param i Target subsystem index
 * @param j Target subsystem index
 * @param shift Vector of control values
 * @param n Number of qubits
 * @return Controlled gate \a A applied to the qubits \a i and \a j of \a state
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] expr_t<Derived1>
apply_ctrl_psi_2q_diag(const Eigen::MatrixBase<Derived1>& state,
                       const Eigen::MatrixBase<Derived2>& A,
                       const std::vector<idx>& ctrl, idx i, idx j,
                       const std::vector<idx>& shift, idx n) {
    // Functional deep copy
    expr_t<Derived1> result = state;

    // Apply transformation in-place on the copy
    apply_ctrl_psi_2q_diag_inplace(result, A, ctrl, i, j, shift, n);

    return result;
}

/**
 * @brief Applies the multi-qubit controlled diagonal gate \a A to the part
 * \a target of the multi-partite state vector \a state in-place
 *
 * @param state Eigen expression (modified in-place)
 * @param A Eigen expression (2^k x 1 vector of diagonal elements)
 * @param ctrl Vector of control qubit indices
 * @param target Subsystem indexes where the gate \a A is applied
 * @param shift Vector of control values (0: positive/|1> control, 1:
 * negative/|0> control)
 * @param n Number of qubits
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] void apply_ctrl_psi_kq_diag_inplace(
    Eigen::MatrixBase<Derived1>& state, const Eigen::MatrixBase<Derived2>& A,
    const std::vector<idx>& ctrl, const std::vector<idx>& target,
    const std::vector<idx>& shift, idx n) {

    const idx k = static_cast<idx>(target.size());
    const idx dim = static_cast<idx>(std::size_t{1} << k);
    const idx outer_dim = static_cast<idx>(std::size_t{1} << (n - k));

    // Input Validation
    assert(static_cast<idx>(A.size()) == dim &&
           "Gate A must have 2^k elements");
    assert(ctrl.size() == shift.size() &&
           "ctrl and shift vectors must have the same size");
    const idx ctrl_size = static_cast<idx>(ctrl.size());

#ifndef NDEBUG
    const idx D = static_cast<idx>(std::size_t{1} << n);
    for (idx t : target) {
        assert(t < n && "Target qubit index must be less than n");
    }
    assert(static_cast<idx>(state.size()) == D &&
           "State must be a 2^n x 1 vector");

    if (k > 1) {
        std::vector<idx> sorted_target = target;
        std::sort(sorted_target.begin(), sorted_target.end());
        assert(std::unique(sorted_target.begin(), sorted_target.end()) ==
                   sorted_target.end() &&
               "Target qubit indices must be distinct");
    }

    std::set<idx> target_set(target.begin(), target.end());
    for (idx c_idx = 0; c_idx < ctrl_size; ++c_idx) {
        const idx c = ctrl[c_idx];
        assert(c < n && target_set.find(c) == target_set.end() &&
               "Control/Target overlap");
        assert((shift[c_idx] == 0 || shift[c_idx] == 1) &&
               "Shift must be 0 or 1");
    }
#endif

    // Control Mask Logic
    idx expected_pattern_for_ones = 0;
    idx expected_zero_mask = 0;
    for (idx c_idx = 0; c_idx < ctrl_size; ++c_idx) {
        const idx bit =
            static_cast<idx>(std::size_t{1} << (n - 1 - ctrl[c_idx]));
        if (shift[c_idx] == 0) {
            expected_pattern_for_ones |= bit;
        } else {
            expected_zero_mask |= bit;
        }
    }

    // Index Mappings (Pre-compute the bitmasks for target configurations)
    std::vector<idx> inner_idx(dim);
    for (idx r = 0; r < dim; ++r) {
        idx index_r = 0;
        for (idx j = 0; j < k; ++j) {
            if ((static_cast<std::size_t>(r) >> (k - 1 - j)) & 1) {
                index_r |= (idx{1} << (n - 1 - target[j]));
            }
        }
        inner_idx[r] = index_r;
    }

    // Identify spectator qubits (qubits neither target nor control)
    std::vector<bool> is_target_or_ctrl(n, false);
    for (idx q : target) {
        is_target_or_ctrl[q] = true;
    }
    for (idx q : ctrl) {
        is_target_or_ctrl[q] = true;
    }

    std::vector<idx> spectator_qubits;
    for (idx i = 0; i < n; ++i) {
        if (!is_target_or_ctrl[i]) {
            spectator_qubits.push_back(i);
        }
    }

    // Pre-calculate outer indices for efficiency
    std::vector<idx> outer_idx(outer_dim);
    for (idx m = 0; m < outer_dim; ++m) {
        idx i_base = 0;
        for (idx j = 0; j < static_cast<idx>(spectator_qubits.size()); ++j) {
            if ((static_cast<std::size_t>(m) >> j) & 1) {
                i_base |= (idx{1} << (n - 1 - spectator_qubits[j]));
            }
        }
        outer_idx[m] = i_base | expected_pattern_for_ones;
    }

#ifdef QPP_OPENMP
// NOLINTNEXTLINE
#pragma omp parallel for
#endif // QPP_OPENMP
    for (idx m = 0; m < outer_dim; ++m) {
        const idx i_base = outer_idx[m];

        // Apply diagonal scaling directly across the target dimension
        for (idx r = 0; r < dim; ++r) {
            state(i_base | inner_idx[r]) *= A.coeff(r);
        }
    }
}

/**
 * @brief Applies the multi-qubit controlled diagonal gate \a A to the part \a
 * target of the multi-partite state vector \a state
 *
 * @param state Eigen expression
 * @param A Eigen expression (2^k x 1 vector of diagonal elements)
 * @param ctrl Vector of control qubit indices
 * @param target Subsystem indexes where the gate \a A is applied
 * @param shift Vector of control values
 * @param n Number of qubits
 * @return Controlled gate \a A applied to the part \a target of \a state
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] qpp::expr_t<Derived1> apply_ctrl_psi_kq_diag(
    const Eigen::MatrixBase<Derived1>& state,
    const Eigen::MatrixBase<Derived2>& A, const std::vector<idx>& ctrl,
    const std::vector<idx>& target, const std::vector<idx>& shift, idx n) {
    // Functional deep copy
    qpp::expr_t<Derived1> result = state;

    // Delegate to in-place implementation
    apply_ctrl_psi_kq_diag_inplace(result, A, ctrl, target, shift, n);

    return result;
}

/**
 * @brief Applies the controlled 1-qubit diagonal gate \a A to the qubit \a i of
 * the multi-partite density matrix \a state in-place
 *
 * @param state Eigen expression (modified in-place)
 * @param A Eigen expression (2 x 1 vector of diagonal elements)
 * @param i Target subsystem index
 * @param ctrl Vector of control qubit indices
 * @param shift Vector of control values (0: positive/|1> control, 1:
 * negative/|0> control)
 * @param n Number of qubits
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] void apply_ctrl_rho_1q_diag_inplace(
    Eigen::MatrixBase<Derived1>& state, const Eigen::MatrixBase<Derived2>& A,
    const std::vector<idx>& ctrl, idx i, const std::vector<idx>& shift, idx n) {
    using Scalar = typename Derived1::Scalar;
    const idx D = static_cast<idx>(std::size_t{1} << n);

    // Input Validation
    assert(i < n && "Target qubit index i must be less than n");
    assert(static_cast<idx>(state.rows()) == D &&
           static_cast<idx>(state.cols()) == D &&
           "State must be a square matrix sized 2^n x 2^n");
    assert(A.size() == 2 && "Diagonal gate A must have 2 elements");
    assert(ctrl.size() == shift.size() &&
           "ctrl and shift vectors must have the same size");
    const idx ctrl_size = static_cast<idx>(ctrl.size());

#ifndef NDEBUG
    for (idx c_idx = 0; c_idx < ctrl_size; ++c_idx) {
        const idx c = ctrl[c_idx];
        assert(c < n && "Control qubit index c must be less than n");
        assert(c != i && "Control qubit cannot also be a target qubit");
        assert((shift[c_idx] == 0 || shift[c_idx] == 1) &&
               "Shift value must be 0 or 1");
    }
#endif

    // Control Mask Logic
    idx expected_pattern_for_ones = 0;
    idx expected_zero_mask = 0;
    for (idx c_idx = 0; c_idx < ctrl_size; ++c_idx) {
        const idx bit =
            static_cast<idx>(std::size_t{1} << (n - 1 - ctrl[c_idx]));
        if (shift[c_idx] == 0) {
            expected_pattern_for_ones |= bit;
        } else {
            expected_zero_mask |= bit;
        }
    }

    auto control_is_met = [&](idx k_base) {
        return ((static_cast<std::size_t>(k_base) &
                 static_cast<std::size_t>(expected_pattern_for_ones)) ==
                static_cast<std::size_t>(expected_pattern_for_ones)) &&
               ((static_cast<std::size_t>(k_base) &
                 static_cast<std::size_t>(expected_zero_mask)) == 0);
    };

    // Indexing Constants
    const idx s_i = static_cast<idx>(std::size_t{1} << (n - 1 - i));
    const Scalar a0 = A.coeff(0);
    const Scalar a1 = A.coeff(1);
    const Scalar a0_conj = std::conj(a0);
    const Scalar a1_conj = std::conj(a1);

#ifdef QPP_OPENMP
// NOLINTNEXTLINE
#pragma omp parallel for collapse(2)
#endif // QPP_OPENMP
    for (idx r0 = 0; r0 < D; r0 += 2 * s_i) {
        for (idx c0 = 0; c0 < D; c0 += 2 * s_i) {
            for (idx r_rest = 0; r_rest < s_i; ++r_rest) {
                for (idx c_rest = 0; c_rest < s_i; ++c_rest) {
                    const idx r = r0 + r_rest;
                    const idx c = c0 + c_rest;

                    const idx r1 = r + s_i;
                    const idx c1 = c + s_i;

                    const bool row_ctrl = control_is_met(r);
                    const bool col_ctrl = control_is_met(c);

                    if (row_ctrl && col_ctrl) {
                        // Both row and column controlled: g_r * conj(g_c)
                        state.coeffRef(r, c) *= a0 * a0_conj;
                        state.coeffRef(r, c1) *= a0 * a1_conj;
                        state.coeffRef(r1, c) *= a1 * a0_conj;
                        state.coeffRef(r1, c1) *= a1 * a1_conj;
                    } else if (row_ctrl) {
                        // Only row controlled: g_r * 1
                        state.coeffRef(r, c) *= a0;
                        state.coeffRef(r, c1) *= a0;
                        state.coeffRef(r1, c) *= a1;
                        state.coeffRef(r1, c1) *= a1;
                    } else if (col_ctrl) {
                        // Only column controlled: 1 * conj(g_c)
                        state.coeffRef(r, c) *= a0_conj;
                        state.coeffRef(r, c1) *= a1_conj;
                        state.coeffRef(r1, c) *= a0_conj;
                        state.coeffRef(r1, c1) *= a1_conj;
                    }
                }
            }
        }
    }
}

/**
 * @brief Applies the controlled 1-qubit diagonal gate \a A to the qubit \a i of
 * the multi-partite density matrix \a state
 *
 * @param state Eigen expression
 * @param A Eigen expression (2 x 1 vector of diagonal elements)
 * @param i Target subsystem index
 * @param ctrl Vector of control qubit indices
 * @param shift Vector of control values
 * @param n Number of qubits
 * @return Controlled gate \a A applied to the qubits \a i of \a state
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] expr_t<Derived1>
apply_ctrl_rho_1q_diag(const Eigen::MatrixBase<Derived1>& state,
                       const Eigen::MatrixBase<Derived2>& A,
                       const std::vector<idx>& ctrl, idx i,
                       const std::vector<idx>& shift, idx n) {
    // Deep copy for functional interface
    expr_t<Derived1> result = state;

    // Apply transformation in-place
    apply_ctrl_rho_1q_diag_inplace(result, A, ctrl, i, shift, n);

    return result;
}

/**
 * @brief Applies the controlled 2-qubit diagonal gate \a A to the qubits \a i
 * and \a j of the multi-partite density matrix \a state in-place
 *
 * @param state Eigen expression (modified in-place)
 * @param A Eigen expression (4 x 1 vector of diagonal elements)
 * @param ctrl Vector of control qubit indices
 * @param i Target subsystem index
 * @param j Target subsystem index
 * @param shift Vector of control values (0: positive/|1> control, 1:
 * negative/|0> control)
 * @param n Number of qubits
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] void
apply_ctrl_rho_2q_diag_inplace(Eigen::MatrixBase<Derived1>& state,
                               const Eigen::MatrixBase<Derived2>& A,
                               const std::vector<idx>& ctrl, idx i, idx j,
                               const std::vector<idx>& shift, idx n) {
    using Scalar = typename Derived1::Scalar;
    const idx D [[maybe_unused]] = static_cast<idx>(std::size_t{1} << n);
    const idx ctrl_size = static_cast<idx>(ctrl.size());

    // Input Validation
    assert(i < n && j < n && i != j &&
           "Target qubits must be distinct and < n");
    assert(static_cast<idx>(state.rows()) == D &&
           static_cast<idx>(state.cols()) == D && "State must be 2^n x 2^n");
    assert(A.size() == 4 && "Diagonal gate A must have 4 elements");
    assert(ctrl.size() == shift.size() && "ctrl/shift size mismatch");

#ifndef NDEBUG
    std::set<idx> target_set = {i, j};
    for (idx c_idx = 0; c_idx < ctrl_size; ++c_idx) {
        const idx c = ctrl[c_idx];
        assert(c < n && target_set.find(c) == target_set.end() &&
               "Control/Target overlap");
    }
#endif

    const idx i_phys = n - i - 1;
    const idx j_phys = n - j - 1;
    const idx P_i = static_cast<idx>(std::size_t{1} << i_phys);
    const idx P_j = static_cast<idx>(std::size_t{1} << j_phys);
    const idx D_rest = static_cast<idx>(std::size_t{1} << (n - 2));

    // Pre-calculate diagonal elements and their conjugates
    const Scalar a[4] = {A.coeff(0), A.coeff(1), A.coeff(2), A.coeff(3)};
    const Scalar a_conj[4] = {std::conj(a[0]), std::conj(a[1]), std::conj(a[2]),
                              std::conj(a[3])};

    // Pre-calculate scaling matrix G(r, c) = a[r] * conj(a[c])
    Eigen::Matrix<Scalar, 4, 4> G;
    for (int r = 0; r < 4; ++r) {
        for (int col = 0; col < 4; ++col) {
            G(r, col) = a[r] * a_conj[col];
        }
    }

    // Control Mask Calculation
    idx expected_pattern_for_ones = 0;
    idx expected_zero_mask = 0;
    for (idx c_idx = 0; c_idx < ctrl_size; ++c_idx) {
        const idx bit =
            static_cast<idx>(std::size_t{1} << (n - 1 - ctrl[c_idx]));
        if (shift[c_idx] == 0) {
            expected_pattern_for_ones |= bit;
        } else {
            expected_zero_mask |= bit;
        }
    }

    auto control_is_met = [&](idx k_base) {
        return ((static_cast<std::size_t>(k_base) &
                 expected_pattern_for_ones) == expected_pattern_for_ones) &&
               ((static_cast<std::size_t>(k_base) & expected_zero_mask) == 0);
    };

#ifdef QPP_OPENMP
// NOLINTNEXTLINE
#pragma omp parallel for collapse(2)
#endif
    for (idx r = 0; r < D_rest; ++r) {
        for (idx c = 0; c < D_rest; ++c) {
            idx r_base = 0;
            idx c_base = 0;
            idx curr_r = r;
            idx curr_c = c;
            for (idx q = 0; q < n; ++q) {
                if (q != i_phys && q != j_phys) {
                    if (curr_r & 1) {
                        r_base |= (idx{1} << q);
                    }
                    if (curr_c & 1) {
                        c_base |= (idx{1} << q);
                    }
                    curr_r >>= 1;
                    curr_c >>= 1;
                }
            }

            const bool row_ctrl = control_is_met(r_base);
            const bool col_ctrl = control_is_met(c_base);

            if (!row_ctrl && !col_ctrl) {
                continue;
            }

            const idx rs[4] = {r_base, r_base + P_j, r_base + P_i,
                               r_base + P_i + P_j};
            const idx cs[4] = {c_base, c_base + P_j, c_base + P_i,
                               c_base + P_i + P_j};

            if (row_ctrl && col_ctrl) {
                for (int m = 0; m < 4; ++m) {
                    for (int n_idx = 0; n_idx < 4; ++n_idx) {
                        state.coeffRef(rs[m], cs[n_idx]) *= G(m, n_idx);
                    }
                }
            } else if (row_ctrl) {
                for (int m = 0; m < 4; ++m) {
                    for (int n_idx = 0; n_idx < 4; ++n_idx) {
                        state.coeffRef(rs[m], cs[n_idx]) *= a[m];
                    }
                }
            } else { // col_ctrl
                for (int m = 0; m < 4; ++m) {
                    for (int n_idx = 0; n_idx < 4; ++n_idx) {
                        state.coeffRef(rs[m], cs[n_idx]) *= a_conj[n_idx];
                    }
                }
            }
        }
    }
}

/**
 * @brief Applies the controlled 2-qubit diagonal gate \a A to the qubits \a i
 * and \a j of the multi-partite density matrix \a state
 *
 * @param state Eigen expression
 * @param A Eigen expression (4 x 1 vector of diagonal elements)
 * @param i Target subsystem index
 * @param j Target subsystem index
 * @param ctrl Vector of control qubit indices
 * @param shift Vector of control values
 * @param n Number of qubits
 * @return Controlled gate \a A applied to the qubits \a i and \a j of \a state
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] expr_t<Derived1>
apply_ctrl_rho_2q_diag(const Eigen::MatrixBase<Derived1>& state,
                       const Eigen::MatrixBase<Derived2>& A,
                       const std::vector<idx>& ctrl, idx i, idx j,
                       const std::vector<idx>& shift, idx n) {
    // Deep copy for the functional return
    expr_t<Derived1> result = state.derived();

    // Delegate to the in-place version
    apply_ctrl_rho_2q_diag_inplace(result, A, ctrl, i, j, shift, n);

    return result;
}

/**
 * @brief Applies the controlled multi-qubit diagonal gate \a A to the part \a
 * target of the multi-partite density matrix \a state in-place
 *
 * @param state Eigen expression (modified in-place)
 * @param A Eigen expression (2^k x 1 vector of diagonal elements)
 * @param target Subsystem indexes where the gate \a A is applied
 * @param ctrl Vector of control qubit indices
 * @param shift Vector of control values (0: positive/|1> control, 1:
 * negative/|0> control)
 * @param n Number of qubits
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] void apply_ctrl_rho_kq_diag_inplace(
    Eigen::MatrixBase<Derived1>& state, const Eigen::MatrixBase<Derived2>& A,
    const std::vector<idx>& ctrl, const std::vector<idx>& target,
    const std::vector<idx>& shift, idx n) {
    using Scalar = typename Derived1::Scalar;

    const idx k = static_cast<idx>(target.size());
    const idx D_k = static_cast<idx>(std::size_t{1} << k);
    const idx D [[maybe_unused]] = static_cast<idx>(std::size_t{1} << n);

    // 1. Validation
    assert(static_cast<idx>(A.size()) == D_k &&
           "Gate A must have 2^k elements");

    // 2. Pre-calculate Target relative offsets
    std::vector<idx> target_phys(k);
    for (idx l = 0; l < k; ++l) {
        target_phys[l] = n - target[l] - 1;
    }

    std::vector<idx> target_rel_idx(D_k, 0);
    for (idx m = 0; m < D_k; ++m) {
        for (idx l = 0; l < k; ++l) {
            if ((m >> (k - 1 - l)) & 1) {
                target_rel_idx[m] |= (idx{1} << target_phys[l]);
            }
        }
    }

    // 3. Pre-calculate Diagonal Scaling Matrix G(m, n) = A[m] * conj(A[n])
    std::vector<Scalar> a_conj(D_k);
    for (idx m = 0; m < D_k; ++m) {
        a_conj[m] = std::conj(A.coeff(m));
    }

    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> G(D_k, D_k);
    for (idx m = 0; m < D_k; ++m) {
        for (idx n_idx = 0; n_idx < D_k; ++n_idx) {
            G(m, n_idx) = A.coeff(m) * a_conj[n_idx];
        }
    }

    // 4. Pre-calculate Spectator (rest) bit patterns to avoid bit-twiddling in
    // loops
    std::vector<idx> rest_phys = fast_complement(target_phys, n);
    const idx D_rest = static_cast<idx>(std::size_t{1} << (n - k));
    std::vector<idx> rest_patterns(D_rest, 0);
    for (idx i = 0; i < D_rest; ++i) {
        for (idx q = 0; q < rest_phys.size(); ++q) {
            if ((i >> q) & 1) {
                rest_patterns[i] |= (idx{1} << rest_phys[q]);
            }
        }
    }

    // 5. Control Masking
    idx expected_pattern_for_ones = 0;
    idx expected_zero_mask = 0;
    for (idx c_idx = 0; c_idx < static_cast<idx>(ctrl.size()); ++c_idx) {
        const idx bit = static_cast<idx>(idx{1} << (n - 1 - ctrl[c_idx]));
        if (shift[c_idx] == 0) {
            expected_pattern_for_ones |= bit;
        } else {
            expected_zero_mask |= bit;
        }
    }

    auto control_is_met = [&](idx k_base) {
        return ((k_base & expected_pattern_for_ones) ==
                expected_pattern_for_ones) &&
               ((k_base & expected_zero_mask) == 0);
    };

    // 6. Execution Loop
#ifdef QPP_OPENMP
#pragma omp parallel for collapse(2)
#endif
    for (idx r = 0; r < D_rest; ++r) {
        for (idx c = 0; c < D_rest; ++c) {
            const idx r_base = rest_patterns[r];
            const idx c_base = rest_patterns[c];

            const bool row_ctrl = control_is_met(r_base);
            const bool col_ctrl = control_is_met(c_base);

            if (!row_ctrl && !col_ctrl) {
                continue;
            }

            if (row_ctrl && col_ctrl) {
                for (idx m = 0; m < D_k; ++m) {
                    const idx row_idx = r_base + target_rel_idx[m];
                    for (idx n_idx = 0; n_idx < D_k; ++n_idx) {
                        state.coeffRef(row_idx,
                                       c_base + target_rel_idx[n_idx]) *=
                            G(m, n_idx);
                    }
                }
            } else if (row_ctrl) {
                for (idx m = 0; m < D_k; ++m) {
                    const Scalar val_m = A.coeff(m);
                    const idx row_idx = r_base + target_rel_idx[m];
                    for (idx n_idx = 0; n_idx < D_k; ++n_idx) {
                        state.coeffRef(row_idx,
                                       c_base + target_rel_idx[n_idx]) *= val_m;
                    }
                }
            } else { // col_ctrl
                for (idx n_idx = 0; n_idx < D_k; ++n_idx) {
                    const Scalar val_n_conj = a_conj[n_idx];
                    const idx col_idx = c_base + target_rel_idx[n_idx];
                    for (idx m = 0; m < D_k; ++m) {
                        state.coeffRef(r_base + target_rel_idx[m], col_idx) *=
                            val_n_conj;
                    }
                }
            }
        }
    }
}

/**
 * @brief Applies the controlled multi-qubit diagonal gate \a A to the part \a
 * target of the multi-partite density matrix \a state
 *
 * @param state Eigen expression
 * @param A Eigen expression (2^k x 1 vector of diagonal elements)
 * @param target Subsystem indexes where the gate \a A is applied
 * @param ctrl Vector of control qubit indices
 * @param shift Vector of control values
 * @param n Number of qubits
 * @return Controlled gate \a A applied to the part \a target of \a state
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] expr_t<Derived1> apply_ctrl_rho_kq_diag(
    const Eigen::MatrixBase<Derived1>& state,
    const Eigen::MatrixBase<Derived2>& A, const std::vector<idx>& ctrl,
    const std::vector<idx>& target, const std::vector<idx>& shift, idx n) {
    // Functional deep copy
    expr_t<Derived1> result = state.derived();

    // Delegate to in-place implementation
    apply_ctrl_rho_kq_diag_inplace(result, A, ctrl, target, shift, n);

    return result;
}
} // namespace qpp::internal::kernels::qubit

#endif /* QPP_INTERNAL_KERNELS_QUBIT_APPLY_CTRL_HPP_ */
