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
 * @file qpp/internal/kernels/qudit/apply_ctrl_diag.hpp
 * @brief Internal highly optimized critical functions for qpp::applyCTRL_diag()
 */

#ifndef QPP_INTERNAL_KERNELS_QUDIT_APPLY_CTRL_DIAG_HPP_
#define QPP_INTERNAL_KERNELS_QUDIT_APPLY_CTRL_DIAG_HPP_

#include <cassert>
#include <vector>
#ifndef NDEBUG
#include <set>
#endif

#include <Eigen/Dense>

#include "qpp/types.hpp"

namespace qpp::internal::kernels::qudit {
/**
 * @brief Applies the controlled 1-qudit diagonal gate \a A to the qudit \a i of
 * the multi-partite state vector \a state in-place
 *
 * @param state Eigen expression (modified in-place)
 * @param A Eigen expression (2 x 1 vector of diagonal elements)
 * @param ctrl Vector of control qudit indices
 * @param i Target subsystem index
 * @param shift Vector of control shifts added to each control qudit value
 * modulo d (e.g., shifts the activating control configuration)
 * @param d Subsystem dimensions
 * @param n Number of qudits
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] void
apply_ctrl_psi_1q_diag_inplace(Eigen::MatrixBase<Derived1>& state,
                               const Eigen::MatrixBase<Derived2>& A,
                               const std::vector<idx>& ctrl, idx i,
                               const std::vector<idx>& shift, idx d, idx n) {

    // Precompute powers of d
    std::vector<idx> d_pow(n + 1, 1);
    for (idx i = 1; i <= n; ++i) {
        d_pow[i] = d_pow[i - 1] * d;
    }

    const idx ctrl_size = static_cast<idx>(ctrl.size());
    const idx spectator_count = n - 1 - ctrl_size; // k = 1

    // Total number of iterations for spectator configurations
    idx spectator_dim = 1;
    for (idx i = 0; i < spectator_count; ++i) {
        spectator_dim *= d;
    }

    // Input Validation
    assert(target.size() == 1 &&
           "This optimized function only supports 1 target qudit");
    assert(d >= 2 && "Qudit dimension must be at least 2");
    assert(static_cast<idx>(A.size()) == d &&
           "Gate A must have exactly d elements");
    assert(ctrl.size() == shift.size() &&
           "ctrl and shift vectors must have the same size");

#ifndef NDEBUG
    const idx D = d_pow[n];
    assert(target[0] < n && "Target qudit index must be less than n");
    assert(static_cast<idx>(state.size()) == D && "State size mismatch");

    std::set<idx> target_set(target.begin(), target.end());
    for (idx c_idx = 0; c_idx < ctrl_size; ++c_idx) {
        const idx c = ctrl[c_idx];
        assert(c < n && target_set.find(c) == target_set.end() &&
               "Control/Target overlap");
        assert(shift[c_idx] < d && "Shift must be less than d");
    }
#endif

    // Identify spectator qudits
    std::vector<bool> is_target_or_ctrl(n, false);
    is_target_or_ctrl[i] = true;
    for (idx q : ctrl) {
        is_target_or_ctrl[q] = true;
    }

    std::vector<idx> spectator_qubits;
    spectator_qubits.reserve(spectator_count);
    for (idx i = 0; i < n; ++i) {
        if (!is_target_or_ctrl[i]) {
            spectator_qubits.push_back(i);
        }
    }

    idx ctrl_dim = 1;
    for (idx i = 0; i < ctrl_size; ++i) {
        ctrl_dim *= d;
    }

    const idx target_stride = d_pow[n - 1 - i];

    // Outer loop splits across all spectator and control combinations
#ifdef QPP_OPENMP
#pragma omp parallel for collapse(2)
#endif // QPP_OPENMP
    for (idx m = 0; m < spectator_dim; ++m) {
        for (idx c_val = 0; c_val < ctrl_dim; ++c_val) {

            idx i_base = 0;

            // 1. Map spectator configurations to their global state index
            // positions
            idx temp_m = m;
            for (idx j = 0; j < spectator_count; ++j) {
                idx q_val = temp_m % d;
                temp_m /= d;
                i_base += q_val * d_pow[n - 1 - spectator_qubits[j]];
            }

            // 2. Map control configurations left-to-right to match global
            // indexing
            idx temp_c = c_val;
            for (idx j = ctrl_size; j-- > 0;) {
                idx q_val = temp_c % d;
                temp_c /= d;
                i_base += q_val * d_pow[n - 1 - ctrl[j]];
            }

            // Reconstruct the actual collective integer value 'c' of the
            // controls
            idx joint_ctrl_power = 0;
            for (idx j = 0; j < ctrl_size; ++j) {
                idx q_val = (c_val / d_pow[ctrl_size - 1 - j]) % d;
                idx effective_val = (q_val + shift[j]) % d;
                joint_ctrl_power = (joint_ctrl_power * d) + effective_val;
            }

            // 3 & 4. Directly update the d coefficients of the target qudit
            // without heap allocation for A_powered.
            for (idx r = 0; r < d; ++r) {
                // target_stride * r maps directly to the global index offset
                // for target qudit state 'r'
                state(i_base + (r * target_stride)) *=
                    std::pow(A.coeff(r), joint_ctrl_power);
            }
        }
    }
}

/**
 * @brief Applies the controlled 1-qudit diagonal gate \a A to the qudit \a i of
 * the multi-partite state vector \a state
 *
 * @param state Eigen expression
 * @param A Eigen expression (2 x 1 vector of diagonal elements)
 * @param i Target subsystem index
 * @param ctrl Vector of control qudit indices
 * @param shift Vector of control shifts added to each control qudit value
 * modulo d (e.g., shifts the activating control configuration)
 * @param n Number of qudits
 * @return Controlled gate \a A applied to the qudit \a i of \a state
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
 * @brief Applies the controlled 2-qudit diagonal gate \a A to the qudits \a i
 * and \a j of the multi-partite state vector \a state in-place
 *
 * @param state Eigen expression (modified in-place)
 * @param A Eigen expression (4 x 1 vector of diagonal elements)
 * @param ctrl Vector of control qudit indices
 * @param i Target subsystem index
 * @param j Target subsystem index
 * @param shift Vector of control shifts added to each control qudit value
 * modulo d (e.g., shifts the activating control configuration)
 * @param d Subsystem dimensions
 * @param n Number of qudits
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] void
apply_ctrl_psi_2q_diag_inplace(Eigen::MatrixBase<Derived1>& state,
                               const Eigen::MatrixBase<Derived2>& A,
                               const std::vector<idx>& ctrl, idx i, idx j,
                               const std::vector<idx>& shift, idx d, idx n) {
    // Precompute powers of d
    std::vector<idx> d_pow(n + 1, 1);
    for (idx q = 1; q <= n; ++q) {
        d_pow[q] = d_pow[q - 1] * d;
    }

    [[maybe_unused]] const idx dimA = d * d; // d^2 since k = 2
    const idx ctrl_size = static_cast<idx>(ctrl.size());
    const idx spectator_count = n - 2 - ctrl_size;

    // Total number of iterations for spectator configurations
    idx spectator_dim = 1;
    for (idx q = 0; q < spectator_count; ++q) {
        spectator_dim *= d;
    }

    // Input Validation
    assert(d >= 2 && "Qudit dimension must be at least 2");
    assert(static_cast<idx>(A.size()) == dimA &&
           "Gate A must have exactly d^2 elements");
    assert(ctrl.size() == shift.size() &&
           "ctrl and shift vectors must have the same size");
    assert(i != j && "Target qudit indices must be distinct");

#ifndef NDEBUG
    const idx D = d_pow[n];
    assert(i < n && j < n && "Target qudit indices must be less than n");
    assert(static_cast<idx>(state.size()) == D && "State size mismatch");

    std::set<idx> target_set{i, j};
    for (idx c_idx = 0; c_idx < ctrl_size; ++c_idx) {
        const idx c = ctrl[c_idx];
        assert(c < n && target_set.find(c) == target_set.end() &&
               "Control/Target overlap");
        assert(shift[c_idx] < d && "Shift must be less than d");
    }
#endif

    // Identify spectator qudits
    std::vector<bool> is_target_or_ctrl(n, false);
    is_target_or_ctrl[i] = true;
    is_target_or_ctrl[j] = true;
    for (idx q : ctrl) {
        is_target_or_ctrl[q] = true;
    }

    std::vector<idx> spectator_qubits;
    spectator_qubits.reserve(spectator_count);
    for (idx q = 0; q < n; ++q) {
        if (!is_target_or_ctrl[q]) {
            spectator_qubits.push_back(q);
        }
    }

    idx ctrl_dim = 1;
    for (idx q = 0; q < ctrl_size; ++q) {
        ctrl_dim *= d;
    }

    // Precompute index strides for target qudits i and j
    // Big-endian mapping consistent with qpp layout: n - 1 - index
    const idx stride_i = d_pow[n - 1 - i];
    const idx stride_j = d_pow[n - 1 - j];

    // Outer loop splits across all spectator and control combinations
#ifdef QPP_OPENMP
#pragma omp parallel for collapse(2)
#endif // QPP_OPENMP
    for (idx m = 0; m < spectator_dim; ++m) {
        for (idx c_val = 0; c_val < ctrl_dim; ++c_val) {

            idx i_base = 0;

            // 1. Map spectator configurations to their global state index
            // positions
            idx temp_m = m;
            for (idx q = 0; q < spectator_count; ++q) {
                idx q_val = temp_m % d;
                temp_m /= d;
                i_base += q_val * d_pow[n - 1 - spectator_qubits[q]];
            }

            // 2. Map control configurations left-to-right to match global
            // indexing
            idx temp_c = c_val;
            for (idx q = ctrl_size; q-- > 0;) {
                idx q_val = temp_c % d;
                temp_c /= d;
                i_base += q_val * d_pow[n - 1 - ctrl[q]];
            }

            // Reconstruct the actual collective integer value 'c' of the
            // controls
            idx joint_ctrl_power = 0;
            for (idx q = 0; q < ctrl_size; ++q) {
                idx q_val = (c_val / d_pow[ctrl_size - 1 - q]) % d;
                idx effective_val = (q_val + shift[q]) % d;
                joint_ctrl_power = (joint_ctrl_power * d) + effective_val;
            }

            // 3 & 4. Unrolled execution across the d^2 element dimensions of
            // targets i and j
            idx r = 0; // Tracks flat index inside matrix A
            for (idx v_i = 0; v_i < d; ++v_i) {
                const idx offset_i = v_i * stride_i;
                for (idx v_j = 0; v_j < d; ++v_j) {
                    const idx offset_j = v_j * stride_j;

                    // State vector position for target states v_i and v_j
                    state(i_base + offset_i + offset_j) *=
                        std::pow(A.coeff(r), joint_ctrl_power);
                    ++r;
                }
            }
        }
    }
}

/**
 * @brief Applies the controlled 2-qudit diagonal gate \a A to the qudits \a i
 * and \a j of the multi-partite state vector \a state
 *
 * @param state Eigen expression
 * @param A Eigen expression (4 x 1 vector of diagonal elements)
 * @param ctrl Vector of control qudit indices
 * @param i Target subsystem index
 * @param j Target subsystem index
 * @param shift Vector of control shifts added to each control qudit value
 * modulo d (e.g., shifts the activating control configuration)
 * @param n Number of qudits
 * @return Controlled gate \a A applied to the qudits \a i and \a j of \a state
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
 * @brief Applies the multi-qudit controlled diagonal gate A^c to the part
 * \a target of the multi-partite state vector \a state in-place.
 * The power 'c' is determined by the collective base-d value of the control
 * qudits.
 *
 * @param state Eigen expression (modified in-place)
 * @param A Eigen expression (d^k x 1 vector of diagonal elements)
 * @param ctrl Vector of control qudit indices
 * @param target Subsystem indexes where the gate \a A is applied
 * @param shift Vector of control shifts added to each control qudit value
 * modulo d (e.g., shifts the activating control configuration)
 * @param d Subsystem dimensions
 * @param n Number of qudits
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] void apply_ctrl_psi_kq_diag_inplace(
    Eigen::MatrixBase<Derived1>& state, const Eigen::MatrixBase<Derived2>& A,
    const std::vector<idx>& ctrl, const std::vector<idx>& target,
    const std::vector<idx>& shift, idx d, idx n) {

    const idx k = static_cast<idx>(target.size());

    // Precompute powers of d
    std::vector<idx> d_pow(n + 1, 1);
    for (idx i = 1; i <= n; ++i) {
        d_pow[i] = d_pow[i - 1] * d;
    }

    const idx dimA = d_pow[k]; // d^k
    const idx ctrl_size = static_cast<idx>(ctrl.size());
    const idx spectator_count = n - k - ctrl_size;

    // Total number of iterations for spectator configurations
    idx spectator_dim = 1;
    for (idx i = 0; i < spectator_count; ++i) {
        spectator_dim *= d;
    }

    // Input Validation
    assert(d >= 2 && "Qudit dimension must be at least 2");
    assert(static_cast<idx>(A.size()) == dimA &&
           "Gate A must have d^k elements");
    assert(ctrl.size() == shift.size() &&
           "ctrl and shift vectors must have the same size");

#ifndef NDEBUG
    const idx D = d_pow[n];
    for (idx t : target) {
        assert(t < n && "Target qudit index must be less than n");
    }
    assert(static_cast<idx>(state.size()) == D && "State size mismatch");

    if (k > 1) {
        std::vector<idx> sorted_target = target;
        std::sort(sorted_target.begin(), sorted_target.end());
        assert(std::unique(sorted_target.begin(), sorted_target.end()) ==
                   sorted_target.end() &&
               "Target qudit indices must be distinct");
    }

    std::set<idx> target_set(target.begin(), target.end());
    for (idx c_idx = 0; c_idx < ctrl_size; ++c_idx) {
        const idx c = ctrl[c_idx];
        assert(c < n && target_set.find(c) == target_set.end() &&
               "Control/Target overlap");
        assert(shift[c_idx] < d && "Shift must be less than d");
    }
#endif

    // Identify spectator qudits
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

    // We must iterate over all possible values of control qudits
    idx ctrl_dim = 1;
    for (idx i = 0; i < ctrl_size; ++i) {
        ctrl_dim *= d;
    }

    // ... [Keep the initial setup, validations, and spectator calculations the
    // same] ...

    // Outer loop splits across all spectator and control combinations
#ifdef QPP_OPENMP
#pragma omp parallel for collapse(2)
#endif // QPP_OPENMP
    for (idx m = 0; m < spectator_dim; ++m) {
        for (idx c_val = 0; c_val < ctrl_dim; ++c_val) {

            idx i_base = 0;

            // 1. Map spectator configurations to their global state index
            // positions
            idx temp_m = m;
            for (idx j = 0; j < static_cast<idx>(spectator_qubits.size());
                 ++j) {
                idx q_val = temp_m % d;
                temp_m /= d;
                i_base += q_val * d_pow[n - 1 - spectator_qubits[j]];
            }

            // 2. Map control configurations left-to-right to match global
            // indexing
            idx temp_c = c_val;
            idx joint_ctrl_power = 0;

            // We unpack temp_c from right to left, but assign it to controls
            // from right to left to preserve standard positional significance.
            for (idx j = ctrl_size; j-- > 0;) {
                idx q_val = temp_c % d;
                temp_c /= d;

                // Add to base index placement
                i_base += q_val * d_pow[n - 1 - ctrl[j]];
            }

            // Reconstruct the actual collective integer value 'c' of the
            // controls based on their literal appearance order in the ctrl
            // vector
            for (idx j = 0; j < ctrl_size; ++j) {
                // Extract what q_val actually was for ctrl[j]
                idx q_val = (c_val / d_pow[ctrl_size - 1 - j]) % d;
                idx effective_val = (q_val + shift[j]) % d;
                joint_ctrl_power = (joint_ctrl_power * d) + effective_val;
            }

            // 3. Pre-extract the component-wise powers of the diagonal matrix
            // elements
            dyn_col_vect<typename Derived2::Scalar> A_powered(dimA);
            for (idx r = 0; r < dimA; ++r) {
                A_powered(r) = std::pow(A.coeff(r), joint_ctrl_power);
            }

            // 4. Apply the scaled diagonal elements directly across the target
            // dimensions
            for (idx r = 0; r < dimA; ++r) {
                idx index_r = 0;
                idx temp_r = r;
                for (idx j = k; j-- > 0;) {
                    idx q_val = temp_r % d;
                    temp_r /= d;
                    index_r += q_val * d_pow[n - 1 - target[j]];
                }

                state(i_base + index_r) *= A_powered(r);
            }
        }
    }
}

/**
 * @brief Applies the multi-qudit controlled diagonal gate \a A to the part \a
 * target of the multi-partite state vector \a state
 *
 * @param state Eigen expression
 * @param A Eigen expression (2^k x 1 vector of diagonal elements)
 * @param ctrl Vector of control qudit indices
 * @param target Subsystem indexes where the gate \a A is applied
 * @param shift Vector of control shifts added to each control qudit value
 * modulo d (e.g., shifts the activating control configuration)
 * @param n Number of qudits
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
 * @brief Applies the controlled 1-qudit diagonal gate \a A to the qudit \a i of
 * the multi-partite density matrix \a state in-place
 *
 * @param state Eigen expression (modified in-place)
 * @param A Eigen expression (2 x 1 vector of diagonal elements)
 * @param ctrl Vector of control qudit indices
 * @param i Target subsystem index
 * @param shift Vector of control shifts added to each control qudit value
 * modulo d (e.g., shifts the activating control configuration)
 * @param d Subsystem dimensions
 * @param n Number of qudits
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] void
apply_ctrl_rho_1q_diag_inplace(Eigen::MatrixBase<Derived1>& state,
                               const Eigen::MatrixBase<Derived2>& A,
                               const std::vector<idx>& ctrl, idx i,
                               const std::vector<idx>& shift, idx d, idx n) {

    // Precompute powers of d
    std::vector<idx> d_pow(n + 1, 1);
    for (idx i = 1; i <= n; ++i) {
        d_pow[i] = d_pow[i - 1] * d;
    }

    const idx ctrl_size = static_cast<idx>(ctrl.size());
    const idx spectator_count = n - 1 - ctrl_size; // k = 1

    idx spectator_dim = 1;
    for (idx i = 0; i < spectator_count; ++i) {
        spectator_dim *= d;
    }

    idx ctrl_dim = 1;
    for (idx i = 0; i < ctrl_size; ++i) {
        ctrl_dim *= d;
    }

    // Input Validation
    assert(target.size() == 1 &&
           "This optimized function only supports 1 target qudit");
    assert(d >= 2 && "Qudit dimension must be at least 2");
    assert(static_cast<idx>(A.size()) == d &&
           "Gate A must have exactly d elements");
    assert(ctrl.size() == shift.size() &&
           "ctrl and shift vectors must have the same size");

#ifndef NDEBUG
    const idx D = d_pow[n];
    assert(target[0] < n && "Target qudit index must be less than n");
    assert(static_cast<idx>(state.rows()) == D &&
           static_cast<idx>(state.cols()) == D &&
           "State must be a d^n x d^n matrix");

    std::set<idx> target_set(target.begin(), target.end());
    for (idx c_idx = 0; c_idx < ctrl_size; ++c_idx) {
        const idx c = ctrl[c_idx];
        assert(c < n && target_set.find(c) == target_set.end() &&
               "Control/Target overlap");
        assert(shift[c_idx] < d && "Shift must be less than d");
    }
#endif

    // Identify spectator qudits
    std::vector<bool> is_target_or_ctrl(n, false);
    is_target_or_ctrl[i] = true;
    for (idx q : ctrl) {
        is_target_or_ctrl[q] = true;
    }

    std::vector<idx> spectator_qubits;
    spectator_qubits.reserve(spectator_count);
    for (idx i = 0; i < n; ++i) {
        if (!is_target_or_ctrl[i]) {
            spectator_qubits.push_back(i);
        }
    }

    // Precalculate powers of diagonal elements for all possible control gate
    // modifications Because k = 1, dimA == d
    dyn_mat<typename Derived2::Scalar> A_powers(d, ctrl_dim);
    for (idx c_val = 0; c_val < ctrl_dim; ++c_val) {
        idx joint_ctrl_power = 0;
        for (idx j = 0; j < ctrl_size; ++j) {
            idx q_val = (c_val / d_pow[ctrl_size - 1 - j]) % d;
            idx effective_val = (q_val + shift[j]) % d;
            joint_ctrl_power = (joint_ctrl_power * d) + effective_val;
        }
        for (idx r = 0; r < d; ++r) {
            A_powers(r, c_val) = std::pow(A.coeff(r), joint_ctrl_power);
        }
    }

    // Precompute the single stride needed for the target qudit
    const idx target_stride = d_pow[n - 1 - i];

    // Loop over row configurations (m_row, c_row) and column configurations
    // (m_col, c_col)
#ifdef QPP_OPENMP
#pragma omp parallel for collapse(2)
#endif // QPP_OPENMP
    for (idx m_row = 0; m_row < spectator_dim; ++m_row) {
        for (idx c_row = 0; c_row < ctrl_dim; ++c_row) {

            // 1. Compute row base global index
            idx row_base = 0;
            idx temp_m = m_row;
            for (idx j = 0; j < spectator_count; ++j) {
                row_base += (temp_m % d) * d_pow[n - 1 - spectator_qubits[j]];
                temp_m /= d;
            }
            idx temp_c = c_row;
            for (idx j = ctrl_size; j-- > 0;) {
                row_base += (temp_c % d) * d_pow[n - 1 - ctrl[j]];
                temp_c /= d;
            }

            // 2. Inner loop over columns to execute the 2D matrix
            // transformation
            for (idx m_col = 0; m_col < spectator_dim; ++m_col) {
                for (idx c_col = 0; c_col < ctrl_dim; ++c_col) {

                    idx col_base = 0;
                    temp_m = m_col;
                    for (idx j = 0; j < spectator_count; ++j) {
                        col_base +=
                            (temp_m % d) * d_pow[n - 1 - spectator_qubits[j]];
                        temp_m /= d;
                    }
                    temp_c = c_col;
                    for (idx j = ctrl_size; j-- > 0;) {
                        col_base += (temp_c % d) * d_pow[n - 1 - ctrl[j]];
                        temp_c /= d;
                    }

                    // 3. Flattened Target Element Application: state(final_row,
                    // final_col) *= U * conj(U) The dynamic nested loops over
                    // r_row and r_col are unrolled into flat stride offsets
                    for (idx r_row = 0; r_row < d; ++r_row) {
                        const std::complex<double> u_row =
                            A_powers(r_row, c_row);
                        const idx final_row =
                            row_base + (r_row * target_stride);

                        for (idx r_col = 0; r_col < d; ++r_col) {
                            const std::complex<double> u_col_conj =
                                std::conj(A_powers(r_col, c_col));
                            const idx final_col =
                                col_base + (r_col * target_stride);

                            state(final_row, final_col) *= (u_row * u_col_conj);
                        }
                    }
                }
            }
        }
    }
}

/**
 * @brief Applies the controlled 1-qudit diagonal gate \a A to the qudit \a i of
 * the multi-partite density matrix \a state
 *
 * @param state Eigen expression
 * @param A Eigen expression (2 x 1 vector of diagonal elements)
 * @param i Target subsystem index
 * @param ctrl Vector of control qudit indices
 * @param shift Vector of control shifts added to each control qudit value
 * modulo d (e.g., shifts the activating control configuration)
 * @param n Number of qudits
 * @return Controlled gate \a A applied to the qudits \a i of \a state
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
 * @brief Applies the controlled 2-qudit diagonal gate \a A to the qudits \a i
 * and \a j of the multi-partite density matrix \a state in-place
 *
 * @param state Eigen expression (modified in-place)
 * @param A Eigen expression (4 x 1 vector of diagonal elements)
 * @param ctrl Vector of control qudit indices
 * @param i Target subsystem index
 * @param j Target subsystem index
 * @param shift Vector of control shifts added to each control qudit value
 * modulo d (e.g., shifts the activating control configuration)
 * @param d Subsystem dimensions
 * @param n Number of qudits
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] void
apply_ctrl_rho_2q_diag_inplace(Eigen::MatrixBase<Derived1>& state,
                               const Eigen::MatrixBase<Derived2>& A,
                               const std::vector<idx>& ctrl, idx i, idx j,
                               const std::vector<idx>& shift, idx d, idx n) {
    // Precompute powers of d
    std::vector<idx> d_pow(n + 1, 1);
    for (idx q = 1; q <= n; ++q) {
        d_pow[q] = d_pow[q - 1] * d;
    }

    const idx dimA = d * d; // d^2 since k = 2
    const idx ctrl_size = static_cast<idx>(ctrl.size());
    const idx spectator_count = n - 2 - ctrl_size;

    idx spectator_dim = 1;
    for (idx q = 0; q < spectator_count; ++q) {
        spectator_dim *= d;
    }

    idx ctrl_dim = 1;
    for (idx q = 0; q < ctrl_size; ++q) {
        ctrl_dim *= d;
    }

    // Input Validation
    assert(d >= 2 && "Qudit dimension must be at least 2");
    assert(static_cast<idx>(A.size()) == dimA &&
           "Gate A must have exactly d^2 elements");
    assert(ctrl.size() == shift.size() &&
           "ctrl and shift vectors must have the same size");
    assert(i != j && "Target qudit indices must be distinct");

#ifndef NDEBUG
    const idx D = d_pow[n];
    assert(i < n && j < n && "Target qudit indices must be less than n");
    assert(static_cast<idx>(state.rows()) == D &&
           static_cast<idx>(state.cols()) == D &&
           "State must be a d^n x d^n matrix");

    std::set<idx> target_set{i, j};
    for (idx c_idx = 0; c_idx < ctrl_size; ++c_idx) {
        const idx c = ctrl[c_idx];
        assert(c < n && target_set.find(c) == target_set.end() &&
               "Control/Target overlap");
        assert(shift[c_idx] < d && "Shift must be less than d");
    }
#endif

    // Identify spectator qudits
    std::vector<bool> is_target_or_ctrl(n, false);
    is_target_or_ctrl[i] = true;
    is_target_or_ctrl[j] = true;
    for (idx q : ctrl) {
        is_target_or_ctrl[q] = true;
    }

    std::vector<idx> spectator_qubits;
    spectator_qubits.reserve(spectator_count);
    for (idx q = 0; q < n; ++q) {
        if (!is_target_or_ctrl[q]) {
            spectator_qubits.push_back(q);
        }
    }

    // Precalculate powers of diagonal elements for all possible control gate
    // configurations
    dyn_mat<typename Derived2::Scalar> A_powers(dimA, ctrl_dim);
    for (idx c_val = 0; c_val < ctrl_dim; ++c_val) {
        idx joint_ctrl_power = 0;
        for (idx q = 0; q < ctrl_size; ++q) {
            idx q_val = (c_val / d_pow[ctrl_size - 1 - q]) % d;
            idx effective_val = (q_val + shift[q]) % d;
            joint_ctrl_power = (joint_ctrl_power * d) + effective_val;
        }
        for (idx r = 0; r < dimA; ++r) {
            A_powers(r, c_val) = std::pow(A.coeff(r), joint_ctrl_power);
        }
    }

    // Precompute index strides for target qudits i and j based on big-endian
    // layout
    const idx stride_i = d_pow[n - 1 - i];
    const idx stride_j = d_pow[n - 1 - j];

    // Loop over row configurations (m_row, c_row) and column configurations
    // (m_col, c_col)
#ifdef QPP_OPENMP
#pragma omp parallel for collapse(2)
#endif // QPP_OPENMP
    for (idx m_row = 0; m_row < spectator_dim; ++m_row) {
        for (idx c_row = 0; c_row < ctrl_dim; ++c_row) {

            // 1. Compute row base global index
            idx row_base = 0;
            idx temp_m = m_row;
            for (idx q = 0; q < spectator_count; ++q) {
                row_base += (temp_m % d) * d_pow[n - 1 - spectator_qubits[q]];
                temp_m /= d;
            }
            idx temp_c = c_row;
            for (idx q = ctrl_size; q-- > 0;) {
                row_base += (temp_c % d) * d_pow[n - 1 - ctrl[q]];
                temp_c /= d;
            }

            // 2. Inner loop over columns to execute the 2D matrix
            // transformation
            for (idx m_col = 0; m_col < spectator_dim; ++m_col) {
                for (idx c_col = 0; c_col < ctrl_dim; ++c_col) {

                    idx col_base = 0;
                    temp_m = m_col;
                    for (idx q = 0; q < spectator_count; ++q) {
                        col_base +=
                            (temp_m % d) * d_pow[n - 1 - spectator_qubits[q]];
                        temp_m /= d;
                    }
                    temp_c = c_col;
                    for (idx q = ctrl_size; q-- > 0;) {
                        col_base += (temp_c % d) * d_pow[n - 1 - ctrl[q]];
                        temp_c /= d;
                    }

                    // 3. Flat Target Mapping over the 2-qudit subsystems (d x d
                    // combinations)
                    idx r_row = 0;
                    for (idx v_row_i = 0; v_row_i < d; ++v_row_i) {
                        const idx offset_row_i = v_row_i * stride_i;
                        for (idx v_row_j = 0; v_row_j < d; ++v_row_j) {

                            const idx final_row =
                                row_base + offset_row_i + (v_row_j * stride_j);
                            const std::complex<double> u_row =
                                A_powers(r_row, c_row);

                            idx r_col = 0;
                            for (idx v_col_i = 0; v_col_i < d; ++v_col_i) {
                                const idx offset_col_i = v_col_i * stride_i;
                                for (idx v_col_j = 0; v_col_j < d; ++v_col_j) {

                                    const idx final_col = col_base +
                                                          offset_col_i +
                                                          (v_col_j * stride_j);
                                    const std::complex<double> u_col_conj =
                                        std::conj(A_powers(r_col, c_col));

                                    state(final_row, final_col) *=
                                        (u_row * u_col_conj);
                                    ++r_col;
                                }
                            }
                            ++r_row;
                        }
                    }
                }
            }
        }
    }
}

/**
 * @brief Applies the controlled 2-qudit diagonal gate \a A to the qudits \a i
 * and \a j of the multi-partite density matrix \a state
 *
 * @param state Eigen expression
 * @param A Eigen expression (4 x 1 vector of diagonal elements)
 * @param i Target subsystem index
 * @param j Target subsystem index
 * @param ctrl Vector of control qudit indices
 * @param shift Vector of control shifts added to each control qudit value
 * modulo d (e.g., shifts the activating control configuration)
 * @param n Number of qudits
 * @return Controlled gate \a A applied to the qudits \a i and \a j of \a state
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
 * @brief Applies the multi-qudit controlled diagonal gate A^c to the part
 * \a target of the multi-partite density matrix \a state in-place.
 * The transformation is \rho -> U * \rho * U^\dagger.
 *
 * @param state Eigen expression (modified in-place, size d^n x d^n)
 * @param A Eigen expression (d^k x 1 vector of diagonal elements)
 * @param ctrl Vector of control qudit indices
 * @param target Subsystem indexes where the gate \a A is applied
 * @param shift Vector of control shifts added to each control qudit value
 * modulo d (e.g., shifts the activating control configuration)
 * @param d The dimension of the qudits (d >= 2)
 * @param n Number of qudits
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] void apply_ctrl_rho_kq_diag_inplace(
    Eigen::MatrixBase<Derived1>& state, const Eigen::MatrixBase<Derived2>& A,
    const std::vector<idx>& ctrl, const std::vector<idx>& target,
    const std::vector<idx>& shift, idx d, idx n) {

    const idx k = static_cast<idx>(target.size());

    // Precompute powers of d
    std::vector<idx> d_pow(n + 1, 1);
    for (idx i = 1; i <= n; ++i) {
        d_pow[i] = d_pow[i - 1] * d;
    }

    const idx dimA = d_pow[k]; // d^k
    const idx ctrl_size = static_cast<idx>(ctrl.size());
    const idx spectator_count = n - k - ctrl_size;

    idx spectator_dim = 1;
    for (idx i = 0; i < spectator_count; ++i) {
        spectator_dim *= d;
    }

    idx ctrl_dim = 1;
    for (idx i = 0; i < ctrl_size; ++i) {
        ctrl_dim *= d;
    }

    // Input Validation
    assert(d >= 2 && "Qudit dimension must be at least 2");
    assert(static_cast<idx>(A.size()) == dimA &&
           "Gate A must have d^k elements");
    assert(ctrl.size() == shift.size() &&
           "ctrl and shift vectors must have the same size");

#ifndef NDEBUG
    const idx D = d_pow[n];
    for (idx t : target) {
        assert(t < n && "Target qudit index must be less than n");
    }
    assert(static_cast<idx>(state.rows()) == D &&
           static_cast<idx>(state.cols()) == D &&
           "State must be a d^n x d^n matrix");

    std::set<idx> target_set(target.begin(), target.end());
    for (idx c_idx = 0; c_idx < ctrl_size; ++c_idx) {
        const idx c = ctrl[c_idx];
        assert(c < n && target_set.find(c) == target_set.end() &&
               "Control/Target overlap");
        assert(shift[c_idx] < d && "Shift must be less than d");
    }
#endif

    // Identify spectator qudits
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

    // Precalculate target configurations index displacements
    std::vector<idx> target_displacements(dimA);
    for (idx r = 0; r < dimA; ++r) {
        idx index_r = 0;
        idx temp_r = r;
        for (idx j = k; j-- > 0;) {
            idx q_val = temp_r % d;
            temp_r /= d;
            index_r += q_val * d_pow[n - 1 - target[j]];
        }
        target_displacements[r] = index_r;
    }

    // Precalculate powers of diagonal elements for all possible control gate
    // modifications Precomputing avoids repeating std::pow inside the nested
    // OMP loops
    dyn_mat<typename Derived2::Scalar> A_powers(dimA, ctrl_dim);
    for (idx c_val = 0; c_val < ctrl_dim; ++c_val) {
        idx joint_ctrl_power = 0;
        for (idx j = 0; j < ctrl_size; ++j) {
            idx q_val = (c_val / d_pow[ctrl_size - 1 - j]) % d;
            idx effective_val = (q_val + shift[j]) % d;
            joint_ctrl_power = (joint_ctrl_power * d) + effective_val;
        }
        for (idx r = 0; r < dimA; ++r) {
            A_powers(r, c_val) = std::pow(A.coeff(r), joint_ctrl_power);
        }
    }

    // Loop over row configurations (m_row, c_row) and column configurations
    // (m_col, c_col)
#ifdef QPP_OPENMP
#pragma omp parallel for collapse(2)
#endif // QPP_OPENMP
    for (idx m_row = 0; m_row < spectator_dim; ++m_row) {
        for (idx c_row = 0; c_row < ctrl_dim; ++c_row) {

            // 1. Compute row base global index
            idx row_base = 0;
            idx temp_m = m_row;
            for (idx j = 0; j < static_cast<idx>(spectator_qubits.size());
                 ++j) {
                row_base += (temp_m % d) * d_pow[n - 1 - spectator_qubits[j]];
                temp_m /= d;
            }
            idx temp_c = c_row;
            for (idx j = ctrl_size; j-- > 0;) {
                row_base += (temp_c % d) * d_pow[n - 1 - ctrl[j]];
                temp_c /= d;
            }

            // 2. Inner loop over columns to execute the 2D matrix
            // transformation
            for (idx m_col = 0; m_col < spectator_dim; ++m_col) {
                for (idx c_col = 0; c_col < ctrl_dim; ++c_col) {

                    idx col_base = 0;
                    temp_m = m_col;
                    for (idx j = 0;
                         j < static_cast<idx>(spectator_qubits.size()); ++j) {
                        col_base +=
                            (temp_m % d) * d_pow[n - 1 - spectator_qubits[j]];
                        temp_m /= d;
                    }
                    temp_c = c_col;
                    for (idx j = ctrl_size; j-- > 0;) {
                        col_base += (temp_c % d) * d_pow[n - 1 - ctrl[j]];
                        temp_c /= d;
                    }

                    // 3. Apply the scaling: state(i, j) *= U(i) * conj(U(j))
                    for (idx r_row = 0; r_row < dimA; ++r_row) {
                        const std::complex<double> u_row =
                            A_powers(r_row, c_row);
                        const idx final_row =
                            row_base + target_displacements[r_row];

                        for (idx r_col = 0; r_col < dimA; ++r_col) {
                            const std::complex<double> u_col_conj =
                                std::conj(A_powers(r_col, c_col));
                            const idx final_col =
                                col_base + target_displacements[r_col];

                            state(final_row, final_col) *= (u_row * u_col_conj);
                        }
                    }
                }
            }
        }
    }
}

/**
 * @brief Applies the controlled multi-qudit diagonal gate \a A to the part \a
 * target of the multi-partite density matrix \a state
 *
 * @param state Eigen expression
 * @param A Eigen expression (2^k x 1 vector of diagonal elements)
 * @param target Subsystem indexes where the gate \a A is applied
 * @param ctrl Vector of control qudit indices
 * @param shift Vector of control shifts added to each control qudit value
 * modulo d (e.g., shifts the activating control configuration)
 * @param n Number of qudits
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
} // namespace qpp::internal::kernels::qudit

#endif /* QPP_INTERNAL_KERNELS_QUDIT_APPLY_CTRL_DIAG_HPP_ */
