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
 * @file qpp/internal/kernels/qudit/apply_ctrl.hpp
 * @brief Internal highly optimized critical functions for qpp::applyCTRL()
 */

#ifndef QPP_INTERNAL_KERNELS_QUDIT_APPLY_CTRL_HPP_
#define QPP_INTERNAL_KERNELS_QUDIT_APPLY_CTRL_HPP_

#include <cassert>
#include <vector>
#ifndef NDEBUG
#include <set>
#endif

#include <Eigen/Dense>

#include "qpp/internal/util.hpp"
#include "qpp/types.hpp"

namespace qpp::internal::kernels::qudit {
/**
 * @brief Applies the controlled 1-qudit gate \a A to the qudit \a i of the
 * multi-partite state vector \a state in-place
 *
 * @param state Eigen expression (modified in-place)
 * @param A Eigen expression (2 x 1 vector of diagonal elements)
 * @param i Target subsystem index
 * @param ctrl Vector of control qudit indices
 * @param shift Vector of control values (0: positive/|1> control, 1:
 * negative/|0> control)
 * @param n Number of qudits
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] void apply_ctrl_psi_1q_inplace(
    Eigen::MatrixBase<Derived1>& state, const Eigen::MatrixBase<Derived2>& A,
    const std::vector<idx>& ctrl, idx i, const std::vector<idx>& shift, idx n) {
    using Scalar = typename Derived1::Scalar;
    const idx D = static_cast<idx>(std::size_t{1} << n);

    // Input Validation
    assert(i < n && "Target qudit index i must be less than n");
    assert(static_cast<idx>(state.size()) == D &&
           "State vector size must be 2^n");
    assert(A.rows() == 2 && A.cols() == 2 && "Gate A must be a 2x2 matrix");
    assert(ctrl.size() == shift.size() &&
           "ctrl and shift vectors must have the same size");
    const idx ctrl_size = static_cast<idx>(ctrl.size());

#ifndef NDEBUG
    for (idx c_idx = 0; c_idx < ctrl_size; ++c_idx) {
        const idx c = ctrl[c_idx];
        assert(c < n && "Control qudit index c must be less than n");
        assert(c != i && "Target qudit i cannot also be a control qudit");
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

    const Scalar a00 = A.coeff(0, 0);
    const Scalar a01 = A.coeff(0, 1);
    const Scalar a10 = A.coeff(1, 0);
    const Scalar a11 = A.coeff(1, 1);

#ifdef QPP_OPENMP
#pragma omp parallel for
#endif // QPP_OPENMP
    for (idx L = 0; L < D; L += jump) {
        for (idx R = 0; R < step; ++R) {
            const idx k0 = L + R;
            const idx k1 = k0 + step;

            if (((static_cast<std::size_t>(k0) &
                  static_cast<std::size_t>(expected_pattern_for_ones)) ==
                 static_cast<std::size_t>(expected_pattern_for_ones)) &&
                ((static_cast<std::size_t>(k0) &
                  static_cast<std::size_t>(expected_zero_mask)) == 0)) {

                const Scalar psi_k0 = state.coeff(k0);
                const Scalar psi_k1 = state.coeff(k1);

                state.coeffRef(k0) = (a00 * psi_k0) + (a01 * psi_k1);
                state.coeffRef(k1) = (a10 * psi_k0) + (a11 * psi_k1);
            }
        }
    }
}

/**
 * @brief Applies the controlled 1-qudit gate \a A to the qudit \a i of the
 * multi-partite state vector \a state
 *
 * @param state Eigen expression
 * @param A Eigen expression (2 x 1 vector of diagonal elements)
 * @param i Target subsystem index
 * @param ctrl Vector of control qudit indices
 * @param shift Vector of control values
 * @param n Number of qudits
 * @return Controlled gate \a A applied to the qudit \a i of \a state
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] expr_t<Derived1>
apply_ctrl_psi_1q(const Eigen::MatrixBase<Derived1>& state,
                  const Eigen::MatrixBase<Derived2>& A,
                  const std::vector<idx>& ctrl, idx i,
                  const std::vector<idx>& shift, idx n) {
    // Deep copy for functional return
    expr_t<Derived1> result = state;

    // Apply transformation in-place
    apply_ctrl_psi_1q_inplace(result, A, ctrl, i, shift, n);

    return result;
}

/**
 * @brief Applies the controlled 2-qudit gate \a A to the qudits \a i and \a j
 * of the multi-partite state vector \a state in-place
 *
 * @param state Eigen expression (modified in-place)
 * @param A Eigen expression (4 x 1 vector of diagonal elements)
 * @param ctrl Vector of control qudit indices
 * @param i Target subsystem index
 * @param j Target subsystem index
 * @param shift Vector of control values (0: positive/|1> control, 1:
 * negative/|0> control)
 * @param n Number of qudits
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] void
apply_ctrl_psi_2q_inplace(Eigen::MatrixBase<Derived1>& state,
                          const Eigen::MatrixBase<Derived2>& A,
                          const std::vector<idx>& ctrl, idx i, idx j,
                          const std::vector<idx>& shift, idx n) {
    using Scalar = typename Derived1::Scalar;
    const idx D = static_cast<idx>(std::size_t{1} << n);

    // Input Validation
    assert(i < n && j < n && i != j &&
           "Target qudit indices i and j must be distinct and less than n");
    assert(static_cast<idx>(state.size()) == D &&
           "State vector size must be 2^n");
    assert(A.rows() == 4 && A.cols() == 4 && "Gate A must be a 4x4 matrix");
    assert(ctrl.size() == shift.size() &&
           "ctrl and shift vectors must have the same size");
    const idx ctrl_size = static_cast<idx>(ctrl.size());

#ifndef NDEBUG
    for (idx c_idx = 0; c_idx < ctrl_size; ++c_idx) {
        const idx c = ctrl[c_idx];
        assert(c < n && "Control qudit index c must be less than n");
        assert(c != i && c != j &&
               "Target qudits cannot also be control qudits");
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

    // Extract A entries
    const Scalar a00 = A.coeff(0, 0), a01 = A.coeff(0, 1), a02 = A.coeff(0, 2),
                 a03 = A.coeff(0, 3);
    const Scalar a10 = A.coeff(1, 0), a11 = A.coeff(1, 1), a12 = A.coeff(1, 2),
                 a13 = A.coeff(1, 3);
    const Scalar a20 = A.coeff(2, 0), a21 = A.coeff(2, 1), a22 = A.coeff(2, 2),
                 a23 = A.coeff(2, 3);
    const Scalar a30 = A.coeff(3, 0), a31 = A.coeff(3, 1), a32 = A.coeff(3, 2),
                 a33 = A.coeff(3, 3);

#ifdef QPP_OPENMP
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

        const idx k01 = static_cast<idx>(static_cast<std::size_t>(k00) |
                                         static_cast<std::size_t>(s_j));
        const idx k10 = static_cast<idx>(static_cast<std::size_t>(k00) |
                                         static_cast<std::size_t>(s_i));
        const idx k11 = static_cast<idx>(static_cast<std::size_t>(k00) |
                                         static_cast<std::size_t>(s_i) |
                                         static_cast<std::size_t>(s_j));

        const Scalar psi00 = state.coeff(k00);
        const Scalar psi01 = state.coeff(k01);
        const Scalar psi10 = state.coeff(k10);
        const Scalar psi11 = state.coeff(k11);

        state.coeffRef(k00) =
            (a00 * psi00) + (a01 * psi01) + (a02 * psi10) + (a03 * psi11);
        state.coeffRef(k01) =
            (a10 * psi00) + (a11 * psi01) + (a12 * psi10) + (a13 * psi11);
        state.coeffRef(k10) =
            (a20 * psi00) + (a21 * psi01) + (a22 * psi10) + (a23 * psi11);
        state.coeffRef(k11) =
            (a30 * psi00) + (a31 * psi01) + (a32 * psi10) + (a33 * psi11);
    }
}

/**
 * @brief Applies the controlled 2-qudit gate \a A to the qudits \a i and \a j
 * of the multi-partite state vector \a state
 *
 * @param state Eigen expression
 * @param A Eigen expression (4 x 1 vector of diagonal elements)
 * @param ctrl Vector of control qudit indices
 * @param i Target subsystem index
 * @param j Target subsystem index
 * @param shift Vector of control values
 * @param n Number of qudits
 * @return Controlled gate \a A applied to the qudits \a i and \a j of \a state
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] expr_t<Derived1>
apply_ctrl_psi_2q(const Eigen::MatrixBase<Derived1>& state,
                  const Eigen::MatrixBase<Derived2>& A,
                  const std::vector<idx>& ctrl, idx i, idx j,
                  const std::vector<idx>& shift, idx n) {
    // Functional deep copy
    expr_t<Derived1> result = state;

    // Apply transformation in-place on the copy
    apply_ctrl_psi_2q_inplace(result, A, ctrl, i, j, shift, n);

    return result;
}

/**
 * @brief Applies the multi-qudit controlled gate \a A to the part \a target of
 * the multi-partite state vector \a state in-place
 *
 * @param state Eigen expression (modified in-place)
 * @param A Eigen expression (2^k x 1 vector of diagonal elements)
 * @param ctrl Vector of control qudit indices
 * @param target Subsystem indexes where the gate \a A is applied
 * @param shift Vector of control values (0: positive/|1> control, 1:
 * negative/|0> control)
 * @param n Number of qudits
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] void apply_ctrl_psi_kq_inplace(
    Eigen::MatrixBase<Derived1>& state, const Eigen::MatrixBase<Derived2>& A,
    const std::vector<idx>& ctrl, const std::vector<idx>& target,
    const std::vector<idx>& shift, idx n) {
    using Scalar = typename Derived1::Scalar;
    using EigenVector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

    const idx k = static_cast<idx>(target.size());
    const idx dim = static_cast<idx>(std::size_t{1} << k);
    const idx outer_dim = static_cast<idx>(std::size_t{1} << (n - k));

    // Input Validation
    assert(static_cast<idx>(A.rows()) == dim &&
           static_cast<idx>(A.cols()) == dim &&
           "Gate A must be a 2^k x 2^k matrix");
    assert(ctrl.size() == shift.size() &&
           "ctrl and shift vectors must have the same size");
    const idx ctrl_size = static_cast<idx>(ctrl.size());

#ifndef NDEBUG
    const idx D = static_cast<idx>(std::size_t{1} << n);
    for (idx t : target) {
        assert(t < n && "Target qudit index must be less than n");
    }
    assert(static_cast<idx>(state.size()) == D &&
           "State must be a 2^n x 1 vector");

    if (k > 1) {
        std::vector<idx> sorted_target = target;
        std::sort(sorted_target.begin(), sorted_target.end());
        assert(std::adjacent_find(sorted_target.begin(), sorted_target.end()) ==
                   sorted_target.end() &&
               "Target qudit indices must be distinct");
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

    // Index Mappings
    std::vector<idx> inner_idx(dim);
    for (idx r = 0; r < dim; ++r) {
        idx index_r = 0;
        for (idx j = 0; j < k; ++j) {
            if ((static_cast<std::size_t>(r) >> j) & std::size_t{1}) {
                index_r += static_cast<idx>(std::size_t{1}
                                            << (n - 1 - target[k - 1 - j]));
            }
        }
        inner_idx[r] = index_r;
    }

    std::vector<bool> is_target(n, false);
    for (idx q : target) {
        is_target[q] = true;
    }
    std::vector<idx> spectator_qudits;
    for (idx i = 0; i < n; ++i) {
        if (!is_target[i]) {
            spectator_qudits.push_back(i);
        }
    }

    std::vector<idx> outer_idx(outer_dim);
    for (idx m = 0; m < outer_dim; ++m) {
        idx i_base = 0;
        for (idx j = 0; j < static_cast<idx>(spectator_qudits.size()); ++j) {
            if ((static_cast<std::size_t>(m) >> j) & std::size_t{1}) {
                i_base += static_cast<idx>(std::size_t{1}
                                           << (n - 1 - spectator_qudits[j]));
            }
        }
        outer_idx[m] = i_base;
    }

    // Parallel Execution
#ifdef QPP_OPENMP
#pragma omp parallel for
#endif // QPP_OPENMP
    for (idx m = 0; m < outer_dim; ++m) {
        const idx i_base = outer_idx[m];

        if (((static_cast<std::size_t>(i_base) & expected_pattern_for_ones) ==
             expected_pattern_for_ones) &&
            ((static_cast<std::size_t>(i_base) & expected_zero_mask) == 0)) {

            EigenVector input_block(dim);
            for (idx c = 0; c < dim; ++c) {
                input_block(c) = state(i_base + inner_idx[c]);
            }

            EigenVector new_block = A * input_block;

            for (idx r = 0; r < dim; ++r) {
                state(i_base + inner_idx[r]) = new_block(r);
            }
        }
    }
}

/**
 * @brief Applies the multi-qudit controlled gate \a A to the part \a target of
 * the multi-partite state vector \a state
 *
 * @param state Eigen expression
 * @param A Eigen expression (2^k x 1 vector of diagonal elements)
 * @param ctrl Vector of control qudit indices
 * @param target Subsystem indexes where the gate \a A is applied
 * @param shift Vector of control values
 * @param n Number of qudits
 * @return Controlled gate \a A applied to the part \a target of \a state
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] qpp::expr_t<Derived1>
apply_ctrl_psi_kq(const Eigen::MatrixBase<Derived1>& state,
                  const Eigen::MatrixBase<Derived2>& A,
                  const std::vector<idx>& ctrl, const std::vector<idx>& target,
                  const std::vector<idx>& shift, idx n) {
    // Functional deep copy
    qpp::expr_t<Derived1> result = state;

    // Delegate to in-place implementation
    apply_ctrl_psi_kq_inplace(result, A, ctrl, target, shift, n);

    return result;
}

/**
 * @brief Applies the controlled 1-qudit gate \a A to the qudit \a i of the
 * multi-partite density matrix \a state in-place
 *
 * @param state Eigen expression (modified in-place)
 * @param A Eigen expression (2x2 matrix)
 * @param i Target subsystem index
 * @param ctrl Vector of control qudit indices
 * @param shift Vector of control values (0: positive/|1> control, 1:
 * negative/|0> control)
 * @param n Number of qudits
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] void apply_ctrl_rho_1q_inplace(
    Eigen::MatrixBase<Derived1>& state, const Eigen::MatrixBase<Derived2>& A,
    const std::vector<idx>& ctrl, idx i, const std::vector<idx>& shift, idx n) {
    using Scalar = typename Derived1::Scalar;
    using Matrix2 = Eigen::Matrix2<Scalar>;
    const idx D = static_cast<idx>(std::size_t{1} << n);

    // Input Validation
    assert(i < n && "Target qudit index i must be less than n");
    assert(static_cast<idx>(state.rows()) == D &&
           static_cast<idx>(state.cols()) == D &&
           "State must be a square matrix sized 2^n x 2^n");
    assert(A.rows() == 2 && A.cols() == 2 && "Gate A must be a 2x2 matrix");
    assert(ctrl.size() == shift.size() &&
           "ctrl and shift vectors must have the same size");
    const idx ctrl_size = static_cast<idx>(ctrl.size());

#ifndef NDEBUG
    for (idx c_idx = 0; c_idx < ctrl_size; ++c_idx) {
        const idx c = ctrl[c_idx];
        assert(c < n && "Control qudit index c must be less than n");
        assert(c != i && "Control qudit cannot also be a target qudit");
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
    const idx p_i = n - 1 - i;
    const idx s_i = static_cast<idx>(std::size_t{1} << p_i);
    const idx D_spec = D / 2;
    const idx p_i_plus_1 = p_i + 1;
    const idx low_mask = s_i - 1;

    const Matrix2 A_block = A;
    const Matrix2 A_dagger = A_block.adjoint();
    const idx total_iterations = D_spec * D_spec;

#ifdef QPP_OPENMP
#pragma omp parallel for
#endif // QPP_OPENMP
    for (idx it = 0; it < total_iterations; ++it) {
        const idx s = it / D_spec;
        const idx s_prime = it % D_spec;

        // Row indices
        const idx s_high_r =
            static_cast<idx>(static_cast<std::size_t>(s) >> p_i);
        const idx s_low_r = static_cast<idx>(
            static_cast<std::size_t>(s) & static_cast<std::size_t>(low_mask));
        const idx r0 = (s_high_r << p_i_plus_1) | s_low_r;
        const idx r1 = r0 + s_i;
        const bool row_ctrl = control_is_met(r0);

        // Column indices
        const idx s_prime_high_c =
            static_cast<idx>(static_cast<std::size_t>(s_prime) >> p_i);
        const idx s_prime_low_c =
            static_cast<idx>(static_cast<std::size_t>(s_prime) &
                             static_cast<std::size_t>(low_mask));
        const idx c0 = (s_prime_high_c << p_i_plus_1) | s_prime_low_c;
        const idx c1 = c0 + s_i;
        const bool col_ctrl = control_is_met(c0);

        if (!row_ctrl && !col_ctrl) {
            continue;
        }

        Matrix2 rho_block;
        rho_block << state.coeff(r0, c0), state.coeff(r0, c1),
            state.coeff(r1, c0), state.coeff(r1, c1);

        Matrix2 result_block;
        if (row_ctrl && col_ctrl) {
            result_block.noalias() = A_block * rho_block * A_dagger;
        } else if (row_ctrl) {
            result_block.noalias() = A_block * rho_block;
        } else {
            result_block.noalias() = rho_block * A_dagger;
        }

        state.coeffRef(r0, c0) = result_block(0, 0);
        state.coeffRef(r0, c1) = result_block(0, 1);
        state.coeffRef(r1, c0) = result_block(1, 0);
        state.coeffRef(r1, c1) = result_block(1, 1);
    }
}

/**
 * @brief Applies the controlled 1-qudit gate \a A to the qudit \a i of the
 * multi-partite density matrix \a state
 *
 * @param state Eigen expression
 * @param A Eigen expression (2x2 matrix)
 * @param i Target subsystem index
 * @param ctrl Vector of control qudit indices
 * @param shift Vector of control values
 * @param n Number of qudits
 * @return Controlled gate \a A applied to the qudits \a i of \a state
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] expr_t<Derived1>
apply_ctrl_rho_1q(const Eigen::MatrixBase<Derived1>& state,
                  const Eigen::MatrixBase<Derived2>& A,
                  const std::vector<idx>& ctrl, idx i,
                  const std::vector<idx>& shift, idx n) {
    // Deep copy for functional interface
    expr_t<Derived1> result = state;

    // Apply transformation in-place
    apply_ctrl_rho_1q_inplace(result, A, ctrl, i, shift, n);

    return result;
}

/**
 * @brief Applies the controlled 2-qudit gate \a A to the qudits \a i and \a j
 * of the multi-partite density matrix \a state in-place
 *
 * @param state Eigen expression (modified in-place)
 * @param A Eigen expression (4x4 matrix)
 * @param i Target subsystem index
 * @param j Target subsystem index
 * @param ctrl Vector of control qudit indices
 * @param shift Vector of control values (0: positive/|1> control, 1:
 * negative/|0> control)
 * @param n Number of qudits
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] void
apply_ctrl_rho_2q_inplace(Eigen::MatrixBase<Derived1>& state,
                          const Eigen::MatrixBase<Derived2>& A,
                          const std::vector<idx>& ctrl, idx i, idx j,
                          const std::vector<idx>& shift, idx n) {
    // Type and Dimension Setup
    using Scalar = typename Derived1::Scalar;
    using ComputeBlockType = Eigen::Matrix<Scalar, 4, 4>;

    const idx D [[maybe_unused]] = static_cast<idx>(std::size_t{1} << n);
    const idx ctrl_size = static_cast<idx>(ctrl.size());

    // Input Validation
    assert(i < n && j < n && i != j &&
           "Target qudits must be distinct and < n");
    assert(static_cast<idx>(state.rows()) == D &&
           static_cast<idx>(state.cols()) == D && "State must be 2^n x 2^n");
    assert(A.rows() == 4 && A.cols() == 4 && "Gate A must be a 4x4 matrix");
    assert(ctrl.size() == shift.size() && "ctrl/shift size mismatch");

#ifndef NDEBUG
    std::set<idx> target_set = {i, j};
    for (idx c_idx = 0; c_idx < ctrl_size; ++c_idx) {
        const idx c = ctrl[c_idx];
        assert(c < n && target_set.find(c) == target_set.end() &&
               "Control/Target overlap");
        assert((shift[c_idx] == 0 || shift[c_idx] == 1) &&
               "Shift must be 0 or 1");
    }
#endif

    const idx i_phys = n - i - 1;
    const idx j_phys = n - j - 1;
    const idx P_i = static_cast<idx>(std::size_t{1} << i_phys);
    const idx P_j = static_cast<idx>(std::size_t{1} << j_phys);
    const idx D_rest =
        (n >= 2) ? static_cast<idx>(std::size_t{1} << (n - 2)) : 1;

    const ComputeBlockType U = A.template cast<Scalar>();
    const ComputeBlockType U_adj = U.adjoint();

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

#ifdef QPP_OPENMP
#pragma omp parallel for default(none)                                         \
    shared(D_rest, n, P_i, P_j, U, U_adj, state, i_phys, j_phys,               \
               expected_pattern_for_ones, expected_zero_mask)
#endif
    for (idx r = 0; r < D_rest; ++r) {
        idx r_base_row = 0;
        idx current_r = r;
        for (idx q = 0; q < n; ++q) {
            if (q != i_phys && q != j_phys) {
                if (current_r & 1) {
                    r_base_row |= (static_cast<idx>(1) << q);
                }
                current_r >>= 1;
            }
        }

        const bool row_ctrl =
            ((static_cast<std::size_t>(r_base_row) &
              expected_pattern_for_ones) == expected_pattern_for_ones) &&
            ((static_cast<std::size_t>(r_base_row) & expected_zero_mask) == 0);

        const idx vec_k_row[4] = {r_base_row, r_base_row + P_j,
                                  r_base_row + P_i, r_base_row + P_i + P_j};

        for (idx c = 0; c < D_rest; ++c) {
            idx r_base_col = 0;
            idx current_c = c;
            for (idx q = 0; q < n; ++q) {
                if (q != i_phys && q != j_phys) {
                    if (current_c & 1) {
                        r_base_col |= (static_cast<idx>(1) << q);
                    }
                    current_c >>= 1;
                }
            }

            const bool col_ctrl =
                ((static_cast<std::size_t>(r_base_col) &
                  expected_pattern_for_ones) == expected_pattern_for_ones) &&
                ((static_cast<std::size_t>(r_base_col) & expected_zero_mask) ==
                 0);

            if (!row_ctrl && !col_ctrl) {
                continue;
            }

            const idx vec_k_col[4] = {r_base_col, r_base_col + P_j,
                                      r_base_col + P_i, r_base_col + P_i + P_j};

            ComputeBlockType M;
            for (int row = 0; row < 4; ++row) {
                for (int col = 0; col < 4; ++col) {
                    M(row, col) = state.coeff(vec_k_row[row], vec_k_col[col]);
                }
            }

            ComputeBlockType M_prime;
            if (row_ctrl && col_ctrl) {
                M_prime.noalias() = U * M * U_adj;
            } else if (row_ctrl) {
                M_prime.noalias() = U * M;
            } else {
                M_prime.noalias() = M * U_adj;
            }

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
 * @brief Applies the controlled 2-qudit gate \a A to the qudits \a i and \a j
 * of the multi-partite density matrix \a state
 *
 * @param state Eigen expression
 * @param A Eigen expression (4x4 matrix)
 * @param i Target subsystem index
 * @param j Target subsystem index
 * @param ctrl Vector of control qudit indices
 * @param shift Vector of control values
 * @param n Number of qudits
 * @return Controlled gate \a A applied to the qudits \a i and \a j of \a state
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] expr_t<Derived1>
apply_ctrl_rho_2q(const Eigen::MatrixBase<Derived1>& state,
                  const Eigen::MatrixBase<Derived2>& A,
                  const std::vector<idx>& ctrl, idx i, idx j,
                  const std::vector<idx>& shift, idx n) {
    // Deep copy for the functional return
    expr_t<Derived1> result = state.derived();

    // Delegate to the in-place version
    apply_ctrl_rho_2q_inplace(result, A, ctrl, i, j, shift, n);

    return result;
}

/**
 * @brief Applies the controlled multi-qudit gate \a A to the part \a target of
 * the multi-partite density matrix \a state in-place
 *
 * @param state Eigen expression (modified in-place)
 * @param A Eigen expression (2^k x 1 vector of diagonal elements)
 * @param target Subsystem indexes where the gate \a A is applied
 * @param ctrl Vector of control qudit indices
 * @param shift Vector of control values (0: positive/|1> control, 1:
 * negative/|0> control)
 * @param n Number of qudits
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] void apply_ctrl_rho_kq_inplace(
    Eigen::MatrixBase<Derived1>& state, const Eigen::MatrixBase<Derived2>& A,
    const std::vector<idx>& ctrl, const std::vector<idx>& target,
    const std::vector<idx>& shift, idx n) {
    using Scalar = typename Derived1::Scalar;
    using ComputeMatrixType = Eigen::Matrix<Scalar, -1, -1>;

    const idx k = static_cast<idx>(target.size());
    const idx D_k = (k == 0) ? 1 : static_cast<idx>(std::size_t{1} << k);
    const idx D [[maybe_unused]] = static_cast<idx>(std::size_t{1} << n);

    // Input Validation
    assert(static_cast<idx>(A.rows()) == D_k &&
           static_cast<idx>(A.cols()) == D_k &&
           "Gate A must be a 2^k x 2^k matrix");
    assert(static_cast<idx>(state.rows()) == D &&
           static_cast<idx>(state.cols()) == D && "State must be 2^n x 2^n");
    assert(static_cast<idx>(ctrl.size()) == static_cast<idx>(shift.size()) &&
           "ctrl/shift size mismatch");

#ifndef NDEBUG
    std::vector<idx> check_set = target;
    check_set.insert(check_set.end(), ctrl.begin(), ctrl.end());
    std::sort(check_set.begin(), check_set.end());
    assert(std::unique(check_set.begin(), check_set.end()) == check_set.end() &&
           "Overlap or duplicate indices detected");
#endif

    // Pre-calculations
    std::vector<idx> P_gate_basis(k);
    std::vector<idx> target_phys(k);
    for (idx l = 0; l < k; ++l) {
        target_phys[l] = n - target[l] - 1;
        P_gate_basis[l] = static_cast<idx>(std::size_t{1} << target_phys[l]);
    }

    std::vector<idx> rest_phys = fast_complement(target_phys, n);
    const idx D_rest =
        (n >= k) ? static_cast<idx>(std::size_t{1} << (n - k)) : 1;

    const ComputeMatrixType U = A.template cast<Scalar>();
    const ComputeMatrixType U_adj = U.adjoint();

    // Control Masks
    idx expected_pattern_for_ones = 0;
    idx expected_zero_mask = 0;
    for (idx c_idx = 0; c_idx < static_cast<idx>(ctrl.size()); ++c_idx) {
        const idx bit =
            static_cast<idx>(std::size_t{1} << (n - 1 - ctrl[c_idx]));
        if (shift[c_idx] == 0) {
            expected_pattern_for_ones |= bit;
        } else {
            expected_zero_mask |= bit;
        }
    }

#ifdef QPP_OPENMP
#pragma omp parallel for default(none)                                         \
    shared(D_rest, n, k, rest_phys, P_gate_basis, state, U, U_adj, D_k,        \
               expected_pattern_for_ones, expected_zero_mask)
#endif
    for (idx r = 0; r < D_rest; ++r) {
        idx r_base_row = 0;
        idx current_r = r;
        for (const auto& q_phys : rest_phys) {
            if (current_r & 1) {
                r_base_row |= (static_cast<idx>(1) << q_phys);
            }
            current_r >>= 1;
        }

        const bool row_ctrl =
            ((static_cast<std::size_t>(r_base_row) &
              expected_pattern_for_ones) == expected_pattern_for_ones) &&
            ((static_cast<std::size_t>(r_base_row) & expected_zero_mask) == 0);

        std::vector<idx> vec_k_row(D_k);
        for (idx m = 0; m < D_k; ++m) {
            idx target_component = 0;
            for (idx l = 0; l < k; ++l) {
                if ((m >> (k - 1 - l)) & 1) {
                    target_component += P_gate_basis[l];
                }
            }
            vec_k_row[m] = r_base_row + target_component;
        }

        for (idx c = 0; c < D_rest; ++c) {
            idx r_base_col = 0;
            idx current_c = c;
            for (const auto& q_phys : rest_phys) {
                if (current_c & 1) {
                    r_base_col |= (static_cast<idx>(1) << q_phys);
                }
                current_c >>= 1;
            }

            const bool col_ctrl =
                ((static_cast<std::size_t>(r_base_col) &
                  expected_pattern_for_ones) == expected_pattern_for_ones) &&
                ((static_cast<std::size_t>(r_base_col) & expected_zero_mask) ==
                 0);

            if (!row_ctrl && !col_ctrl) {
                continue;
            }

            ComputeMatrixType M(D_k, D_k);
            std::vector<idx> vec_k_col(D_k);
            for (idx m = 0; m < D_k; ++m) {
                idx target_component = 0;
                for (idx l = 0; l < k; ++l) {
                    if ((m >> (k - 1 - l)) & 1) {
                        target_component += P_gate_basis[l];
                    }
                }
                vec_k_col[m] = r_base_col + target_component;
            }

            for (idx row = 0; row < D_k; ++row) {
                for (idx col = 0; col < D_k; ++col) {
                    M(row, col) = state.coeff(vec_k_row[row], vec_k_col[col]);
                }
            }

            ComputeMatrixType M_prime;
            if (row_ctrl && col_ctrl) {
                M_prime.noalias() = U * M * U_adj;
            } else if (row_ctrl) {
                M_prime.noalias() = U * M;
            } else {
                M_prime.noalias() = M * U_adj;
            }

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
 * @brief Applies the controlled multi-qudit gate \a A to the part \a target of
 * the multi-partite density matrix \a state
 *
 * @param state Eigen expression
 * @param A Eigen expression
 * @param target Subsystem indexes where the gate \a A is applied
 * @param ctrl Vector of control qudit indices
 * @param shift Vector of control values
 * @param n Number of qudits
 * @return Controlled gate \a A applied to the part \a target of \a state
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] expr_t<Derived1>
apply_ctrl_rho_kq(const Eigen::MatrixBase<Derived1>& state,
                  const Eigen::MatrixBase<Derived2>& A,
                  const std::vector<idx>& ctrl, const std::vector<idx>& target,
                  const std::vector<idx>& shift, idx n) {
    // Functional deep copy
    expr_t<Derived1> result = state.derived();

    // Delegate to in-place implementation
    apply_ctrl_rho_kq_inplace(result, A, ctrl, target, shift, n);

    return result;
}
} // namespace qpp::internal::kernels::qudit

#endif /* QPP_INTERNAL_KERNELS_QUDIT_APPLY_CTRL_HPP_ */
