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
 * @file qpp/internal/kernels/qudit/apply_ctrl_fan.hpp
 * @brief Internal highly optimized critical functions for qpp::applyCTRL_fan()
 */

#ifndef QPP_INTERNAL_KERNELS_QUDIT_APPLY_CTRL_FAN_HPP_
#define QPP_INTERNAL_KERNELS_QUDIT_APPLY_CTRL_FAN_HPP_

#include <cassert>
#include <optional>
#include <vector>
#include "qpp/internal/util.hpp"
#ifndef NDEBUG
#include <set>
#endif

#include <Eigen/Dense>

#include "qpp/types.hpp"

namespace qpp::internal::kernels::qudit {
/**
 * @brief Applies the single qudit controlled-gate \a A with multiple
 * control qudits listed in \a ctrl to the part \a target of the multi-partite
 * state vector \a state in-place, i.e., CTRL-A-A-...-A
 *
 * @param state Eigen expression (modified in-place)
 * @param A Eigen expression (2x2 matrix)
 * @param ctrl Vector of control qudit indexes
 * @param target Vector of target subsystem indexes
 * @param shift Vector of control shifts added to each control qudit value
 * modulo d (e.g., shifts the activating control configuration)
 * @param d Subsystem dimensions
 * @param n Number of qudits
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] void apply_ctrl_fan_psi_inplace(
    Eigen::MatrixBase<Derived1>& state, const Eigen::MatrixBase<Derived2>& A,
    const std::vector<idx>& ctrl, const std::vector<idx>& target,
    const std::vector<idx>& shift, idx d, idx n) {

    expr_t<Derived1>& rstate = state.derived();
    const idx D = static_cast<idx>(rstate.rows());

    // Early exit condition
    if (d <= 1 || D <= 1 || target.empty()) {
        return;
    }

    // Assertions
    assert(static_cast<idx>(A.rows()) == d && static_cast<idx>(A.cols()) == 1 &&
           "Gate A must be a d x 1 vector of diagonal elements");
    assert(static_cast<idx>(rstate.cols()) == 1 && "State must be d^n x 1");
    assert(shift.size() == ctrl.size() && "Shift size must match ctrl size");
    assert(target.size() == 1 && "Target must contain exactly 1 qudit index");
    assert(target[0] < n && "Target qudit index must be less than n");

    const dyn_mat<typename Derived2::Scalar>& rA = A.derived();
    std::vector<idx> dims(n, d);

    // Precompute the table of A^k matrix powers
    std::vector<dyn_mat<typename Derived1::Scalar>> Ak;
    Ak.reserve(d);
    for (idx k = 0; k < d; ++k) {
        Ak.emplace_back(powm(rA, k));
    }

    // Seed the result vector with the original state to optimize away the
    // baseline (psi / d) mathematical overhead from the loop.
    dyn_col_vect<typename Derived1::Scalar> result = rstate;

#ifdef QPP_OPENMP
#pragma omp parallel for
#endif
    for (idx r = 0; r < d; ++r) {
        std::vector<idx> local_shift = shift;
        std::transform(local_shift.begin(), local_shift.end(),
                       local_shift.begin(),
                       [r, d](idx elem) { return (elem + r) % d; });

        // Project the original state vector onto the control configurations
        dyn_col_vect<typename Derived1::Scalar> chopped_psi =
            internal::project_ket_on_dits(rstate, local_shift, ctrl, dims, D);

        if (chopped_psi.norm() > 0) {
            dyn_col_vect<typename Derived1::Scalar> transformed_psi =
                chopped_psi;
            for (idx elem : target) {
                transformed_psi = apply(transformed_psi, Ak[r], {elem}, dims);
            }

            // Difference vector: modified sub-state minus original sub-state
            dyn_col_vect<typename Derived1::Scalar> delta =
                transformed_psi - chopped_psi;

#ifdef QPP_OPENMP
#pragma omp critical
#endif
            {
                result += delta;
            }
        }
    }

    rstate = std::move(result);
}

/**
 * @brief Applies the single qudit controlled-gate \a A with multiple
 * control qudits listed in \a ctrl to the part \a target of the multi-partite
 * state vector \a state
 *
 * @return Controlled gate \a A applied to the targets of \a state
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] expr_t<Derived1>
apply_ctrl_fan_psi(const Eigen::MatrixBase<Derived1>& state,
                   const Eigen::MatrixBase<Derived2>& A,
                   const std::vector<idx>& ctrl, const std::vector<idx>& target,
                   const std::vector<idx>& shift, idx n) {
    // Functional deep copy
    expr_t<Derived1> result = state.derived();

    // Delegate to in-place implementation
    apply_ctrl_fan_psi_inplace(result, A, ctrl, target, shift, n);

    return result;
}

/**
 * @brief Applies the single qudit controlled-gate \a A with multiple
 * control qudits listed in \a ctrl to the part \a target of the multi-partite
 * density matrix \a state in-place, i.e., CTRL-A-A-...-A
 *
 * @param state Eigen expression (modified in-place)
 * @param A Eigen expression (2x2 matrix)
 * @param ctrl Vector of control qudit indexes
 * @param target Vector of target subsystem indexes
 * @param shift Vector of control shifts added to each control qudit value
 * modulo d (e.g., shifts the activating control configuration)
 * @param d Subsystem dimensions
 * @param n Number of qudits
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] void apply_ctrl_fan_rho_inplace(
    Eigen::MatrixBase<Derived1>& state, const Eigen::MatrixBase<Derived2>& A,
    const std::vector<idx>& ctrl, const std::vector<idx>& target,
    const std::vector<idx>& shift, idx n) {
    const dyn_mat<typename Derived2::Scalar>& rA = A.derived();
    {
        // Density matrix: compute into temporary then move
        dyn_mat<typename Derived1::Scalar> result =
            dyn_mat<typename Derived1::Scalar>::Zero(D, D);

#ifdef QPP_OPENMP
#pragma omp parallel for
#endif
        for (idx i = 0; i < D; ++i) {
            dyn_col_vect<typename Derived1::Scalar> psi_i =
                dyn_col_vect<typename Derived1::Scalar>::Zero(D);
            psi_i(i) = 1;

            dyn_col_vect<typename Derived1::Scalar> phi_i_ket =
                rstate.row(i).adjoint();

            psi_i = applyCTRL_fan_ket(psi_i);
            dyn_row_vect<typename Derived1::Scalar> phi_i_bra =
                applyCTRL_fan_ket(phi_i_ket).adjoint();

#ifdef QPP_OPENMP
#pragma omp critical
#endif
            {
                result += psi_i * phi_i_bra;
            }
        }
        rstate = std::move(result);
    }
}

/**
 * @brief Applies the single qudit controlled-gate \a A with multiple
 * control qudits listed in \a ctrl to the part \a target of the multi-partite
 * density matrix \a state
 *
 * @return Controlled gate \a A applied to the targets of \a state
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] expr_t<Derived1>
apply_ctrl_fan_rho(const Eigen::MatrixBase<Derived1>& state,
                   const Eigen::MatrixBase<Derived2>& A,
                   const std::vector<idx>& ctrl, const std::vector<idx>& target,
                   const std::vector<idx>& shift, idx n) {
    // Functional deep copy
    expr_t<Derived1> result = state.derived();

    // Delegate to in-place implementation
    apply_ctrl_fan_rho_inplace(result, A, ctrl, target, shift, n);

    return result;
}
} // namespace qpp::internal::kernels::qudit

#endif /* QPP_INTERNAL_KERNELS_QUDIT_APPLY_CTRL_FAN_HPP_ */
