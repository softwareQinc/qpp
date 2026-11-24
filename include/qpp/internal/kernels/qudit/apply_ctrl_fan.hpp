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
#include <vector>

#include <Eigen/Dense>

#include "qpp/types.hpp"

#include "qpp/internal/kernels/qudit/apply.hpp"
#include "qpp/internal/util.hpp"

namespace qpp::internal::kernels::qudit {
/**
 * @brief Applies a controlled-fan gate to a multi-partite state vector
 *
 * @param psi Input state vector
 * @param Ak Precomputed powers of the gate
 * @param ctrl Control subsystem indexes
 * @param target Target subsystem indexes
 * @param dims Dimensions of the multi-partite system
 * @param shift Inverted control shifts
 * @return Transformed state vector
 */
template <typename Scalar>
[[qpp::critical, qpp::parallel]]
dyn_col_vect<Scalar> apply_ctrl_fan_psi_worker(
    const dyn_col_vect<Scalar>& psi, const std::vector<dyn_mat<Scalar>>& Ak,
    const std::vector<idx>& ctrl, const std::vector<idx>& target,
    const std::vector<idx>& dims, const std::vector<idx>& shift) {
    using RealScalar = typename Eigen::NumTraits<Scalar>::Real;

    const idx D = static_cast<idx>(psi.rows());
    const idx d = dims[ctrl[0]];

    dyn_col_vect<Scalar> result = dyn_col_vect<Scalar>::Zero(D);

#ifdef QPP_ENABLE_OPENMP
#pragma omp parallel for
#endif // QPP_ENABLE_OPENMP
    for (idx r = 0; r < d; ++r) {
        std::vector<idx> ctrl_shift = shift;
        std::transform(ctrl_shift.begin(), ctrl_shift.end(), ctrl_shift.begin(),
                       [r, d](idx elem) { return (elem + r) % d; });

        dyn_col_vect<Scalar> chopped_psi =
            internal::project_ket_on_dits(psi, ctrl_shift, ctrl, dims, D);

        // Explicit conversion avoids MSVC C4244
        dyn_col_vect<Scalar> chopped_psi_bar =
            (psi / Scalar{static_cast<RealScalar>(d)}) - chopped_psi;

        if (!chopped_psi.isZero(RealScalar{0})) {
            for (idx elem : target) {
                internal::kernels::qudit::apply_psi_kq_inplace(
                    chopped_psi, Ak[r], {elem}, dims);
            }
        }

#ifdef QPP_ENABLE_OPENMP
#pragma omp critical
#endif // QPP_ENABLE_OPENMP
        {
            result += chopped_psi + chopped_psi_bar;
        }
    }

    return result;
}

/**
 * @brief Applies the single qudit controlled-gate \a A with multiple
 * control qudits listed in \a ctrl to the part \a target of the
 * multi-partite state vector \a state in-place, i.e., CTRL-\a A-\a A-...-\a A
 *
 * @param state Eigen expression (modified in-place)
 * @param A Eigen expression, single qudit quantum gate
 * @param ctrl Control subsystem indexes
 * @param target Target qudit indexes; the gate \a A is applied on every
 * one of them depending on the values of the control qudits
 * @param dims Dimensions of the multi-partite system
 * @param shift Optional, performs the control as if the \a ctrl qudits were
 * \f$X\f$-incremented component-wise by \a shift
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] void apply_ctrl_fan_psi_kq_inplace(
    Eigen::MatrixBase<Derived1>& state, const Eigen::MatrixBase<Derived2>& A,
    const std::vector<idx>& ctrl, const std::vector<idx>& target,
    const std::vector<idx>& dims, const std::vector<idx>& shift) {
    auto& rstate = state.derived();
    const dyn_mat<typename Derived2::Scalar>& rA = A.derived();

    const idx D = static_cast<idx>(rstate.rows());
    if (D == 1) {
        return;
    }

    const idx d = dims[ctrl[0]];

    // construct the table of A^k
    std::vector<dyn_mat<typename Derived1::Scalar>> Ak;
    for (idx k = 0; k < d; ++k) {
        Ak.emplace_back(powm(rA, k));
    }

    rstate = apply_ctrl_fan_psi_worker(
        dyn_col_vect<typename Derived1::Scalar>(rstate), Ak, ctrl, target, dims,
        shift);
}

/**
 * @brief Applies the single qudit controlled-gate \a A with multiple
 * control qudits listed in \a ctrl to the part \a target of the
 * multi-partite density matrix \a state in-place, i.e., CTRL-\a A-\a A-...-\a A
 *
 * @param state Eigen expression (modified in-place)
 * @param A Eigen expression, single qudit quantum gate
 * @param ctrl Control subsystem indexes
 * @param target Target qudit indexes
 * @param dims Dimensions of the multi-partite system
 * @param shift Optional, performs the control as if the \a ctrl qudits were
 * \f$X\f$-incremented component-wise by \a shift
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] void apply_ctrl_fan_rho_kq_inplace(
    Eigen::MatrixBase<Derived1>& state, const Eigen::MatrixBase<Derived2>& A,
    const std::vector<idx>& ctrl, const std::vector<idx>& target,
    const std::vector<idx>& dims, const std::vector<idx>& shift) {
    auto& rstate = state.derived();
    const dyn_mat<typename Derived2::Scalar>& rA = A.derived();

    using Scalar = typename Derived1::Scalar;

    const idx D = static_cast<idx>(rstate.rows());
    if (D == 1) {
        return;
    }

    const idx d = dims[ctrl[0]];

    // construct the table of A^k
    std::vector<dyn_mat<Scalar>> Ak;
    for (idx k = 0; k < d; ++k) {
        Ak.emplace_back(powm(rA, k));
    }

    dyn_mat<Scalar> result = dyn_mat<Scalar>::Zero(D, D);

#ifdef QPP_ENABLE_OPENMP
#pragma omp parallel for
#endif // QPP_ENABLE_OPENMP
    for (idx i = 0; i < D; ++i) {
        dyn_col_vect<Scalar> psi_i = dyn_col_vect<Scalar>::Zero(D);
        psi_i(i) = 1;

        dyn_col_vect<Scalar> phi_i_ket = rstate.row(i).adjoint();

        psi_i = apply_ctrl_fan_psi_worker(psi_i, Ak, ctrl, target, dims, shift);

        dyn_row_vect<Scalar> phi_i_bra =
            apply_ctrl_fan_psi_worker(phi_i_ket, Ak, ctrl, target, dims, shift)
                .adjoint();

#ifdef QPP_ENABLE_OPENMP
#pragma omp critical
#endif // QPP_ENABLE_OPENMP
        {
            result += psi_i * phi_i_bra;
        }
    }

    rstate = std::move(result);
}
} // namespace qpp::internal::kernels::qudit

#endif /* QPP_INTERNAL_KERNELS_QUDIT_APPLY_CTRL_FAN_HPP_ */
