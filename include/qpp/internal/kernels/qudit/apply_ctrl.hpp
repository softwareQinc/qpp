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

#include <Eigen/Dense>

#include "qpp/types.hpp"

#include "qpp/internal/kernels/qudit/apply.hpp"
#include "qpp/internal/util.hpp"

namespace qpp::internal::kernels::qudit {
/**
 * @brief Applies the controlled-gate \a A to the part \a target of the
 * multi-partite state vector \a state in-place
 *
 * @param state Eigen expression (modified in-place)
 * @param A Eigen expression
 * @param ctrl Control subsystem indexes
 * @param target Subsystem indexes where the gate \a A is applied
 * @param dims Dimensions of the multi-partite system
 * @param shift Optional, performs the control as if the \a ctrl qudits were
 * \f$X\f$-incremented component-wise by \a shift
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] void apply_ctrl_psi_kq_inplace(
    Eigen::MatrixBase<Derived1>& state, const Eigen::MatrixBase<Derived2>& A,
    const std::vector<idx>& ctrl, const std::vector<idx>& target,
    const std::vector<idx>& dims, const std::vector<idx>& shift) {
    auto& rstate = state.derived();
    const auto& rA = A.derived();

    const idx D = static_cast<idx>(rstate.rows());
    const idx d = dims[ctrl[0]];

    using Scalar = typename Derived1::Scalar;
    using Ket = dyn_col_vect<Scalar>;
    using Mat = dyn_mat<Scalar>;

    std::vector<Mat> Ak;
    Ak.reserve(d);
    for (idx k = 0; k < d; ++k) {
        Ak.emplace_back(powm(rA, k));
    }

    Ket original = rstate;
    Ket result = original;

    for (idx r = 0; r < d; ++r) {
        std::vector<idx> ctrl_shift = shift;
        std::transform(ctrl_shift.begin(), ctrl_shift.end(), ctrl_shift.begin(),
                       [r, d](idx e) { return (e + r) % d; });

        Ket projected =
            internal::project_ket_on_dits(original, ctrl_shift, ctrl, dims, D);

        if (projected.isZero()) {
            continue;
        }

        Ket transformed = projected;
        internal::kernels::qudit::apply_psi_kq_inplace(transformed, Ak[r],
                                                       target, dims);

        result += transformed - projected;
    }

    rstate = std::move(result);
}

/**
 * @brief Applies the controlled-gate \a A to the part \a target of the
 * multi-partite density matrix \a state in-place
 *
 * @param state Eigen expression (modified in-place)
 * @param A Eigen expression
 * @param ctrl Control subsystem indexes
 * @param target Subsystem indexes where the gate \a A is applied
 * @param dims Dimensions of the multi-partite system
 * @param shift Optional, performs the control as if the \a ctrl qudits were
 * \f$X\f$-incremented component-wise by \a shift
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] void apply_ctrl_rho_kq_inplace(
    Eigen::MatrixBase<Derived1>& state, const Eigen::MatrixBase<Derived2>& A,
    const std::vector<idx>& ctrl, const std::vector<idx>& target,
    const std::vector<idx>& dims, const std::vector<idx>& shift) {
    auto& rstate = state.derived();
    const auto& rA = A.derived();

    const idx D = static_cast<idx>(rstate.rows());
    const idx d = dims[ctrl[0]];

    using Scalar = typename Derived1::Scalar;
    using Ket = dyn_col_vect<Scalar>;
    using Mat = dyn_mat<Scalar>;

    std::vector<Mat> Ak;
    Ak.reserve(d);
    for (idx k = 0; k < d; ++k) {
        Ak.emplace_back(powm(rA, k));
    }

    std::vector<std::vector<idx>> ctrl_shifts(d);
    for (idx r = 0; r < d; ++r) {
        ctrl_shifts[r] = shift;
        std::transform(ctrl_shifts[r].begin(), ctrl_shifts[r].end(),
                       ctrl_shifts[r].begin(),
                       [r, d](idx e) { return (e + r) % d; });
    }

    auto applyCTRL_ket_worker = [&](const Ket& v) -> Ket {
        Ket result = v;

        for (idx r = 0; r < d; ++r) {
            Ket projected =
                internal::project_ket_on_dits(v, ctrl_shifts[r], ctrl, dims, D);

            if (projected.isZero()) {
                continue;
            }

            Ket transformed = projected;
            internal::kernels::qudit::apply_psi_kq_inplace(transformed, Ak[r],
                                                           target, dims);

            result += transformed - projected;
        }

        return result;
    };

    // Apply V from the left
    Mat tmp(D, D);

#ifdef QPP_ENABLE_OPENMP
#pragma omp parallel for
#endif // QPP_ENABLE_OPENMP
    for (idx i = 0; i < D; ++i) {
        Ket column = rstate.col(i);
        tmp.col(i) = applyCTRL_ket_worker(column);
    }

    // Apply V^\dagger from the right
    Mat result(D, D);

#ifdef QPP_ENABLE_OPENMP
#pragma omp parallel for
#endif // QPP_ENABLE_OPENMP
    for (idx i = 0; i < D; ++i) {
        Ket row = tmp.row(i).adjoint();
        result.row(i) = applyCTRL_ket_worker(row).adjoint();
    }

    rstate = std::move(result);
}
} // namespace qpp::internal::kernels::qudit

#endif /* QPP_INTERNAL_KERNELS_QUDIT_APPLY_CTRL_HPP_ */
