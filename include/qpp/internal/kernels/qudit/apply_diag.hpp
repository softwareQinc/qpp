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
 * @file qpp/internal/kernels/qudit/apply_diag.hpp
 * @brief Internal highly optimized critical functions for qpp::apply_diag()
 */

#ifndef QPP_INTERNAL_KERNELS_QUDIT_APPLY_DIAG_HPP_
#define QPP_INTERNAL_KERNELS_QUDIT_APPLY_DIAG_HPP_

#include <cassert>
#include <vector>

#include <Eigen/Dense>

#include "qpp/types.hpp"

#include "qpp/internal/util.hpp"

namespace qpp::internal::kernels::qudit {
/**
 * @brief Applies the diagonal gate \a A to the part \a target of the
 * multi-partite state vector \a state in-place
 *
 * @param state Eigen expression (modified in-place)
 * @param A Eigen expression (D x 1 vector of diagonal elements)
 * @param target Subsystem indexes where the gate A is applied
 * @param dims Dimensions of the multi-partite system
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] void apply_psi_kq_diag_inplace(
    Eigen::MatrixBase<Derived1>& state, const Eigen::MatrixBase<Derived2>& A,
    const std::vector<idx>& target, const std::vector<idx>& dims) {
    auto& rstate = state.derived();
    const auto& rA = A.derived();

    const idx n = dims.size();
    const idx D = static_cast<idx>(rstate.rows());
    const idx DA = static_cast<idx>(rA.size());
    const idx gate_size = target.size();

    if (D == 1) {
        rstate(0) *= rA(0);
        return;
    }

    idx Cdims[internal::maxn];
    idx CdimsA[internal::maxn];
    idx CdimsA_bar[internal::maxn];

    const std::vector<idx> gate_bar = fast_complement(target, n);
    const idx gate_bar_size = gate_bar.size();
    const idx DA_bar = D / DA;

    for (idx k = 0; k < n; ++k) {
        Cdims[k] = dims[k];
    }
    for (idx k = 0; k < gate_size; ++k) {
        CdimsA[k] = dims[target[k]];
    }
    for (idx k = 0; k < gate_bar_size; ++k) {
        CdimsA_bar[k] = dims[gate_bar[k]];
    }

#ifdef QPP_ENABLE_OPENMP
#pragma omp parallel for collapse(2)
#endif // QPP_ENABLE_OPENMP
    for (idx m = 0; m < DA; ++m) {
        for (idx r = 0; r < DA_bar; ++r) {
            idx Cmidx[internal::maxn];
            idx CmidxA[internal::maxn];
            idx CmidxA_bar[internal::maxn];

            // complement multi-index
            internal::n2multiidx(r, gate_bar_size, CdimsA_bar, CmidxA_bar);
            for (idx k = 0; k < gate_bar_size; ++k) {
                Cmidx[gate_bar[k]] = CmidxA_bar[k];
            }

            // gate multi-index
            internal::n2multiidx(m, gate_size, CdimsA, CmidxA);
            for (idx k = 0; k < gate_size; ++k) {
                Cmidx[target[k]] = CmidxA[k];
            }

            // apply diagonal element
            rstate(internal::multiidx2n(Cmidx, n, Cdims)) *= rA(m);
        }
    }
}

/**
 * @brief Applies the diagonal gate \a A to the part \a target of the
 * multi-partite density matrix \a state in-place
 *
 * @param state Eigen expression (modified in-place)
 * @param A Eigen expression (D x 1 vector of diagonal elements)
 * @param target Subsystem indexes where the gate A is applied
 * @param dims Dimensions of the multi-partite system
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] void apply_rho_kq_diag_inplace(
    Eigen::MatrixBase<Derived1>& state, const Eigen::MatrixBase<Derived2>& A,
    const std::vector<idx>& target, const std::vector<idx>& dims) {
    auto& rstate = state.derived();
    const auto& rA = A.derived();

    const idx n = dims.size();
    const idx D = static_cast<idx>(rstate.rows());
    const idx DA = static_cast<idx>(rA.size());
    const idx gate_size = target.size();

    if (D == 1) {
        rstate(0, 0) *= rA(0) * std::conj(rA(0));
        return;
    }

    idx Cdims[internal::maxn];
    idx CdimsA[internal::maxn];
    idx CdimsA_bar[internal::maxn];

    const std::vector<idx> gate_bar = fast_complement(target, n);
    const idx gate_bar_size = gate_bar.size();
    const idx DA_bar = D / DA;

    for (idx k = 0; k < n; ++k) {
        Cdims[k] = dims[k];
    }
    for (idx k = 0; k < gate_size; ++k) {
        CdimsA[k] = dims[target[k]];
    }
    for (idx k = 0; k < gate_bar_size; ++k) {
        CdimsA_bar[k] = dims[gate_bar[k]];
    }

#ifdef QPP_ENABLE_OPENMP
#pragma omp parallel for collapse(4)
#endif // QPP_ENABLE_OPENMP
    for (idx m1 = 0; m1 < DA; ++m1) {
        for (idx r1 = 0; r1 < DA_bar; ++r1) {
            for (idx m2 = 0; m2 < DA; ++m2) {
                for (idx r2 = 0; r2 < DA_bar; ++r2) {
                    idx Cmidxrow[internal::maxn];
                    idx Cmidxcol[internal::maxn];
                    idx CmidxArow[internal::maxn];
                    idx CmidxAcol[internal::maxn];
                    idx CmidxA_barrow[internal::maxn];
                    idx CmidxA_barcol[internal::maxn];

                    // complement multi-indices
                    internal::n2multiidx(r1, gate_bar_size, CdimsA_bar,
                                         CmidxA_barrow);
                    internal::n2multiidx(r2, gate_bar_size, CdimsA_bar,
                                         CmidxA_barcol);
                    for (idx k = 0; k < gate_bar_size; ++k) {
                        Cmidxrow[gate_bar[k]] = CmidxA_barrow[k];
                        Cmidxcol[gate_bar[k]] = CmidxA_barcol[k];
                    }

                    // gate multi-indices
                    internal::n2multiidx(m1, gate_size, CdimsA, CmidxArow);
                    internal::n2multiidx(m2, gate_size, CdimsA, CmidxAcol);
                    for (idx k = 0; k < gate_size; ++k) {
                        Cmidxrow[target[k]] = CmidxArow[k];
                        Cmidxcol[target[k]] = CmidxAcol[k];
                    }

                    const idx row = internal::multiidx2n(Cmidxrow, n, Cdims);
                    const idx col = internal::multiidx2n(Cmidxcol, n, Cdims);

                    // apply A on left and A^\dagger on right
                    rstate(row, col) *= rA(m1) * std::conj(rA(m2));
                }
            }
        }
    }
}
} // namespace qpp::internal::kernels::qudit

#endif /* QPP_INTERNAL_KERNELS_QUDIT_APPLY_DIAG_HPP_ */
