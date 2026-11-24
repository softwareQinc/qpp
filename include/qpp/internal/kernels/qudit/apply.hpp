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
 * @file qpp/internal/kernels/qudit/apply.hpp
 * @brief Internal highly optimized critical functions for qpp::apply()
 */

#ifndef QPP_INTERNAL_KERNELS_QUDIT_APPLY_HPP_
#define QPP_INTERNAL_KERNELS_QUDIT_APPLY_HPP_

#include <cassert>
#include <vector>

#include <Eigen/Dense>

#include "qpp/types.hpp"

#include "qpp/internal/util.hpp"

namespace qpp::internal::kernels::qudit {
/**
 * @brief Applies the gate \a A to the part \a target of the multi-partite
 * state vector \a state in-place
 *
 * @param state Eigen expression (modified in-place)
 * @param A Eigen expression
 * @param target Subsystem indexes where the gate \a A is applied
 * @param dims Dimensions of the multi-partite system
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] void apply_psi_kq_inplace(
    Eigen::MatrixBase<Derived1>& state, const Eigen::MatrixBase<Derived2>& A,
    const std::vector<idx>& target, const std::vector<idx>& dims) {
    auto& rstate = state.derived();
    const auto& rA = A.derived();
    const idx n = dims.size();
    const idx gate_size = target.size();

    const idx D = static_cast<idx>(rstate.rows()); // total dimension
    const idx DA = static_cast<idx>(rA.rows()); // dimension of gate subsystem

    idx Cdims[internal::maxn];      // local dimensions total
    idx CdimsA[internal::maxn];     // local dimensions gate
    idx CdimsA_bar[internal::maxn]; // local dimensions complement

    // compute the complementary subsystem of gate w.r.t. dims
    std::vector<idx> gate_bar = fast_complement(target, n);
    // number of subsystems that are complementary to the gate
    idx gate_bar_size = gate_bar.size();

    idx DA_bar = 1; // dimension of the complement
    for (idx i = 0; i < gate_bar_size; ++i) {
        DA_bar *= dims[gate_bar[i]];
    }

    for (idx k = 0; k < n; ++k) {
        Cdims[k] = dims[k];
    }
    for (idx k = 0; k < gate_size; ++k) {
        CdimsA[k] = dims[target[k]];
    }
    for (idx k = 0; k < gate_bar_size; ++k) {
        CdimsA_bar[k] = dims[gate_bar[k]];
    }

    if (D == 1) {
        return;
    }

    // worker, computes the coefficient and the index for the ket case, used in
    // #pragma omp parallel for collapse
    auto coeff_idx_ket =
        [&](idx m_,
            idx r_) noexcept -> std::pair<typename Derived1::Scalar, idx> {
        idx indx = 0;
        typename Derived1::Scalar coeff = 0;

        idx Cmidx[internal::maxn];      // the total multi-index
        idx CmidxA[internal::maxn];     // the gate multi-index
        idx CmidxA_bar[internal::maxn]; // the complement multi-index

        // compute the index

        // set the complement multi-index
        internal::n2multiidx(r_, gate_bar_size, CdimsA_bar, CmidxA_bar);
        for (idx k = 0; k < gate_bar_size; ++k) {
            Cmidx[gate_bar[k]] = CmidxA_bar[k];
        }

        // set the gate multi-index
        internal::n2multiidx(m_, gate_size, CdimsA, CmidxA);
        for (idx k = 0; k < gate_size; ++k) {
            Cmidx[target[k]] = CmidxA[k];
        }

        // we now got the total index
        indx = internal::multiidx2n(Cmidx, n, Cdims);

        // compute the coefficient
        // reads from the original unmodified rstate
        for (idx n_ = 0; n_ < DA; ++n_) {
            internal::n2multiidx(n_, gate_size, CdimsA, CmidxA);
            for (idx k = 0; k < gate_size; ++k) {
                Cmidx[target[k]] = CmidxA[k];
            }
            coeff += rA(m_, n_) * rstate(internal::multiidx2n(Cmidx, n, Cdims));
        }

        return std::make_pair(coeff, indx);
    }; /* end coeff_idx_ket */

    // Must write to a separate result matrix first to avoid aliasing
    // corruption
    auto result = rstate;

#ifdef QPP_ENABLE_OPENMP
#pragma omp parallel for collapse(2)
#endif // QPP_ENABLE_OPENMP
    for (idx m = 0; m < DA; ++m) {
        for (idx r = 0; r < DA_bar; ++r) {
            auto [coeff, index] = coeff_idx_ket(m, r);
            result(index) = coeff;
        }
    }

    // Apply back to the referenced state matrix
    rstate = std::move(result);
}

/**
 * @brief Applies the gate \a A to the part \a target of the multi-partite
 * density matrix \a state in-place
 *
 * @param state Eigen expression (modified in-place)
 * @param A Eigen expression
 * @param target Subsystem indexes where the gate \a A is applied
 * @param dims Dimensions of the multi-partite system
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] void apply_rho_kq_inplace(
    Eigen::MatrixBase<Derived1>& state, const Eigen::MatrixBase<Derived2>& A,
    const std::vector<idx>& target, const std::vector<idx>& dims) {
    auto& rstate = state.derived();
    const auto& rA = A.derived();
    const idx n = dims.size();
    const idx gate_size = target.size();

    const idx D = static_cast<idx>(rstate.rows()); // total dimension
    const idx DA = static_cast<idx>(rA.rows()); // dimension of gate subsystem

    idx Cdims[internal::maxn];      // local dimensions total
    idx CdimsA[internal::maxn];     // local dimensions gate
    idx CdimsA_bar[internal::maxn]; // local dimensions complement

    // compute the complementary subsystem of gate w.r.t. dims
    std::vector<idx> gate_bar = fast_complement(target, n);
    // number of subsystems that are complementary to the gate
    idx gate_bar_size = gate_bar.size();

    idx DA_bar = 1; // dimension of the complement
    for (idx i = 0; i < gate_bar_size; ++i) {
        DA_bar *= dims[gate_bar[i]];
    }

    for (idx k = 0; k < n; ++k) {
        Cdims[k] = dims[k];
    }
    for (idx k = 0; k < gate_size; ++k) {
        CdimsA[k] = dims[target[k]];
    }
    for (idx k = 0; k < gate_bar_size; ++k) {
        CdimsA_bar[k] = dims[gate_bar[k]];
    }

    if (D == 1) {
        return;
    }

    // worker, computes the coefficient and the index for the density matrix
    // case, used in #pragma omp parallel for collapse
    auto coeff_idx_rho = [&](idx m1_, idx r1_, idx m2_, idx r2_) noexcept
        -> std::tuple<typename Derived1::Scalar, idx, idx> {
        idx idxrow = 0;
        idx idxcol = 0;
        typename Derived1::Scalar coeff = 0, lhs = 1, rhs = 1;

        idx Cmidxrow[internal::maxn];      // the total row multi-index
        idx Cmidxcol[internal::maxn];      // the total col multi-index
        idx CmidxArow[internal::maxn];     // the gate row multi-index
        idx CmidxAcol[internal::maxn];     // the gate col multi-index
        idx CmidxA_barrow[internal::maxn]; // the complement row multi-index
        idx CmidxA_barcol[internal::maxn]; // the complement col multi-index

        // compute the ket/bra indexes

        // set the complement multi-index
        internal::n2multiidx(r1_, gate_bar_size, CdimsA_bar, CmidxA_barrow);
        internal::n2multiidx(r2_, gate_bar_size, CdimsA_bar, CmidxA_barcol);
        for (idx k = 0; k < gate_bar_size; ++k) {
            Cmidxrow[gate_bar[k]] = CmidxA_barrow[k];
            Cmidxcol[gate_bar[k]] = CmidxA_barcol[k];
        }

        // set the gate multi-index
        internal::n2multiidx(m1_, gate_size, CdimsA, CmidxArow);
        internal::n2multiidx(m2_, gate_size, CdimsA, CmidxAcol);
        for (idx k = 0; k < gate_size; ++k) {
            Cmidxrow[target[k]] = CmidxArow[k];
            Cmidxcol[target[k]] = CmidxAcol[k];
        }

        // we now got the total row/col indexes
        idxrow = internal::multiidx2n(Cmidxrow, n, Cdims);
        idxcol = internal::multiidx2n(Cmidxcol, n, Cdims);

        // compute the coefficient
        // reads from the original unmodified rstate
        for (idx n1_ = 0; n1_ < DA; ++n1_) {
            internal::n2multiidx(n1_, gate_size, CdimsA, CmidxArow);
            for (idx k = 0; k < gate_size; ++k) {
                Cmidxrow[target[k]] = CmidxArow[k];
            }
            idx idxrowtmp = internal::multiidx2n(Cmidxrow, n, Cdims);

            lhs = rA(m1_, n1_);

            for (idx n2_ = 0; n2_ < DA; ++n2_) {
                internal::n2multiidx(n2_, gate_size, CdimsA, CmidxAcol);
                for (idx k = 0; k < gate_size; ++k) {
                    Cmidxcol[target[k]] = CmidxAcol[k];
                }

                rhs = adjoint(rA)(n2_, m2_);

                idx idxcoltmp = internal::multiidx2n(Cmidxcol, n, Cdims);

                coeff += lhs * rstate(idxrowtmp, idxcoltmp) * rhs;
            }
        }

        return std::make_tuple(coeff, idxrow, idxcol);
    }; /* end coeff_idx_rho */

    // Must write to a separate result matrix first to avoid aliasing
    // corruption
    auto result = rstate;

#ifdef QPP_ENABLE_OPENMP
#pragma omp parallel for collapse(4)
#endif // QPP_ENABLE_OPENMP
    for (idx m1 = 0; m1 < DA; ++m1) {
        for (idx r1 = 0; r1 < DA_bar; ++r1) {
            for (idx m2 = 0; m2 < DA; ++m2) {
                for (idx r2 = 0; r2 < DA_bar; ++r2) {
                    auto coeff_idxes = coeff_idx_rho(m1, r1, m2, r2);
                    result(std::get<1>(coeff_idxes), std::get<2>(coeff_idxes)) =
                        std::get<0>(coeff_idxes);
                }
            }
        }
    }

    // Apply back to the referenced state matrix
    rstate = std::move(result);
}
} // namespace qpp::internal::kernels::qudit

#endif /* QPP_INTERNAL_KERNELS_QUDIT_APPLY_HPP_ */
