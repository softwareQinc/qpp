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
 * @file qpp/internal/kernels/qudit/ptranspose.hpp
 * @brief Internal highly optimized critical functions for qpp::ptranspose()
 */

#ifndef QPP_INTERNAL_KERNELS_QUDIT_PTRANSPOSE_HPP_
#define QPP_INTERNAL_KERNELS_QUDIT_PTRANSPOSE_HPP_

#include <cassert>
#include <vector>

#include <Eigen/Dense>

#include "qpp/options.hpp"
#include "qpp/types.hpp"

#include "qpp/internal/util.hpp"

namespace qpp::internal::kernels::qudit {
/**
 * @brief Partial transpose
 *
 * Partial transpose of the multi-partite state vector \a A
 * over the list \a target of subsystems
 *
 * @param A Eigen expression
 * @param target Subsystem indexes
 * @param dims Dimensions of the multi-partite system
 * @return Partial transpose \f$(\cdot)^{T_{subsys}}\f$ over the subsystems
 * \a target in a multi-partite system, as a dynamic matrix over the same
 * scalar field as \a A
 */
template <typename Derived>
[[qpp::critical, qpp::parallel]] dyn_mat<typename Derived::Scalar>
ptranspose_psi_kq(const Eigen::MatrixBase<Derived>& A,
                  const std::vector<idx>& target,
                  const std::vector<idx>& dims) {
    using Scalar = typename Derived::Scalar;
    const auto& rA = A.derived();

    if (target.empty()) {
        return rA * rA.adjoint();
    }

    if (target.size() == dims.size()) {
        return (rA * rA.adjoint()).transpose();
    }

    const idx n = dims.size();
    const idx D = static_cast<idx>(rA.rows());
    const idx n_subsys = target.size();

    idx Cdims[internal::maxn];
    idx Csubsys[internal::maxn];

    std::copy(dims.begin(), dims.end(), Cdims);
    std::copy(target.begin(), target.end(), Csubsys);

    dyn_mat<Scalar> result(D, D);

    for (idx j = 0; j < D; ++j) {
        idx Cmidxcol[internal::maxn];

        internal::n2multiidx(j, n, Cdims, Cmidxcol);

#ifdef QPP_ENABLE_OPENMP
#pragma omp parallel for
#endif // QPP_ENABLE_OPENMP
        for (idx i = 0; i < D; ++i) {
            idx midxcoltmp[internal::maxn];
            idx midxrow[internal::maxn];

            std::copy_n(Cmidxcol, n, midxcoltmp);

            internal::n2multiidx(i, n, Cdims, midxrow);

            for (idx k = 0; k < n_subsys; ++k) {
                std::swap(midxcoltmp[Csubsys[k]], midxrow[Csubsys[k]]);
            }

            result(i, j) = rA(internal::multiidx2n(midxrow, n, Cdims)) *
                           Eigen::numext::conj(
                               rA(internal::multiidx2n(midxcoltmp, n, Cdims)));
        }
    }

    return result;
}

/**
 * @brief Partial transpose
 *
 * Partial transpose of the multi-partite density matrix \a A
 * over the list \a target of subsystems
 *
 * @param A Eigen expression
 * @param target Subsystem indexes
 * @param dims Dimensions of the multi-partite system
 * @return Partial transpose \f$(\cdot)^{T_{subsys}}\f$ over the subsystems
 * \a target in a multi-partite system, as a dynamic matrix over the same
 * scalar field as \a A
 */
template <typename Derived>
[[qpp::critical, qpp::parallel]] dyn_mat<typename Derived::Scalar>
ptranspose_rho_kq(const Eigen::MatrixBase<Derived>& A,
                  const std::vector<idx>& target,
                  const std::vector<idx>& dims) {
    using Scalar = typename Derived::Scalar;
    const auto& rA = A.derived();

    if (target.empty()) {
        return rA;
    }

    if (target.size() == dims.size()) {
        return rA.transpose();
    }

    const idx n = dims.size();
    const idx D = static_cast<idx>(rA.rows());
    const idx n_subsys = target.size();

    idx Cdims[internal::maxn];
    idx Csubsys[internal::maxn];

    std::copy(dims.begin(), dims.end(), Cdims);
    std::copy(target.begin(), target.end(), Csubsys);

    dyn_mat<Scalar> result(D, D);

    for (idx j = 0; j < D; ++j) {
        idx Cmidxcol[internal::maxn];

        internal::n2multiidx(j, n, Cdims, Cmidxcol);

#ifdef QPP_ENABLE_OPENMP
#pragma omp parallel for
#endif // QPP_ENABLE_OPENMP
        for (idx i = 0; i < D; ++i) {
            idx midxcoltmp[internal::maxn];
            idx midxrow[internal::maxn];

            std::copy_n(Cmidxcol, n, midxcoltmp);

            internal::n2multiidx(i, n, Cdims, midxrow);

            for (idx k = 0; k < n_subsys; ++k) {
                std::swap(midxcoltmp[Csubsys[k]], midxrow[Csubsys[k]]);
            }

            const idx row = internal::multiidx2n(midxrow, n, Cdims);
            const idx col = internal::multiidx2n(midxcoltmp, n, Cdims);

            result(i, j) = rA(row, col);
        }
    }

    return result;
}

} // namespace qpp::internal::kernels::qudit

#endif /* QPP_INTERNAL_KERNELS_QUDIT_PTRANSPOSE_HPP_ */
