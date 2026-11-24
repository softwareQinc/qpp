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
 * @file qpp/internal/kernels/qudit/syspermute.hpp
 * @brief Internal highly optimized critical functions for qpp::syspermute()
 */

#ifndef QPP_INTERNAL_KERNELS_QUDIT_SYSPERMUTE_HPP_
#define QPP_INTERNAL_KERNELS_QUDIT_SYSPERMUTE_HPP_

#include <cassert>
#include <vector>

#include <Eigen/Dense>

#include "qpp/options.hpp"
#include "qpp/types.hpp"

#include "qpp/internal/util.hpp"

namespace qpp::internal::kernels::qudit {
/**
 * @brief Qudit state vector subsystem permutation
 *
 * Permutes the subsystems of a qudit state vector. The qudit \a perm[\a i] is
 * permuted to the location \a i.
 *
 * @param A Eigen expression
 * @param perm Permutation
 * @param n Number of qudits
 * @return Permuted qudit system, as a dynamic matrix over the same scalar field
 * as \a A
 */
template <typename Derived>
[[qpp::critical, qpp::parallel]] dyn_mat<typename Derived::Scalar>
syspermute_psi_kq(const Eigen::MatrixBase<Derived>& A,
                  const std::vector<idx>& perm, const std::vector<idx>& dims) {
    using Scalar = typename Derived::Scalar;
    const auto& rA = A.derived();

    const idx n = dims.size();
    const idx D = static_cast<idx>(rA.rows());

    idx Cdims[internal::maxn];
    idx Cperm[internal::maxn];

    std::copy(dims.begin(), dims.end(), Cdims);
    std::copy(perm.begin(), perm.end(), Cperm);

    dyn_mat<Scalar> result(D, 1);

#ifdef QPP_ENABLE_OPENMP
#pragma omp parallel for
#endif // QPP_ENABLE_OPENMP
    for (idx i = 0; i < D; ++i) {
        idx midx[internal::maxn];
        idx midxtmp[internal::maxn];
        idx permdims[internal::maxn];

        internal::n2multiidx(i, n, Cdims, midx);

        for (idx k = 0; k < n; ++k) {
            permdims[k] = Cdims[Cperm[k]];
            midxtmp[k] = midx[Cperm[k]];
        }

        const idx pos = internal::multiidx2n(midxtmp, n, permdims);

        result(pos) = rA.data()[i];
    }

    return result;
}

/**
 * @brief Qudit density matrix subsystem permutation
 *
 * Permutes the subsystems of a qudit density matrix. The qudit \a perm[\a i]
 * is permuted to the location \a i.
 *
 * @param A Eigen expression
 * @param perm Permutation
 * @param n Number of qudits
 * @return Permuted qudit system, as a dynamic matrix over the same scalar field
 * as \a A
 */
template <typename Derived>
[[qpp::critical, qpp::parallel]] dyn_mat<typename Derived::Scalar>
syspermute_rho_kq(const Eigen::MatrixBase<Derived>& A,
                  const std::vector<idx>& perm, const std::vector<idx>& dims) {
    using Scalar = typename Derived::Scalar;
    const auto& rA = A.derived();

    const idx n = dims.size();
    const idx D = static_cast<idx>(rA.rows());
    const idx total_systems = 2 * n;
    const idx total_elements = D * D;

    idx Cdims[2 * internal::maxn];
    idx Cperm[2 * internal::maxn];

    for (idx i = 0; i < n; ++i) {
        Cdims[i] = dims[i];
        Cperm[i] = perm[i];

        Cdims[i + n] = dims[i];
        Cperm[i + n] = perm[i] + n;
    }

    dyn_mat<Scalar> result(total_elements, 1);

#ifdef QPP_ENABLE_OPENMP
#pragma omp parallel for
#endif // QPP_ENABLE_OPENMP
    for (idx i = 0; i < total_elements; ++i) {
        idx midx[2 * internal::maxn];
        idx midxtmp[2 * internal::maxn];
        idx permdims[2 * internal::maxn];

        internal::n2multiidx(i, total_systems, Cdims, midx);

        for (idx k = 0; k < total_systems; ++k) {
            permdims[k] = Cdims[Cperm[k]];
            midxtmp[k] = midx[Cperm[k]];
        }

        const idx pos = internal::multiidx2n(midxtmp, total_systems, permdims);

        result(pos) = rA.data()[i];
    }

    result.resize(D, D);

    return result;
}
} // namespace qpp::internal::kernels::qudit

#endif /* QPP_INTERNAL_KERNELS_QUDIT_SYSPERMUTE_HPP_ */
