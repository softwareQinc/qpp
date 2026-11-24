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
 * @file qpp/internal/kernels/qudit/ptrace.hpp
 * @brief Internal highly optimized critical functions for qpp::ptrace()
 */

#ifndef QPP_INTERNAL_KERNELS_QUDIT_PTRACE_HPP_
#define QPP_INTERNAL_KERNELS_QUDIT_PTRACE_HPP_

#include <cassert>
#include <vector>

#include <Eigen/Dense>

#include "qpp/types.hpp"

#include "qpp/internal/util.hpp"

namespace qpp::internal::kernels::qudit {
/**
 * @brief Partial trace
 *
 * Partial trace of the multi-partite state vector \a A over the list \a target
 * of subsystems
 *
 * @param A Eigen expression
 * @param target Subsystem indexes
 * @param dims Dimensions of the multi-partite system
 * @return Partial trace \f$Tr_{subsys}(\cdot)\f$ over the subsystems \a
 * target in a multi-partite system, as a dynamic matrix over the same
 * scalar field as \a A
 */
template <typename Derived>
[[qpp::critical, qpp::parallel]] dyn_mat<typename Derived::Scalar>
ptrace_psi_kq(const Eigen::MatrixBase<Derived>& A,
              const std::vector<idx>& target, const std::vector<idx>& dims) {
    using Scalar = typename Derived::Scalar;
    const auto& rA = A.derived();

    const idx n = dims.size();
    const idx n_subsys = target.size();
    const idx n_subsys_bar = n - n_subsys;

    if (target.size() == dims.size()) {
        dyn_mat<Scalar> result(1, 1);
        result(0, 0) = (rA.adjoint() * rA).value();
        return result;
    }

    if (target.empty()) {
        return rA * rA.adjoint();
    }

    const idx D = static_cast<idx>(rA.rows());

    idx Dsubsys = 1;
    for (idx i = 0; i < n_subsys; ++i) {
        Dsubsys *= dims[target[i]];
    }
    const idx Dsubsys_bar = D / Dsubsys;

    idx Cdims[internal::maxn];
    idx Csubsys[internal::maxn];
    idx Cdimssubsys[internal::maxn];
    idx Csubsys_bar[internal::maxn];
    idx Cdimssubsys_bar[internal::maxn];

    const std::vector<idx> subsys_bar = fast_complement(target, n);

    std::copy(dims.begin(), dims.end(), Cdims);
    std::copy(target.begin(), target.end(), Csubsys);
    std::copy(subsys_bar.begin(), subsys_bar.end(), Csubsys_bar);

    for (idx i = 0; i < n_subsys; ++i) {
        Cdimssubsys[i] = dims[target[i]];
    }
    for (idx i = 0; i < n_subsys_bar; ++i) {
        Cdimssubsys_bar[i] = dims[subsys_bar[i]];
    }

    dyn_mat<Scalar> result(Dsubsys_bar, Dsubsys_bar);

    for (idx j = 0; j < Dsubsys_bar; ++j) {
        idx Cmidxcolsubsys_bar[internal::maxn];

        internal::n2multiidx(j, n_subsys_bar, Cdimssubsys_bar,
                             Cmidxcolsubsys_bar);

#ifdef QPP_ENABLE_OPENMP
#pragma omp parallel for
#endif // QPP_ENABLE_OPENMP
        for (idx i = 0; i < Dsubsys_bar; ++i) {
            idx Cmidxrow[internal::maxn];
            idx Cmidxcol[internal::maxn];
            idx Cmidxrowsubsys_bar[internal::maxn];
            idx Cmidxsubsys[internal::maxn];

            internal::n2multiidx(i, n_subsys_bar, Cdimssubsys_bar,
                                 Cmidxrowsubsys_bar);

            for (idx k = 0; k < n_subsys_bar; ++k) {
                Cmidxrow[Csubsys_bar[k]] = Cmidxrowsubsys_bar[k];
                Cmidxcol[Csubsys_bar[k]] = Cmidxcolsubsys_bar[k];
            }

            Scalar sm = 0;

            for (idx a = 0; a < Dsubsys; ++a) {
                internal::n2multiidx(a, n_subsys, Cdimssubsys, Cmidxsubsys);

                for (idx k = 0; k < n_subsys; ++k) {
                    Cmidxrow[Csubsys[k]] = Cmidxcol[Csubsys[k]] =
                        Cmidxsubsys[k];
                }

                sm += rA(internal::multiidx2n(Cmidxrow, n, Cdims)) *
                      Eigen::numext::conj(
                          rA(internal::multiidx2n(Cmidxcol, n, Cdims)));
            }

            result(i, j) = sm;
        }
    }

    return result;
}

/**
 * @brief Partial trace
 *
 * Partial trace of the multi-partite density matrix \a A over the list
 * \a target of subsystems
 *
 * @param A Eigen expression
 * @param target Subsystem indexes
 * @param dims Dimensions of the multi-partite system
 * @return Partial trace \f$Tr_{subsys}(\cdot)\f$ over the subsystems \a
 * target in a multi-partite system, as a dynamic matrix over the same
 * scalar field as \a A
 */
template <typename Derived>
[[qpp::critical, qpp::parallel]] dyn_mat<typename Derived::Scalar>
ptrace_rho_kq(const Eigen::MatrixBase<Derived>& A,
              const std::vector<idx>& target, const std::vector<idx>& dims) {
    using Scalar = typename Derived::Scalar;
    const auto& rA = A.derived();

    const idx n = dims.size();
    const idx n_subsys = target.size();
    const idx n_subsys_bar = n - n_subsys;

    if (target.size() == dims.size()) {
        dyn_mat<Scalar> result(1, 1);
        result(0, 0) = rA.trace();
        return result;
    }

    if (target.empty()) {
        return rA;
    }

    const idx D = static_cast<idx>(rA.rows());

    idx Dsubsys = 1;
    for (idx i = 0; i < n_subsys; ++i) {
        Dsubsys *= dims[target[i]];
    }
    const idx Dsubsys_bar = D / Dsubsys;

    idx Cdims[internal::maxn];
    idx Csubsys[internal::maxn];
    idx Cdimssubsys[internal::maxn];
    idx Csubsys_bar[internal::maxn];
    idx Cdimssubsys_bar[internal::maxn];

    const std::vector<idx> subsys_bar = fast_complement(target, n);

    std::copy(dims.begin(), dims.end(), Cdims);
    std::copy(target.begin(), target.end(), Csubsys);
    std::copy(subsys_bar.begin(), subsys_bar.end(), Csubsys_bar);

    for (idx i = 0; i < n_subsys; ++i) {
        Cdimssubsys[i] = dims[target[i]];
    }
    for (idx i = 0; i < n_subsys_bar; ++i) {
        Cdimssubsys_bar[i] = dims[subsys_bar[i]];
    }

    dyn_mat<Scalar> result(Dsubsys_bar, Dsubsys_bar);

    for (idx j = 0; j < Dsubsys_bar; ++j) {
        idx Cmidxcolsubsys_bar[internal::maxn];

        internal::n2multiidx(j, n_subsys_bar, Cdimssubsys_bar,
                             Cmidxcolsubsys_bar);

#ifdef QPP_ENABLE_OPENMP
#pragma omp parallel for
#endif // QPP_ENABLE_OPENMP
        for (idx i = 0; i < Dsubsys_bar; ++i) {
            idx Cmidxrow[internal::maxn];
            idx Cmidxcol[internal::maxn];
            idx Cmidxrowsubsys_bar[internal::maxn];
            idx Cmidxsubsys[internal::maxn];

            internal::n2multiidx(i, n_subsys_bar, Cdimssubsys_bar,
                                 Cmidxrowsubsys_bar);

            for (idx k = 0; k < n_subsys_bar; ++k) {
                Cmidxrow[Csubsys_bar[k]] = Cmidxrowsubsys_bar[k];
                Cmidxcol[Csubsys_bar[k]] = Cmidxcolsubsys_bar[k];
            }

            Scalar sm = 0;

            for (idx a = 0; a < Dsubsys; ++a) {
                internal::n2multiidx(a, n_subsys, Cdimssubsys, Cmidxsubsys);

                for (idx k = 0; k < n_subsys; ++k) {
                    Cmidxrow[Csubsys[k]] = Cmidxcol[Csubsys[k]] =
                        Cmidxsubsys[k];
                }

                sm += rA(internal::multiidx2n(Cmidxrow, n, Cdims),
                         internal::multiidx2n(Cmidxcol, n, Cdims));
            }

            result(i, j) = sm;
        }
    }

    return result;
}
} // namespace qpp::internal::kernels::qudit

#endif /* QPP_INTERNAL_KERNELS_QUDIT_PTRACE_HPP_ */
