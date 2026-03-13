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
 * @file qpp/internal/kernels/qubit/ptrace.hpp
 * @brief Internal highly optimized critical functions for qpp::ptrace()
 */

#ifndef QPP_INTERNAL_KERNELS_QUBIT_PTRACE_HPP_
#define QPP_INTERNAL_KERNELS_QUBIT_PTRACE_HPP_

#include <cassert>
#include <vector>
#ifndef NDEBUG
#include <set>
#endif

#include <Eigen/Dense>

#include "qpp/internal/util.hpp"
#include "qpp/types.hpp"

namespace qpp::internal::kernels::qubit {
/**
 * @brief Qubit state vector partial trace
 * @see qpp::ptrace1(), qpp::ptrace2()
 *
 * Partial trace of the multi-partite qubit state vector over the subsystems
 * \a target
 *
 * @param A Eigen expression
 * @param target Subsystem indexes
 * @param n Number of qubits
 * @return Partial trace \f$Tr_{subsys}(\cdot)\f$ over the subsystems \a target
 * in a multi-partite system, as a dynamic matrix over the same scalar field as
 * \a A
 */
template <typename Derived>
[[qpp::critical, qpp::parallel]] dyn_mat<typename Derived::Scalar>
ptrace_psi_kq(const Eigen::MatrixBase<Derived>& A,
              const std::vector<idx>& target, idx n) {
    using Scalar = typename Derived::Scalar;
    const typename Eigen::MatrixBase<Derived>::EvalReturnType& rA = A.derived();

    idx D = static_cast<idx>(rA.rows());
    idx n_subsys = target.size();

    // Fast return for full trace
    if (n_subsys == n) {
        dyn_mat<Scalar> result(1, 1);
        result(0, 0) = (adjoint(rA) * rA).value();
        return result;
    }

    // Fast return for empty trace
    if (n_subsys == 0) {
        return rA * adjoint(rA);
    }

    // --- Optimization Setup ---
    std::vector<idx> subsys_bar = fast_complement(target, n);
    idx n_subsys_bar = subsys_bar.size();

    idx Dsubsys = internal::safe_pow<idx>(2, target.size());
    idx Dsubsys_bar = D / Dsubsys;

    dyn_mat<Scalar> result(Dsubsys_bar, Dsubsys_bar);

    // Precompute global strides: 2^(n - 1 - k)
    std::vector<idx> global_strides(n);
    for (idx k = 0; k < n; ++k) {
        global_strides[k] = static_cast<idx>(1) << (n - 1 - k);
    }

    // Extract strides for Bar (kept) subsystems
    std::vector<idx> bar_strides(n_subsys_bar);
    for (idx k = 0; k < n_subsys_bar; ++k) {
        bar_strides[k] = global_strides[subsys_bar[k]];
    }

    // Pre-calculate Trace Offsets (Sum of strides for target bits)
    std::vector<idx> trace_offsets(Dsubsys);
    for (idx a = 0; a < Dsubsys; ++a) {
        idx offset = 0;
        for (idx k = 0; k < n_subsys; ++k) {
            if ((a >> (n_subsys - 1 - k)) & 1) {
                offset += global_strides[target[k]];
            }
        }
        trace_offsets[a] = offset;
    }
    const idx* trace_offsets_ptr = trace_offsets.data();

    // Lambda: Expand compressed index (0..Dbar-1) to global base offset
    // using bitmasks
    auto expand_bits = [&](idx idx_compressed) -> idx {
        idx expanded = 0;
        for (idx k = 0; k < n_subsys_bar; ++k) {
            if ((idx_compressed >> (n_subsys_bar - 1 - k)) & 1) {
                expanded += bar_strides[k];
            }
        }
        return expanded;
    };

    // Execution Loop
    for (idx j = 0; j < Dsubsys_bar; ++j) {
        idx col_base = expand_bits(j);
#ifdef QPP_OPENMP
// NOLINTNEXTLINE
#pragma omp parallel for
#endif // QPP_OPENMP
        for (idx i = 0; i < Dsubsys_bar; ++i) {
            idx row_base = expand_bits(i);
            Scalar sm = 0;
            // Tight inner loop: iterating over pre-computed offsets
            for (idx k = 0; k < Dsubsys; ++k) {
                idx offset = trace_offsets_ptr[k];
                sm += rA(row_base + offset) * std::conj(rA(col_base + offset));
            }
            result(i, j) = sm;
        }
    }

    return result;
}

/**
 * @brief Qubit density matrix partial trace
 * @see qpp::ptrace1(), qpp::ptrace2()
 *
 * Partial trace of the multi-partite density matrix over the subsystems
 * \a target
 *
 * @param A Eigen expression
 * @param target Subsystem indexes
 * @param n Number of qubits
 * @return Partial trace \f$Tr_{subsys}(\cdot)\f$ over the subsystems \a target
 * in a multi-partite system, as a dynamic matrix over the same scalar field as
 * \a A
 */
template <typename Derived>
[[qpp::critical, qpp::parallel]] dyn_mat<typename Derived::Scalar>
ptrace_rho_kq(const Eigen::MatrixBase<Derived>& A,
              const std::vector<idx>& target, idx n) {
    using Scalar = typename Derived::Scalar;
    const typename Eigen::MatrixBase<Derived>::EvalReturnType& rA = A.derived();

    idx D = static_cast<idx>(rA.rows());
    idx n_subsys = target.size();

    // Fast return for full trace
    if (n_subsys == n) {
        dyn_mat<Scalar> result(1, 1);
        result(0, 0) = rA.trace();
        return result;
    }

    // Fast return for empty trace
    if (n_subsys == 0) {
        return rA;
    }

    // --- Optimization Setup ---
    std::vector<idx> subsys_bar = fast_complement(target, n);
    idx n_subsys_bar = subsys_bar.size();

    idx Dsubsys = internal::safe_pow<idx>(2, target.size());
    idx Dsubsys_bar = D / Dsubsys;

    dyn_mat<Scalar> result(Dsubsys_bar, Dsubsys_bar);

    // Precompute global strides: 2^(n - 1 - k)
    std::vector<idx> global_strides(n);
    for (idx k = 0; k < n; ++k) {
        global_strides[k] = static_cast<idx>(1) << (n - 1 - k);
    }

    // Extract strides for Bar (kept) subsystems
    std::vector<idx> bar_strides(n_subsys_bar);
    for (idx k = 0; k < n_subsys_bar; ++k) {
        bar_strides[k] = global_strides[subsys_bar[k]];
    }

    // Pre-calculate Trace Offsets
    std::vector<idx> trace_offsets(Dsubsys);
    for (idx a = 0; a < Dsubsys; ++a) {
        idx offset = 0;
        for (idx k = 0; k < n_subsys; ++k) {
            if ((a >> (n_subsys - 1 - k)) & 1) {
                offset += global_strides[target[k]];
            }
        }
        trace_offsets[a] = offset;
    }
    const idx* trace_offsets_ptr = trace_offsets.data();

    auto expand_bits = [&](idx idx_compressed) -> idx {
        idx expanded = 0;
        for (idx k = 0; k < n_subsys_bar; ++k) {
            if ((idx_compressed >> (n_subsys_bar - 1 - k)) & 1) {
                expanded += bar_strides[k];
            }
        }
        return expanded;
    };

    // Execution Loop
    for (idx j = 0; j < Dsubsys_bar; ++j) {
        idx col_base = expand_bits(j);
#ifdef QPP_OPENMP
// NOLINTNEXTLINE
#pragma omp parallel for
#endif // QPP_OPENMP
        for (idx i = 0; i < Dsubsys_bar; ++i) {
            idx row_base = expand_bits(i);
            Scalar sm = 0;
            // Tight inner loop: iterating over pre-computed offsets
            for (idx k = 0; k < Dsubsys; ++k) {
                idx offset = trace_offsets_ptr[k];
                sm += rA(row_base + offset, col_base + offset);
            }
            result(i, j) = sm;
        }
    }

    return result;
}
} // namespace qpp::internal::kernels::qubit

#endif /* QPP_INTERNAL_KERNELS_QUBIT_PTRACE_HPP_ */
