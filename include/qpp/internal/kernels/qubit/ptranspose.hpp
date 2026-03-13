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
 * @file qpp/internal/kernels/qubit/ptranspose.hpp
 * @brief Internal highly optimized critical functions for qpp::ptranspose()
 */

#ifndef QPP_INTERNAL_KERNELS_QUBIT_PTRANSPOSE_HPP_
#define QPP_INTERNAL_KERNELS_QUBIT_PTRANSPOSE_HPP_

#include <cassert>
#include <vector>
#ifndef NDEBUG
#include <set>
#endif

#include <Eigen/Dense>

#include "qpp/types.hpp"

namespace qpp::internal::kernels::qubit {
/**
 * @brief Qubit state vector partial transpose
 *
 * Partial transpose of the multi-partite qubit state vector over the list
 * \a target of subsystems
 *
 * @param A Eigen expression
 * @param target Subsystem indexes
 * @return Partial transpose \f$(\cdot)^{T_{subsys}}\f$ over the subsytems
 * \a target in a qubit multi-partite system, as a dynamic matrix over the same
 * scalar field as \a A
 */
template <typename Derived>
dyn_mat<typename Derived::Scalar> [[qpp::critical, qpp::parallel]]
ptranspose_psi_kq(const Eigen::MatrixBase<Derived>& A,
                  const std::vector<idx>& target, idx n) {
    using scalar_t = typename Derived::Scalar;
    const auto& rA = A.derived();

    // Input Validation
    assert(internal::check_nonzero_size(rA) && "A"); // Zero-size
    assert(internal::check_cvector(rA) && "A must be a column vector (ket)");
    assert(n > 0 && "n must be > 0");

    // target must be valid w.r.t. n and contain unique indices
#ifndef NDEBUG
    for (idx k : target) {
        assert(k < n && "target must be valid w.r.t. n");
    }

    auto tmp = target;
    std::sort(tmp.begin(), tmp.end());
    auto it = std::adjacent_find(tmp.begin(), tmp.end());
    assert(it == tmp.end() && "target contains duplicate subsystem indices");

    // D must equal 2^n for qubits
    idx D_expected = (1ULL << n);
    assert(static_cast<idx>(rA.rows()) == D_expected && "A/n size mismatch");
#endif

    idx D = static_cast<idx>(rA.rows());
    idx n_subsys = target.size();

    // Trivial cases
    if (n_subsys == 0) {
        // no partial transpose requested
        return rA * rA.adjoint();
    }
    if (n_subsys == n) {
        // partial transpose over all subsystems == full transpose of density
        // matrix
        return (rA * rA.adjoint()).transpose();
    }

    // build mask M: bit k set means subsystem k is transposed
    idx M = 0;
    for (idx k : target) {
        M |= (1ULL << k);
    }

    dyn_mat<scalar_t> result(D, D);
    result.setZero(); // ensure initialized

    // alias psi accessor
    auto psi = [&](idx t) -> scalar_t const& { return rA(t); };

    // worker lambda computes rho_PT(i,j) = psi[i_T] * conj(psi[j_T])
    auto worker = [&](idx i, idx j) noexcept -> scalar_t {
        idx i_T = (i & (~M)) | (j & M);
        idx j_T = (j & (~M)) | (i & M);
        return conj(psi(i_T)) * psi(j_T);
    };

#ifdef QPP_OPENMP
// NOLINTNEXTLINE
#pragma omp parallel for collapse(2)
#endif // QPP_OPENMP
    for (idx i = 0; i < D; ++i) {
        for (idx j = 0; j < D; ++j) {
            result(i, j) = worker(i, j);
        }
    }

    return result;
}

/**
 * @brief Qubit density matrix partial transpose
 *
 * Partial transpose of the multi-partite qubit density matrix over the list
 * \a target of subsystems
 *
 * @param A Eigen expression
 * @param target Subsystem indexes
 * @return Partial transpose \f$(\cdot)^{T_{subsys}}\f$ over the subsytems
 * \a target in a qubit multi-partite system, as a dynamic matrix over the same
 * scalar field as \a A
 */
template <typename Derived>
dyn_mat<typename Derived::Scalar> [[qpp::critical, qpp::parallel]]
ptranspose_rho_kq(const Eigen::MatrixBase<Derived>& A,
                  const std::vector<idx>& target, idx n) {
    const typename Eigen::MatrixBase<Derived>::EvalReturnType& rA = A.derived();

    // Input Validation
    assert(internal::check_nonzero_size(rA) &&
           "A"); // Corresponds to exception::ZeroSize
    // check_dims(dims) removed, replaced by n > 0 check if necessary, but D
    // implies this check_subsys_match_dims simplified to check target is valid
    // w.r.t. n
    assert(std::all_of(target.cbegin(), target.cend(),
                       [&](idx k) { return k < n; }) &&
           "target must be valid w.r.t. n");
    assert(internal::check_square_mat(rA) &&
           "A must be a square matrix (density "
           "matrix)"); // Corresponds to
                       // check_square_mat/MatrixNotSquareNorCvector
    // check_dims_match_mat replaced by check that D = 2^n
#ifndef NDEBUG
    idx D_expected = (1ULL << n);
    assert(static_cast<idx>(rA.rows()) == D_expected && "A/n size mismatch");
#endif

    idx D = static_cast<idx>(rA.rows());
    idx n_subsys = target.size();

    // Trivial cases
    if (n_subsys == n) { // Use n instead of dims.size()
        return rA.transpose();
    }
    if (target.empty()) {
        return rA;
    }

    dyn_mat<typename Derived::Scalar> result(D, D);

    // OPTIMIZATION: Precompute the partial transpose mask (M)
    // The mask M has a '1' at the bit positions (subsystem indices) specified
    // in target. For qubits, the subsystem index is directly the bit position.
    idx M = 0;
    for (idx k : target) {
        M |= (1ULL << k);
    }

    auto worker = [&](idx i, idx j) noexcept -> typename Derived::Scalar {
        // i' = (i & ~M) | (j & M)
        idx i_T = (i & (~M)) | (j & M);

        // j' = (j & ~M) | (i & M)
        idx j_T = (j & (~M)) | (i & M);

        // result(i, j) = A_j'i' (Partial transpose swaps row and column indices
        // of non-target blocks)
        return rA(j_T, i_T);
    }; /* end worker */

#ifdef QPP_OPENMP
// NOLINTNEXTLINE
#pragma omp parallel for collapse(2)
#endif // QPP_OPENMP
    for (idx i = 0; i < D; ++i) {
        for (idx j = 0; j < D; ++j) {
            result(i, j) = worker(i, j);
        }
    }

    return result;
}
} // namespace qpp::internal::kernels::qubit

#endif /* QPP_INTERNAL_KERNELS_QUBIT_PTRANSPOSE_HPP_ */
