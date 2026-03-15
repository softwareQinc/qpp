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
 * @file qpp/internal/kernels/qubit/syspermute.hpp
 * @brief Internal highly optimized critical functions for qpp::syspermute()
 */

#ifndef QPP_INTERNAL_KERNELS_QUBIT_SYSPERMUTE_HPP_
#define QPP_INTERNAL_KERNELS_QUBIT_SYSPERMUTE_HPP_

#include <cassert>
#include <vector>
#ifndef NDEBUG
#include <set>
#endif

#include <Eigen/Dense>

#include "qpp/types.hpp"

namespace qpp::internal::kernels::qubit {
/**
 * @brief Qubit state vector subsystem permutation
 *
 * Permutes the subsystems of a qubit state vector. The qubit \a perm[\a i] is
 * permuted to the location \a i.
 *
 * @param A Eigen expression
 * @param perm Permutation
 * @param n Number of qubits
 * @return Permuted qubit system, as a dynamic matrix over the same scalar field
 * as \a A
 */
template <typename Derived>
[[qpp::critical, qpp::parallel]] dyn_mat<typename Derived::Scalar>
syspermute_psi_kq(const Eigen::MatrixBase<Derived>& A,
                  const std::vector<idx>& perm, idx n) {
    const typename Eigen::MatrixBase<Derived>::EvalReturnType& rA = A.derived();

    idx D = static_cast<idx>(rA.rows());

    dyn_mat<typename Derived::Scalar> result(D, 1);

    // Calculate the **inverse permutation** (new_pos -> old_pos).
    // This allows us to construct the new index bit by bit.
    std::vector<idx> inv_perm(n);
    for (idx k = 0; k < n; ++k) {
        inv_perm[perm[k]] = k;
    }

    // Define the bit-permutation logic directly in a block to be used inside
    // the loop to avoid lambda overhead.
    auto permute_bits = [&inv_perm, n](idx old_index) noexcept -> idx {
        std::size_t new_index = 0;
        const std::size_t old_index_sz = static_cast<std::size_t>(old_index);
        for (idx k = 0; k < n; ++k) {
            // k is the new qubit position (new bit position, 0 to n-1)
            // inv_perm[k] is the old qubit position (old bit position, 0 to
            // n-1)
            idx old_pos = inv_perm[k];

            // 1. Extract the bit at old_pos in the old_index
            std::size_t bit = (old_index_sz >> old_pos) & std::size_t{1};

            // 2. Set this bit at the new_pos (k) in the new_index
            new_index |= (bit << k);
        }
        return static_cast<idx>(new_index);
    };

#ifdef QPP_OPENMP
// NOLINTNEXTLINE
#pragma omp parallel for
#endif // QPP_OPENMP
    for (idx i = 0; i < D; ++i) {
        // Assign the original element to the calculated permuted index
        result(permute_bits(i)) = rA.data()[i];
    }

    return result;
}

/**
 * @brief Qubit density matrix subsystem permutation
 *
 * Permutes the subsystems of a qubit density matrix. The qubit \a perm[\a i] is
 * permuted to the location \a i.
 *
 * @param A Eigen expression
 * @param perm Permutation
 * @param n Number of qubits
 * @return Permuted qubit system, as a dynamic matrix over the same scalar field
 * as \a A
 */
template <typename Derived>
[[qpp::critical, qpp::parallel]] dyn_mat<typename Derived::Scalar>
syspermute_rho_kq(const Eigen::MatrixBase<Derived>& A,
                  const std::vector<idx>& perm, idx n) {
    const typename Eigen::MatrixBase<Derived>::EvalReturnType& rA = A.derived();

    idx D = static_cast<idx>(rA.rows());

    dyn_mat<typename Derived::Scalar> result(D, D);

    // Calculate the **inverse permutation** (new_pos -> old_pos).
    // This allows us to construct the new index bit by bit.
    std::vector<idx> inv_perm(n);
    for (idx k = 0; k < n; ++k) {
        inv_perm[perm[k]] = k;
    }

    // Define the bit-permutation logic directly in a block to be used inside
    // the loop to avoid lambda overhead.
    auto permute_bits = [&inv_perm, n](idx old_index) noexcept -> idx {
        std::size_t new_index = 0;
        const std::size_t old_index_sz = static_cast<std::size_t>(old_index);
        for (idx k = 0; k < n; ++k) {
            // k is the new qubit position (new bit position, 0 to n-1)
            // inv_perm[k] is the old qubit position (old bit position, 0 to
            // n-1)
            idx old_pos = inv_perm[k];

            // 1. Extract the bit at old_pos in the old_index
            std::size_t bit = (old_index_sz >> old_pos) & std::size_t{1};

            // 2. Set this bit at the new_pos (k) in the new_index
            new_index |= (bit << k);
        }
        return static_cast<idx>(new_index);
    };

#ifdef QPP_OPENMP
// NOLINTNEXTLINE
#pragma omp parallel for
#endif // QPP_OPENMP
    for (idx i = 0; i < D * D; ++i) {
        idx permuted_linear_index;

        // Density matrix: i = row * D + col
        // Row index is i / D. Col index is i % D.
        idx row_index = i / D;
        idx col_index = i % D;

        // Permute row and column indices separately
        idx new_row_index = permute_bits(row_index);
        idx new_col_index = permute_bits(col_index);

        // The linear index in the permuted system is new_row * D + new_col
        permuted_linear_index = (new_row_index * D) + new_col_index;

        // Assign the original element to the calculated permuted index
        result(permuted_linear_index) = rA.data()[i];
    }

    return result;
}
} // namespace qpp::internal::kernels::qubit

#endif /* QPP_INTERNAL_KERNELS_QUBIT_SYSPERMUTE_HPP_ */
