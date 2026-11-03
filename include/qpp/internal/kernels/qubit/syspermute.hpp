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
#include <cstddef>
#include <limits>
#include <vector>

#include <Eigen/Dense>

#include "qpp/types.hpp"

namespace qpp::internal::kernels::qubit {
namespace detail {

// Maps each old computational-basis index to its new index. Qubit 0 is the
// most-significant bit, and perm[k] names the old qubit placed at position k.
inline std::vector<idx> syspermute_index_map(const std::vector<idx>& perm,
                                             idx n) {
    assert(n < std::numeric_limits<std::size_t>::digits);
    assert(perm.size() == n);

#ifndef NDEBUG
    std::vector<bool> seen(n, false);
    for (idx k : perm) {
        assert(k < n);
        assert(!seen[k] && "perm must contain distinct subsystem indices");
        seen[k] = true;
    }
#endif

    const idx D = static_cast<idx>(std::size_t{1} << n);
    std::vector<idx> index_map(D);

    for (idx old_index = 0; old_index < D; ++old_index) {
        std::size_t new_index = 0;
        for (idx k = 0; k < n; ++k) {
            const std::size_t old_bit =
                (static_cast<std::size_t>(old_index) >> (n - 1 - perm[k])) &
                std::size_t{1};
            new_index |= old_bit << (n - 1 - k);
        }
        index_map[old_index] = static_cast<idx>(new_index);
    }

    return index_map;
}

} // namespace detail

/**
 * @brief Qubit state vector subsystem permutation
 *
 * Permutes the subsystems of a qubit state vector. The qubit perm[i] is
 * permuted to location i.
 *
 * @param A Eigen column-vector expression
 * @param perm Permutation
 * @param n Number of qubits
 * @return Permuted qubit system, as a dynamic matrix over the same scalar field
 * as A
 */
template <typename Derived>
[[qpp::critical, qpp::parallel]] dyn_mat<typename Derived::Scalar>
syspermute_psi_kq(const Eigen::MatrixBase<Derived>& A,
                  const std::vector<idx>& perm, idx n) {
    const auto& rA = A.derived();
    const auto index_map = detail::syspermute_index_map(perm, n);
    const idx D = static_cast<idx>(index_map.size());

    assert(rA.rows() == static_cast<Eigen::Index>(D));
    assert(rA.cols() == 1);

    dyn_mat<typename Derived::Scalar> result(D, 1);

#ifdef QPP_ENABLE_OPENMP
#pragma omp parallel for
#endif // QPP_ENABLE_OPENMP
    for (idx i = 0; i < D; ++i) {
        result(index_map[i]) = rA(i);
    }

    return result;
}

/**
 * @brief Qubit density matrix subsystem permutation
 *
 * Permutes the subsystems of a qubit density matrix. The qubit perm[i] is
 * permuted to location i.
 *
 * @param A Eigen matrix expression
 * @param perm Permutation
 * @param n Number of qubits
 * @return Permuted qubit system, as a dynamic matrix over the same scalar field
 * as A
 */
template <typename Derived>
[[qpp::critical, qpp::parallel]] dyn_mat<typename Derived::Scalar>
syspermute_rho_kq(const Eigen::MatrixBase<Derived>& A,
                  const std::vector<idx>& perm, idx n) {
    const auto& rA = A.derived();
    const auto index_map = detail::syspermute_index_map(perm, n);
    const idx D = static_cast<idx>(index_map.size());

    assert(rA.rows() == static_cast<Eigen::Index>(D));
    assert(rA.cols() == static_cast<Eigen::Index>(D));

    dyn_mat<typename Derived::Scalar> result(D, D);

#ifdef QPP_ENABLE_OPENMP
#pragma omp parallel for collapse(2)
#endif // QPP_ENABLE_OPENMP
    for (idx i = 0; i < D; ++i) {
        for (idx j = 0; j < D; ++j) {
            result(index_map[i], index_map[j]) = rA(i, j);
        }
    }

    return result;
}

} // namespace qpp::internal::kernels::qubit

#endif /* QPP_INTERNAL_KERNELS_QUBIT_SYSPERMUTE_HPP_ */
