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
 * @file qpp/internal/kernels/qubit/apply_ctrl_fan.hpp
 * @brief Internal highly optimized critical functions for qpp::applyCTRL_fan()
 */

#ifndef QPP_INTERNAL_KERNELS_QUBIT_APPLY_CTRL_FAN_HPP_
#define QPP_INTERNAL_KERNELS_QUBIT_APPLY_CTRL_FAN_HPP_

#include <cassert>
#include <optional>
#include <vector>
#ifndef NDEBUG
#include <set>
#endif

#include <Eigen/Dense>

#include "qpp/types.hpp"

namespace qpp::internal::kernels::qubit {
/**
 * @brief Applies the single qubit controlled-gate \a A with multiple
 * control qubits listed in \a ctrl to the part \a target of the multi-partite
 * state vector \a state in-place, i.e., CTRL-A-A-...-A
 *
 * @param state Eigen expression (modified in-place)
 * @param A Eigen expression (2x2 matrix)
 * @param ctrl Vector of control qubit indexes
 * @param target Vector of target subsystem indexes
 * @param shift Binary vector (0/1) where a value of 1 indicates a negated
 * control
 * @param n Number of qubits
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] void apply_ctrl_fan_psi_inplace(
    Eigen::MatrixBase<Derived1>& state, const Eigen::MatrixBase<Derived2>& A,
    const std::vector<idx>& ctrl, const std::vector<idx>& target,
    const std::vector<idx>& shift, idx n) {
    const dyn_mat<typename Derived2::Scalar>& rA = A.derived();

    // Assertions
    assert(rA.rows() == 2 && rA.cols() == 2 && "Gate A must be a 2x2 matrix");
    assert(ctrl.size() == shift.size() && "Shift size must match ctrl size");

    idx D = static_cast<idx>(std::size_t{1} << n);
    if (D <= 1 || target.empty()) {
        return;
    }

    // Bitwise Control Mask Setup
    idx ctrl_mask = 0;
    idx ctrl_val = 0;
    const idx ctrl_size = static_cast<idx>(ctrl.size());
    for (idx i = 0; i < ctrl_size; ++i) {
        idx bit = static_cast<idx>(std::size_t{1} << (n - 1 - ctrl[i]));
        ctrl_mask |= bit;
        if (shift[i] == 0) {
            ctrl_val |= bit; // Triggers on |1>
        }
    }

    // Gate matrix elements
    const auto a00 = rA(0, 0), a01 = rA(0, 1);
    const auto a10 = rA(1, 0), a11 = rA(1, 1);

    // Apply gate A to each target qubit sequentially under the same control
    for (idx t : target) {
        idx target_bit = static_cast<idx>(std::size_t{1} << (n - 1 - t));
        idx fixed_mask = ctrl_mask | target_bit;
        idx free_mask = static_cast<idx>(((std::size_t{1} << n) - 1) &
                                         ~static_cast<std::size_t>(fixed_mask));

        idx num_t_active =
            static_cast<idx>(std::size_t{1} << (n - ctrl_size - 1));

#ifdef QPP_OPENMP
// NOLINTNEXTLINE
#pragma omp parallel for
#endif // QPP_OPENMP
        for (idx i = 0; i < num_t_active; ++i) {
            // Software PDEP bit-scattering to find indices that match the rest
            // bits
            std::size_t res = 0;
            std::size_t val = static_cast<std::size_t>(i);
            std::size_t m = static_cast<std::size_t>(free_mask);
            std::size_t pos = 0;
            while (m) {
                if (val & (std::size_t{1} << pos)) {
                    res |= (m & (~m + 1));
                }
                m &= m - 1;
                pos++;
            }

            idx i0 = static_cast<idx>(res) | ctrl_val;
            idx i1 = i0 | target_bit;

            auto v0 = state(i0);
            auto v1 = state(i1);
            state(i0) = (a00 * v0) + (a01 * v1);
            state(i1) = (a10 * v0) + (a11 * v1);
        }
    }
}

/**
 * @brief Applies the single qubit controlled-gate \a A with multiple
 * control qubits listed in \a ctrl to the part \a target of the multi-partite
 * state vector \a state
 *
 * @return Controlled gate \a A applied to the targets of \a state
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] expr_t<Derived1>
apply_ctrl_fan_psi(const Eigen::MatrixBase<Derived1>& state,
                   const Eigen::MatrixBase<Derived2>& A,
                   const std::vector<idx>& ctrl, const std::vector<idx>& target,
                   const std::vector<idx>& shift, idx n) {
    // Functional deep copy
    expr_t<Derived1> result = state.derived();

    // Delegate to in-place implementation
    apply_ctrl_fan_psi_inplace(result, A, ctrl, target, shift, n);

    return result;
}
/**
 * @brief Applies the single qubit controlled-gate \a A with multiple
 * control qubits listed in \a ctrl to the part \a target of the multi-partite
 * density matrix \a state in-place, i.e., CTRL-A-A-...-A
 *
 * @param state Eigen expression (modified in-place)
 * @param A Eigen expression (2x2 matrix)
 * @param ctrl Vector of control qubit indexes
 * @param target Vector of target subsystem indexes
 * @param shift Binary vector (0/1) where a value of 1 indicates a negated
 * control
 * @param n Number of qubits
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] void apply_ctrl_fan_rho_inplace(
    Eigen::MatrixBase<Derived1>& state, const Eigen::MatrixBase<Derived2>& A,
    const std::vector<idx>& ctrl, const std::vector<idx>& target,
    const std::vector<idx>& shift, idx n) {
    const dyn_mat<typename Derived2::Scalar>& rA = A.derived();
    idx D = static_cast<idx>(std::size_t{1} << n);

    // Assertions
    assert(rA.rows() == 2 && rA.cols() == 2 && "Gate A must be a 2x2 matrix");
    assert(static_cast<idx>(state.rows()) == D &&
           static_cast<idx>(state.cols()) == D && "State must be 2^n x 2^n");
    assert(shift.size() == ctrl.size() && "Shift size must match ctrl size");

    if (D <= 1 || target.empty()) {
        return;
    }

    // Bitwise Control Mask Setup
    idx ctrl_mask = 0;
    idx ctrl_val = 0;
    const idx ctrl_size = static_cast<idx>(ctrl.size());
    for (idx i = 0; i < ctrl_size; ++i) {
        idx bit = static_cast<idx>(std::size_t{1} << (n - 1 - ctrl[i]));
        ctrl_mask |= bit;
        if (shift[i] == 0) {
            ctrl_val |= bit; // Triggers on |1>
        }
    }

    // Gate matrix elements and their conjugates
    const auto a00 = rA(0, 0), a01 = rA(0, 1);
    const auto a10 = rA(1, 0), a11 = rA(1, 1);
    const auto b00 = std::conj(a00), b01 = std::conj(a01);
    const auto b10 = std::conj(a10), b11 = std::conj(a11);

    for (idx t : target) {
        idx target_bit = static_cast<idx>(std::size_t{1} << (n - 1 - t));
        idx fixed_mask = ctrl_mask | target_bit;
        idx free_mask = static_cast<idx>(((std::size_t{1} << n) - 1) &
                                         ~static_cast<std::size_t>(fixed_mask));
        idx num_t_active =
            static_cast<idx>(std::size_t{1} << (n - ctrl_size - 1));

        // Precalculate address pairs
        std::vector<std::pair<idx, idx>> active_pairs(num_t_active);
#ifdef QPP_OPENMP
// NOLINTNEXTLINE
#pragma omp parallel for
#endif
        for (idx i = 0; i < num_t_active; ++i) {
            std::size_t res = 0;
            std::size_t val = static_cast<std::size_t>(i);
            std::size_t m = static_cast<std::size_t>(free_mask);
            std::size_t pos = 0;
            while (m) {
                if (val & (std::size_t{1} << pos)) {
                    res |= (m & (~m + 1));
                }
                m &= m - 1;
                pos++;
            }
            idx i0 = static_cast<idx>(res) | ctrl_val;
            active_pairs[i] = {i0, i0 | target_bit};
        }

        // 1. Row operations (Left application: A * rho)
#ifdef QPP_OPENMP
// NOLINTNEXTLINE
#pragma omp parallel for
#endif
        for (idx c = 0; c < D; ++c) {
            for (idx i = 0; i < num_t_active; ++i) {
                idx i0 = active_pairs[i].first;
                idx i1 = active_pairs[i].second;
                auto v0 = state(i0, c);
                auto v1 = state(i1, c);
                state(i0, c) = (a00 * v0) + (a01 * v1);
                state(i1, c) = (a10 * v0) + (a11 * v1);
            }
        }

        // 2. Column operations (Right application: rho * A_dagger)
#ifdef QPP_OPENMP
// NOLINTNEXTLINE
#pragma omp parallel for
#endif
        for (idx i = 0; i < num_t_active; ++i) {
            idx k0 = active_pairs[i].first;
            idx k1 = active_pairs[i].second;
            for (idx r = 0; r < D; ++r) {
                auto v0 = state(r, k0);
                auto v1 = state(r, k1);
                state(r, k0) = (v0 * b00) + (v1 * b01);
                state(r, k1) = (v0 * b10) + (v1 * b11);
            }
        }
    }
}

/**
 * @brief Applies the single qubit controlled-gate \a A with multiple
 * control qubits listed in \a ctrl to the part \a target of the multi-partite
 * density matrix \a state
 *
 * @return Controlled gate \a A applied to the targets of \a state
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] expr_t<Derived1>
apply_ctrl_fan_rho(const Eigen::MatrixBase<Derived1>& state,
                   const Eigen::MatrixBase<Derived2>& A,
                   const std::vector<idx>& ctrl, const std::vector<idx>& target,
                   const std::vector<idx>& shift, idx n) {
    // Functional deep copy
    expr_t<Derived1> result = state.derived();

    // Delegate to in-place implementation
    apply_ctrl_fan_rho_inplace(result, A, ctrl, target, shift, n);

    return result;
}
} // namespace qpp::internal::kernels::qubit

#endif /* QPP_INTERNAL_KERNELS_QUBIT_APPLY_CTRL_FAN_HPP_ */
