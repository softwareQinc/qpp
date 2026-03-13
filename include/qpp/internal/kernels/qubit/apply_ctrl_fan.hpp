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
 * @brief Applies the controlled 1-qubit gate \a A to the qubit \a i of the
 * multi-partite state vector \a state
 *
 * @param state Eigen expression
 * @param A Eigen expression (2x2 matrix)
 * @param i Target subsystem index
 * @param ctrl Vector of control qubit indexes
 * @param shift Binary vector (0/1) where a value of 1 indicates a negated
 * control qubit
 * @param n Number of qubits
 * @return Controlled gate \a A applied to the qubit \a i of \a state
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] expr_t<Derived1>
apply_ctrl_fan_psi(const Eigen::MatrixBase<Derived1>& state,
                   const Eigen::MatrixBase<Derived2>& A,
                   const std::vector<idx>& ctrl, const std::vector<idx>& target,
                   const std::vector<idx>& shift, idx n) {

    const expr_t<Derived1>& rstate = state.derived();
    const dyn_mat<typename Derived2::Scalar>& rA = A.derived();

    // ==========================================
    // ASSERTIONS (Zero runtime overhead in Release)
    // ==========================================
    assert((std::is_same_v<typename Derived1::Scalar,
                           typename Derived2::Scalar>) &&
           "Type mismatch between A and state");
    assert(rA.size() > 0 && rstate.size() > 0 && ctrl.size() > 0 &&
           target.size() > 0 && "Zero size inputs not allowed");
    assert(rA.rows() == rA.cols() && "Gate A must be a square matrix");
    assert(rA.rows() == 2 && "Aggressively optimized for qubits (d=2) only");

    idx D = 1ULL << n; // Total dimension 2^n
    assert(static_cast<idx>(rstate.rows()) == D &&
           "State vector dimension mismatch with 2^n");

    assert(shift.size() == ctrl.size() && "Shift size must match ctrl size");
#ifndef NDEBUG
    for (idx s : shift) {
        assert(s < 2 && "Shift values must be < 2 for qubits");
    }
#endif

    if (D == 1) {
        return rstate; // Trivial 1D state
    }

    // ==========================================
    // BITWISE CONTROL MASK SETUP
    // ==========================================
    idx ctrl_mask = 0;
    idx ctrl_val = 0;
    const idx ctrl_size = static_cast<idx>(ctrl.size());
    for (idx i = 0; i < ctrl_size; ++i) {
        idx bit = 1ULL << (n - 1 - ctrl[i]); // MSB convention index
        ctrl_mask |= bit;
        if (shift[i] == 0) {
            ctrl_val |= bit; // Default to triggering on |1>
        }
    }

    // Gate matrix elements
    auto a00 = rA(0, 0), a01 = rA(0, 1);
    auto a10 = rA(1, 0), a11 = rA(1, 1);

    // ==========================================
    // STATE VECTOR (KET) EVALUATION
    // ==========================================
    dyn_col_vect<typename Derived1::Scalar> result = rstate;

    for (idx t : target) {
        idx target_bit = 1ULL << (n - 1 - t);
        idx fixed_mask = ctrl_mask | target_bit;

        idx free_mask = ((1ULL << n) - 1) & ~fixed_mask;
        idx num_t_active = 1ULL << (n - ctrl_size - 1);

#ifdef QPP_OPENMP
// NOLINTNEXTLINE
#pragma omp parallel for
#endif // QPP_OPENMP
        for (idx i = 0; i < num_t_active; ++i) {
            // Software PDEP bit-scattering
            idx res = 0, val = i, m = free_mask, pos = 0;
            while (m) {
                if (val & (1ULL << pos)) {
                    res |= (m & -m);
                }
                m &= m - 1;
                pos++;
            }
            idx i0 = res | ctrl_val;
            idx i1 = i0 | target_bit;

            auto v0 = result(i0);
            auto v1 = result(i1);
            result(i0) = (a00 * v0) + (a01 * v1);
            result(i1) = (a10 * v0) + (a11 * v1);
        }
    }
    return result;
}

/**
 * @brief Applies the controlled 1-qubit gate \a A to the qubit \a i of the
 * multi-partite density matrix \a state
 *
 * @param state Eigen expression
 * @param A Eigen expression (2x2 matrix)
 * @param i Target subsystem index
 * @param ctrl Vector of control qubit indexes
 * @param shift Binary vector (0/1) where a value of 1 indicates a negated
 * control qubit
 * @param n Number of qubits
 * @return Controlled gate \a A applied to the qubit \a i of \a state
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] expr_t<Derived1>
apply_ctrl_fan_rho(const Eigen::MatrixBase<Derived1>& state,
                   const Eigen::MatrixBase<Derived2>& A,
                   const std::vector<idx>& ctrl, const std::vector<idx>& target,
                   const std::vector<idx>& shift, idx n) {

    const expr_t<Derived1>& rstate = state.derived();
    const dyn_mat<typename Derived2::Scalar>& rA = A.derived();
    idx D = 1ULL << n; // Total dimension 2^n

    // ==========================================
    // ASSERTIONS (Zero runtime overhead in Release)
    // ==========================================
    assert((std::is_same_v<typename Derived1::Scalar,
                           typename Derived2::Scalar>) &&
           "Type mismatch between A and state");
    assert(rA.size() > 0 && rstate.size() > 0 && ctrl.size() > 0 &&
           target.size() > 0 && "Zero size inputs not allowed");
    assert(rA.rows() == rA.cols() && "Gate A must be a square matrix");
    assert(rA.rows() == 2 && "Aggressively optimized for qubits (d=2) only");
    assert(shift.size() == ctrl.size() && "Shift size must match ctrl size");

#ifndef NDEBUG
    assert(static_cast<idx>(rstate.rows()) == D &&
           static_cast<idx>(rstate.cols()) == D &&
           "Density matrix dimension mismatch with 2^n");

    for (idx s : shift) {
        assert(s < 2 && "Shift values must be < 2 for qubits");
    }
#endif

    if (D == 1) {
        return rstate; // Trivial 1D state
    }

    // ==========================================
    // BITWISE CONTROL MASK SETUP
    // ==========================================
    idx ctrl_mask = 0;
    idx ctrl_val = 0;
    const idx ctrl_size = static_cast<idx>(ctrl.size());
    for (idx i = 0; i < ctrl_size; ++i) {
        idx bit = 1ULL << (n - 1 - ctrl[i]); // MSB convention index
        ctrl_mask |= bit;
        if (shift[i] == 0) {
            ctrl_val |= bit; // Default to triggering on |1>
        }
    }

    // Gate matrix elements
    auto a00 = rA(0, 0), a01 = rA(0, 1);
    auto a10 = rA(1, 0), a11 = rA(1, 1);

    // ==========================================
    // DENSITY MATRIX EVALUATION
    // ==========================================
    dyn_mat<typename Derived1::Scalar> result = rstate;

    auto b00 = std::conj(a00), b01 = std::conj(a01);
    auto b10 = std::conj(a10), b11 = std::conj(a11);

    for (idx t : target) {
        idx target_bit = 1ULL << (n - 1 - t);
        idx fixed_mask = ctrl_mask | target_bit;

        idx free_mask = ((1ULL << n) - 1) & ~fixed_mask;
        idx num_t_active = 1ULL << (n - ctrl_size - 1);

        // Precalculate address pairs to share between Row/Col operations
        std::vector<std::pair<idx, idx>> active_pairs(num_t_active);

#ifdef QPP_OPENMP
// NOLINTNEXTLINE
#pragma omp parallel for
#endif // QPP_OPENMP
        for (idx i = 0; i < num_t_active; ++i) {
            idx res = 0, val = i, m = free_mask, pos = 0;
            while (m) {
                if (val & (1ULL << pos)) {
                    res |= (m & -m);
                }
                m &= m - 1;
                pos++;
            }
            idx i0 = res | ctrl_val;
            active_pairs[i] = {i0, i0 | target_bit};
        }

        // 1. Row operations (Column-major cache optimized)
#ifdef QPP_OPENMP
// NOLINTNEXTLINE
#pragma omp parallel for
#endif // QPP_OPENMP
        for (idx c = 0; c < D; ++c) {
            for (idx i = 0; i < num_t_active; ++i) {
                idx i0 = active_pairs[i].first;
                idx i1 = active_pairs[i].second;
                auto v0 = result(i0, c);
                auto v1 = result(i1, c);
                result(i0, c) = (a00 * v0) + (a01 * v1);
                result(i1, c) = (a10 * v0) + (a11 * v1);
            }
        }

        // 2. Column operations
#ifdef QPP_OPENMP
// NOLINTNEXTLINE
#pragma omp parallel for
#endif // QPP_OPENMP
        for (idx i = 0; i < num_t_active; ++i) {
            idx k0 = active_pairs[i].first;
            idx k1 = active_pairs[i].second;
            for (idx r = 0; r < D; ++r) {
                auto v0 = result(r, k0);
                auto v1 = result(r, k1);
                result(r, k0) = (v0 * b00) + (v1 * b01);
                result(r, k1) = (v0 * b10) + (v1 * b11);
            }
        }
    }
    return result;
}
} // namespace qpp::internal::kernels::qubit

#endif /* QPP_INTERNAL_KERNELS_QUBIT_APPLY_CTRL_FAN_HPP_ */
