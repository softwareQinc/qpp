/*
 * This file is part of Quantum++.
 *
 * Copyright (c) 2017 - 2025 softwareQ Inc. All rights reserved.
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
 * \file qpp/internal/speedup.hpp
 * \brief Internal highly optimized functions
 */

#ifndef QPP_INTERNAL_SPEEDUP_HPP_
#define QPP_INTERNAL_SPEEDUP_HPP_

#include <vector>

#include <Eigen/Dense>

#include "qpp/types.hpp"

namespace qpp {
namespace internal {
/**
 * \brief Applies the 1-qubit gate \a A to the qubit \a i of the
 * multi-partite state vector \a state
 *
 * \param state Eigen expression
 * \param A Eigen expression
 * \param i Subsystem index where the gate \a A is applied
 * \paran n Number of qubits
 * \return Gate \a A applied to the qubit \a target of \a state
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] expr_t<Derived1>
apply_psi_1q(const Eigen::MatrixBase<Derived1>& state,
             const Eigen::MatrixBase<Derived2>& A, idx i, idx n) {}

/**
 * \brief Applies the 2-qubit gate \a A to the qubits \a i and \a j of the
 * multi-partite state vector \a state
 *
 * \param state Eigen expression
 * \param A Eigen expression
 * \param i Subsystem index where the gate \a A is applied
 * \param j Subsystem index where the gate \a A is applied
 * \paran n Number of qubits
 * \return Gate \a A applied to the qubit \a target of \a state
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] expr_t<Derived1>
apply_psi_2q(const Eigen::MatrixBase<Derived1>& state,
             const Eigen::MatrixBase<Derived2>& A, idx i, idx j, idx n) {}

/**
 * \brief Applies the 3-qubit gate \a A to the qubit \a i, \a j, and \a k of the
 * multi-partite state vector \a state
 *
 * \param state Eigen expression
 * \param A Eigen expression
 * \param i Subsystem index where the gate \a A is applied
 * \param j Subsystem index where the gate \a A is applied
 * \param k Subsystem index where the gate \a A is applied
 * \paran n Number of qubits
 * \return Gate \a A applied to the qubit \a target of \a state
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] expr_t<Derived1>
apply_psi_3q(const Eigen::MatrixBase<Derived1>& state,
             const Eigen::MatrixBase<Derived2>& A, idx i, idx j, idx k, idx n) {
}

/**
 * \brief Applies the multi-qubit gate \a A to the part \a target of the
 * multi-partite state vector \a state
 *
 * \param state Eigen expression
 * \param A Eigen expression
 * \param target Subsystem indexes where the gate \a A is applied
 * \return Gate \a A applied to the part \a target of \a state
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] expr_t<Derived1>
apply_psi_multi(const Eigen::MatrixBase<Derived1>& state,
                const Eigen::MatrixBase<Derived2>& A,
                const std::vector<idx>& target) {}

/**
 * \brief Applies the 1-qubit gate \a A to the qubit \a i of the
 * multi-partite density matrix \a state
 *
 * \param state Eigen expression
 * \param A Eigen expression
 * \param i Subsystem index where the gate \a A is applied
 * \paran n Number of qubits
 * \return Gate \a A applied to the qubit \a target of \a state
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] expr_t<Derived1>
apply_rho_1q(const Eigen::MatrixBase<Derived1>& state,
             const Eigen::MatrixBase<Derived2>& A, idx i, idx n) {}

/**
 * \brief Applies the 2-qubit gate \a A to the qubits \a i and \a j of the
 * multi-partite density matrix \a state
 *
 * \param state Eigen expression
 * \param A Eigen expression
 * \param i Subsystem index where the gate \a A is applied
 * \param j Subsystem index where the gate \a A is applied
 * \paran n Number of qubits
 * \return Gate \a A applied to the qubit \a target of \a state
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] expr_t<Derived1>
apply_rho_2q(const Eigen::MatrixBase<Derived1>& state,
             const Eigen::MatrixBase<Derived2>& A, idx i, idx j, idx n) {}

/**
 * \brief Applies the 3-qubit gate \a A to the qubit \a i, \a j, and \a k of the
 * multi-partite density matrix \a state
 *
 * \param state Eigen expression
 * \param A Eigen expression
 * \param i Subsystem index where the gate \a A is applied
 * \param j Subsystem index where the gate \a A is applied
 * \param k Subsystem index where the gate \a A is applied
 * \paran n Number of qubits
 * \return Gate \a A applied to the qubit \a target of \a state
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] expr_t<Derived1>
apply_rho_3q(const Eigen::MatrixBase<Derived1>& state,
             const Eigen::MatrixBase<Derived2>& A, idx i, idx j, idx k, idx n) {
}

/**
 * \brief Applies the multi-qubit gate \a A to the part \a target of the
 * multi-partite density matrix \a state
 *
 * \param state Eigen expression
 * \param A Eigen expression
 * \param target Subsystem indexes where the gate \a A is applied
 * \return Gate \a A applied to the part \a target of \a state
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] expr_t<Derived1>
apply_rho_multi(const Eigen::MatrixBase<Derived1>& state,
                const Eigen::MatrixBase<Derived2>& A,
                const std::vector<idx>& target) {}

} /* namespace internal */
} /* namespace qpp */

#endif /* QPP_INTERNAL_SPEEDUP_HPP_ */
