/*
 * This file is part of Quantum++.
 *
 * Copyright (c) 2013 - 2023 softwareQ Inc. All rights reserved.
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
 * \file types.hpp
 * \brief Type aliases
 */

#ifndef QPP_TYPES_HPP_
#define QPP_TYPES_HPP_

namespace qpp {
/**
 * \brief Non-negative integer index (we use an unsigned type)
 */
using idx = long long int;

/**
 * \brief Big integer
 */
using bigint = long long int;

/**
 * \brief Complex number in double precision
 */
using cplx = std::complex<double>;

/**
 * \brief Complex (double precision) dynamic Eigen column vector
 */
using ket = Eigen::VectorXcd;

/**
 * \brief Complex (double precision) dynamic Eigen row vector
 */
using bra = Eigen::RowVectorXcd;

/**
 * \brief Complex (double precision) dynamic Eigen matrix
 */
using cmat = Eigen::MatrixXcd;

/**
 * \brief Real (double precision) dynamic Eigen matrix
 */
using dmat = Eigen::MatrixXd;

/**
 * \brief Dynamic Eigen matrix over the field specified by \a Scalar
 *
 * Example:
 * \code
 * // type of mat is Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>
 * dyn_mat<float> mat(2, 3);
 * \endcode
 */
template <typename Scalar> // Eigen::MatrixX_type (where type = Scalar)
using dyn_mat = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;

/**
 * \brief Dynamic Eigen column vector over the field specified by \a Scalar
 *
 * Example:
 * \code
 * // type of col_vect is Eigen::Matrix<float, Eigen::Dynamic, 1>
 * dyn_col_vect<float> col_vect(2);
 * \endcode
 */
template <typename Scalar> // Eigen::VectorX_type (where type = Scalar)
using dyn_col_vect = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

/**
 * \brief Dynamic Eigen row vector over the field specified by \a Scalar
 *
 * Example:
 * \code
 * // type of row_vect is Eigen::Matrix<float, 1, Eigen::Dynamic>
 * dyn_row_vect<float> row_vect(3);
 * \endcode
 */
template <typename Scalar> // Eigen::RowVectorX_type (where type = Scalar)
using dyn_row_vect = Eigen::Matrix<Scalar, 1, Eigen::Dynamic>;

namespace internal {
/**
 * \brief Eigen type (ket/density matrix) deduced from the expression Derived
 */
template <typename Derived>
using eval_t =
    std::decay_t<typename Eigen::MatrixBase<Derived>::EvalReturnType>;

/**
 * \brief Detect if the expression Derived is a bra at compile time
 */
template <typename Derived>
bool constexpr is_bra() {
    return (eval_t<Derived>::RowsAtCompileTime == 1);
}

/**
 * \brief Detect if the expression Derived is a ket at compile time
 */
template <typename Derived>
bool constexpr is_ket() {
    return (eval_t<Derived>::ColsAtCompileTime == 1);
}
} /* namespace internal */

/**
 * \brief Eigen type (ket/density matrix) deduced from the expression Derived
 */
// thanks @antoine-bussy for the suggestion
// https://github.com/softwareQinc/qpp/issues/132#issuecomment-1258360069
template <typename Derived>
using expr_t = Eigen::Matrix<typename internal::eval_t<Derived>::Scalar,
                             internal::is_bra<Derived>() ? 1 : Eigen::Dynamic,
                             internal::is_ket<Derived>() ? 1 : Eigen::Dynamic,
                             internal::eval_t<Derived>::Options,
                             internal::is_bra<Derived>() ? 1 : Eigen::Dynamic,
                             internal::is_ket<Derived>() ? 1 : Eigen::Dynamic>;

/**
 * \brief Quantumly-accessible Random Access Memory (qRAM)
 */
using qram = std::vector<idx>;

/**
 * \brief Variant type-matching utility for std::visit
 * \tparam Ts Type list
 */
template <class... Ts>
struct overloaded : Ts... {
    using Ts::operator()...;
};

/**
 * \brief Template deduction rule
 */
template <class... Ts>
overloaded(Ts...) -> overloaded<Ts...>;

} /* namespace qpp */

#endif /* QPP_TYPES_HPP_ */
