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

// Fundamental types, can be changed at compile time

/**
 * \brief Integer index
 */
#if defined(QPP_IDX_DEFAULT)
using idx = std::size_t;
#elif defined(QPP_IDX_SHORT)
using idx = short int;
#elif defined(QPP_IDX_INT)
using idx = int;
#elif defined(QPP_IDX_LONG)
using idx = long int;
#elif defined(QPP_IDX_LONG_LONG)
using idx = long long int;
#elif defined(QPP_IDX_USHORT)
using idx = unsigned short int;
#elif defined(QPP_IDX_UINT)
using idx = unsigned int;
#elif defined(QPP_IDX_ULONG)
using idx = unsigned long int;
#elif defined(QPP_IDX_ULONG_LONG)
using idx = unsigned long long int;
#endif
static_assert(std::is_integral_v<idx>, "Type must be integral");
static_assert(sizeof(idx) > 1, "Type must be at least 2 bytes long");

/**
 * \brief Signed big integer
 */
#if defined(QPP_BIGINT_DEFAULT)
using bigint = long long int;
#elif defined(QPP_BIGINT_SHORT)
using bigint = short int;
#elif defined(QPP_BIGINT_INT)
using bigint = int;
#elif defined(QPP_BIGINT_LONG)
using bigint = long int;
#elif defined(QPP_BIGINT_LONG_LONG)
using bigint = long long int;
#endif
static_assert(std::is_integral_v<bigint>, "Type must be integral");
static_assert(std::is_signed_v<bigint>, "Type must be signed");
static_assert(sizeof(bigint) > 1, "Type must be at least 2 bytes long");

/**
 * \brief Floating-point type
 */
#if defined(QPP_FP_DEFAULT)
using realT = double; // default floating-point type
#elif defined(QPP_FP_FLOAT)
using realT = float;
#elif defined(QPP_FP_DOUBLE)
using realT = double;
#elif defined(QPP_FP_LONG_DOUBLE)
using realT = long double;
#endif
static_assert(std::is_floating_point_v<realT>, "Type myst be floating-point");

// The types below are dependent types, please do not change anything below this
// line

/**
 * \brief Unsigned big integer
 */
using ubigint = std::make_unsigned<bigint>::type;
static_assert(std::is_integral_v<ubigint>, "Type must be integral");
static_assert(std::is_unsigned_v<ubigint>, "Type must be unsigned");
static_assert(sizeof(ubigint) > 1, "Type must be at least 2 bytes long");

/**
 * \brief Complex number in realT precision
 */
using cplx = std::complex<realT>;

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

/**
 * \brief Complex (realT precision) dynamic Eigen column vector
 */
using ket = dyn_col_vect<cplx>;

/**
 * \brief Complex (realT precision) dynamic Eigen row vector
 */
using bra = dyn_row_vect<cplx>;

/**
 * \brief Complex (realT precision) dynamic Eigen matrix
 */
using cmat = dyn_mat<cplx>;

/**
 * \brief Real (realT precision) dynamic Eigen matrix
 */
using rmat = dyn_mat<realT>;

/**
 * \brief Quantumly-accessible Random Access Memory (qRAM)
 */
using qram = std::vector<idx>;

namespace internal {
/**
 * \brief Eigen type (ket/density matrix) deduced from the expression Derived
 */
template <typename Derived>
using eval_t =
    std::decay_t<typename Eigen::MatrixBase<Derived>::EvalReturnType>;
} /* namespace internal */

/**
 * \brief Eigen type (ket/density matrix) deduced from the expression Derived
 */
// thanks @antoine-bussy for the suggestion
// https://github.com/softwareQinc/qpp/issues/132#issuecomment-1258360069
template <typename Derived>
using expr_t = Eigen::Matrix<
    typename internal::eval_t<Derived>::Scalar,
    internal::eval_t<Derived>::RowsAtCompileTime == 1 ? 1 : Eigen::Dynamic,
    internal::eval_t<Derived>::ColsAtCompileTime == 1 ? 1 : Eigen::Dynamic,
    internal::eval_t<Derived>::Options,
    internal::eval_t<Derived>::RowsAtCompileTime == 1 ? 1 : Eigen::Dynamic,
    internal::eval_t<Derived>::ColsAtCompileTime == 1 ? 1 : Eigen::Dynamic>;

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
