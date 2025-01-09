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
 * \file qpp/traits.hpp
 * \brief Type traits
 */

#ifndef QPP_TRAITS_HPP_
#define QPP_TRAITS_HPP_

#include <type_traits>

#include <Eigen/Dense>

namespace qpp {
/**
 * \brief Checks whether the type is compatible with an STL-like iterable
 * container
 * \see qpp::is_iterable_v
 *
 * Provides the constant member \a value which is equal to \a true, if \a T is
 * compatible with an iterable container, i.e., provides at least \a begin()
 * and \a end() member functions and allows de-referencing. Otherwise, \a value
 * is equal to \a false.
 */
// silence g++4.8.x warning about non-virtual destructor in inherited class
#if defined(__GNUC__) && !defined(__clang__) && !defined(__INTEL_COMPILER) &&  \
    (__GNUC__ == 4) && (__GNUC_MINOR__ == 8)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Weffc++"
#endif
template <typename T, typename = void>
struct is_iterable : std::false_type {};
#if defined(__GNUC__) && !defined(__clang__) && !defined(__INTEL_COMPILER) &&  \
    (__GNUC__ == 4) && (__GNUC_MINOR__ == 8)
#pragma GCC diagnostic pop
#endif

/**
 * \brief Checks whether the type is compatible with an STL-like iterable
 * container, specialization for STL-like iterable containers
 */
// silence g++4.8.x warning about non-virtual destructor in inherited class
#if defined(__GNUC__) && !defined(__clang__) && !defined(__INTEL_COMPILER) &&  \
    (__GNUC__ == 4) && (__GNUC_MINOR__ == 8)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Weffc++"
#endif
template <typename T>
struct is_iterable<T, std::void_t<decltype(std::declval<T>().begin()),
                                  decltype(std::declval<T>().end()),
                                  decltype(*(std::declval<T>().begin()))>>
    : std::true_type {};
#if defined(__GNUC__) && !defined(__clang__) && !defined(__INTEL_COMPILER) &&  \
    (__GNUC__ == 4) && (__GNUC_MINOR__ == 8)
#pragma GCC diagnostic pop
#endif

/**
 * \brief Checks whether the type is compatible with an STL-like iterable
 * container, helper variable template
 * \see qpp::is_iterable
 */
template <typename T>
inline constexpr bool is_iterable_v = is_iterable<T>::value;

/**
 * \brief Checks whether the type is an Eigen matrix expression
 * \see qpp::is_matrix_expression_v
 *
 * Provides the constant member \a value which is equal to \a true, if the type
 * is an Eigen matrix expression of type \a Eigen::MatrixBase<Derived>.
 * Otherwise, \a value is equal to \a false.
 */
// thanks to @davidhigh https://stackoverflow.com/a/40293333/3093378
// silence g++4.8.x warning about non-virtual destructor in inherited class
#if defined(__GNUC__) && !defined(__clang__) && !defined(__INTEL_COMPILER) &&  \
    (__GNUC__ == 4) && (__GNUC_MINOR__ == 8)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Weffc++"
#endif
template <typename Derived>
struct is_matrix_expression
    : std::is_base_of<Eigen::MatrixBase<std::decay_t<Derived>>,
                      std::decay_t<Derived>> {};
#if defined(__GNUC__) && !defined(__clang__) && !defined(__INTEL_COMPILER) &&  \
    (__GNUC__ == 4) && (__GNUC_MINOR__ == 8)
#pragma GCC diagnostic pop
#endif

/**
 * \brief Checks whether the type is an Eigen matrix expression, helper variable
 * template
 * \see qpp::is_matrix_expression
 */
template <typename T>
inline constexpr bool is_matrix_expression_v = is_matrix_expression<T>::value;

/**
 * \brief Checks whether the type is a complex type
 * \see qpp::is_complex_v
 *
 * Provides the constant member \a value which is equal to \a true, if the type
 * is a complex type, i.e., \a std::complex<T>
 */
// silence g++4.8.x warning about non-virtual destructor in inherited class
#if defined(__GNUC__) && !defined(__clang__) && !defined(__INTEL_COMPILER) &&  \
    (__GNUC__ == 4) && (__GNUC_MINOR__ == 8)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Weffc++"
#endif
template <typename T>
struct is_complex : std::false_type {};
#if defined(__GNUC__) && !defined(__clang__) && !defined(__INTEL_COMPILER) &&  \
    (__GNUC__ == 4) && (__GNUC_MINOR__ == 8)
#pragma GCC diagnostic pop
#endif

/**
 * \brief Checks whether the type is a complex number type, specialization for
 * complex types
 */
// silence g++4.8.x warning about non-virtual destructor in inherited class
#if defined(__GNUC__) && !defined(__clang__) && !defined(__INTEL_COMPILER) &&  \
    (__GNUC__ == 4) && (__GNUC_MINOR__ == 8)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Weffc++"
#endif
template <typename T>
struct is_complex<std::complex<T>> : std::true_type {};
#if defined(__GNUC__) && !defined(__clang__) && !defined(__INTEL_COMPILER) &&  \
    (__GNUC__ == 4) && (__GNUC_MINOR__ == 8)
#pragma GCC diagnostic pop
#endif

/**
 * \brief Checks whether the type is a complex number type, helper variable
 * template
 * \see qpp::is_complex
 */
template <typename T>
inline constexpr bool is_complex_v = is_complex<T>::value;

namespace internal {
/**
 * \brief Eigen type (ket/bra/density matrix) deduced from the expression
 * Derived
 */
template <typename Derived>
using eval_t =
    std::decay_t<typename Eigen::MatrixBase<Derived>::EvalReturnType>;
} /* namespace internal */

/**
 * \brief Detect if the expression Derived is a row vector (bra) at compile time
 */
template <typename Derived>
bool constexpr is_bra_v() {
    return (internal::eval_t<Derived>::RowsAtCompileTime == 1);
}

/**
 * \brief Detect if the expression Derived is a column vector(ket) at compile
 * time
 */
template <typename Derived>
bool constexpr is_ket_v() {
    return (internal::eval_t<Derived>::ColsAtCompileTime == 1);
}

} /* namespace qpp */

#endif /* QPP_TRAITS_HPP_ */
