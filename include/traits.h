/*
 * This file is part of Quantum++.
 *
 * MIT License
 *
 * Copyright (c) 2013 - 2018 Vlad Gheorghiu (vgheorgh@gmail.com)
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
* \file traits.h
* \brief Type traits
*/

#ifndef TRAITS_H_
#define TRAITS_H_

namespace qpp {
// Citing from http://en.cppreference.com/w/cpp/types/void_t:
// "Until CWG 1558 (a C++14 defect), unused parameters in alias templates were
// not guaranteed to ensure SFINAE and could be ignored, so earlier compilers
// require a more complex definition of void_t, such as:"
/**
* \brief Helper for qpp::to_void<> alias template
* \see qpp::to_void<>
*/
template <typename... Ts>
struct make_void {
    typedef void type;
};

/**
* \brief Alias template that implements the proposal for void_t
*
* \see http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2014/n3911
*/
template <typename... Ts>
using to_void = typename make_void<Ts...>::type;

/**
* \brief Checks whether \a T is compatible with an STL-like iterable container
*
* Provides the constant member \a value which is equal to \a true,
* if \a T is compatible with an iterable container, i.e. provides at least
* \a begin() and \a end() member functions.
* Otherwise, \a value is equal to \a false.
*/
// silence g++4.8.x warning about non-virtual destructor in inherited class
#if ((__GNUC__ == 4) && (__GNUC_MINOR__ == 8) && !__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Weffc++"
#endif
template <typename T, typename = void>
struct is_iterable : std::false_type {};
#if ((__GNUC__ == 4) && (__GNUC_MINOR__ == 8) && !__clang__)
#pragma GCC diagnostic pop
#endif

/**
* \brief Checks whether \a T is compatible with an STL-like iterable container,
* specialization for STL-like iterable containers
*/
// silence g++4.8.x warning about non-virtual destructor in inherited class
#if ((__GNUC__ == 4) && (__GNUC_MINOR__ == 8) && !__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Weffc++"
#endif
template <typename T>
struct is_iterable<
    T, to_void<decltype(std::declval<T>().begin()),
               decltype(std::declval<T>().end()), typename T::value_type>>
    : std::true_type {};
#if ((__GNUC__ == 4) && (__GNUC_MINOR__ == 8) && !__clang__)
#pragma GCC diagnostic pop
#endif

/**
* \brief Checks whether the type is an Eigen matrix expression
*
* Provides the constant member \a value which is equal to \a true,
* if the type is an Eigen matrix expression of type \a
* Eigen::MatrixBase<Derived>.
* Otherwise, \a value is equal to \a false.
*/
// thanks to @davidhigh http://stackoverflow.com/a/40293333/3093378
// silence g++4.8.x warning about non-virtual destructor in inherited class
#if ((__GNUC__ == 4) && (__GNUC_MINOR__ == 8) && !__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Weffc++"
#endif
template <typename Derived>
struct is_matrix_expression
    : std::is_base_of<Eigen::MatrixBase<typename std::decay<Derived>::type>,
                      typename std::decay<Derived>::type> {};
#if ((__GNUC__ == 4) && (__GNUC_MINOR__ == 8) && !__clang__)
#pragma GCC diagnostic pop
#endif

/**
* \brief Checks whether the type is a complex type
*
* Provides the constant member \a value which is equal to \a true,
* if the type is a complex type, i.e. \a std::complex<T>
*/
// silence g++4.8.x warning about non-virtual destructor in inherited class
#if ((__GNUC__ == 4) && (__GNUC_MINOR__ == 8) && !__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Weffc++"
#endif
template <typename T>
struct is_complex : std::false_type {};
#if ((__GNUC__ == 4) && (__GNUC_MINOR__ == 8) && !__clang__)
#pragma GCC diagnostic pop
#endif

/**
* \brief \brief Checks whether the type is a complex number type,
* specialization for complex types
*/
// silence g++4.8.x warning about non-virtual destructor in inherited class
#if ((__GNUC__ == 4) && (__GNUC_MINOR__ == 8) && !__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Weffc++"
#endif
template <typename T>
struct is_complex<std::complex<T>> : std::true_type {};
#if ((__GNUC__ == 4) && (__GNUC_MINOR__ == 8) && !__clang__)
#pragma GCC diagnostic pop
#endif

} /* namespace qpp */

#endif /* TRAITS_H_ */
