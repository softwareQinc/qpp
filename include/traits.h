/*
 * Quantum++
 *
 * Copyright (c) 2013 - 2017 Vlad Gheorghiu (vgheorgh@gmail.com)
 *
 * This file is part of Quantum++.
 *
 * Quantum++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Quantum++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Quantum++.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
* \file traits.h
* \brief Type traits
*/

#ifndef TRAITS_H_
#define TRAITS_H_

namespace qpp
{
/**
* \brief Alias template that implements the proposal for void_t
*
* \see http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2014/n3911
*/
// Citing from http://en.cppreference.com/w/cpp/types/void_t:
// "Until CWG 1558 (a C++14 defect), unused parameters in alias templates were
// not guaranteed to ensure SFINAE and could be ignored, so earlier compilers
// require a more complex definition of void_t, such as:" 
template<typename... Ts> struct make_void { typedef void type;};
template<typename... Ts> using to_void = typename make_void<Ts...>::type;

/**
* \brief Checks whether \a T is compatible with an STL-like iterable container
*
* Provides the constant member \a value which is equal to \a true,
* if \a T is compatible with an iterable container, i.e. provides at least
* \a begin() and \a end() member functions.
* Otherwise, \a value is equal to \a false.
*/
// silence g++4.8.x warning about non-virtual destructor in inherited class
#if ((__GNUC__ == 4) && (__GNUC_MINOR__ == 8)  && !__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Weffc++"
#endif
template<typename T, typename = void>
struct is_iterable : std::false_type
{
};
#if ((__GNUC__ == 4) && (__GNUC_MINOR__ == 8)  && !__clang__)
#pragma GCC diagnostic pop
#endif

/**
* \brief Checks whether \a T is compatible with an STL-like iterable container,
* specialization for STL-like iterable containers
*/
// silence g++4.8.x warning about non-virtual destructor in inherited class
#if ((__GNUC__ == 4) && (__GNUC_MINOR__ == 8)  && !__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Weffc++"
#endif
template<typename T>
struct is_iterable<T,
        to_void<decltype(std::declval<T>().begin()),
                decltype(std::declval<T>().end()),
                typename T::value_type
        >> : std::true_type
{
};
#if ((__GNUC__ == 4) && (__GNUC_MINOR__ == 8)  && !__clang__)
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
#if ((__GNUC__ == 4) && (__GNUC_MINOR__ == 8)  && !__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Weffc++"
#endif
template<typename Derived>
struct is_matrix_expression : std::is_base_of
                <
                    Eigen::MatrixBase<typename std::decay<Derived>::type>,
                    typename std::decay<Derived>::type
                >
{
};
#if ((__GNUC__ == 4) && (__GNUC_MINOR__ == 8)  && !__clang__)
#pragma GCC diagnostic pop
#endif

/**
* \brief Checks whether the type is a complex type
*
* Provides the constant member \a value which is equal to \a true,
* if the type is a complex type, i.e. \a std::complex<T>
*/
// silence g++4.8.x warning about non-virtual destructor in inherited class
#if ((__GNUC__ == 4) && (__GNUC_MINOR__ == 8)  && !__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Weffc++"
#endif
template<typename T>
struct is_complex : std::false_type
{
};
#if ((__GNUC__ == 4) && (__GNUC_MINOR__ == 8)  && !__clang__)
#pragma GCC diagnostic pop
#endif

/**
* \brief \brief Checks whether the type is a complex number type,
* specialization for complex types
*/
// silence g++4.8.x warning about non-virtual destructor in inherited class
#if ((__GNUC__ == 4) && (__GNUC_MINOR__ == 8)  && !__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Weffc++"
#endif
template<typename T>
struct is_complex<std::complex<T>> : std::true_type
{
};
#if ((__GNUC__ == 4) && (__GNUC_MINOR__ == 8)  && !__clang__)
#pragma GCC diagnostic pop
#endif


} /* namespace qpp */

#endif /* TRAITS_H_ */

