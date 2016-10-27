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

// Collection of some useful type traits

namespace qpp
{

/**
* \brief Alias template that implements the proposal for void_t
*
* \see http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2014/n3911
*/
template<typename ...>
using to_void = void;

/**
* \brief Checks whether \a T is compatible with an STL-like iterable container
*
* Provides the member constant \a value which is equal to \a true,
* if \a T is compatible with an iterable container, i.e. provides at least
* \a begin() and \a end() member functions.
* Otherwise, \a value is equal to \a false.
*/
template<typename T, typename = void>
struct is_iterable : std::false_type
{
};

/**
* \brief Checks whether \a T is compatible with an STL-like iterable container,
* specialization for STL-like iterable containers
*/
template<typename T>
struct is_iterable<T,
        to_void<decltype(std::declval<T>().begin()),
                decltype(std::declval<T>().end()),
                typename T::value_type
        >> : std::true_type
{
};

/**
* \brief Checks whether the type is an Eigen matrix expression
*
* Provides the member constant \a value which is equal to \a true,
* if the type is an Eigen matrix expression of type \a
* Eigen::MatrixBase<Derived>.
* Otherwise, \a value is equal to \a false.
*/
template<typename... Derived>
struct is_matrix_expression : std::false_type
{
};

/**
* \brief Checks whether the type is an Eigen matrix expression,
* specialization for Eigen matrix expressions
*/
// thanks to @davidhigh http://stackoverflow.com/a/40293333/3093378
template<typename Derived>
struct is_matrix_expression<Derived>
        : std::is_base_of
                  <
                      Eigen::MatrixBase<typename std::decay<Derived>::type>,
                      typename std::decay<Derived>
                  >
{
};

/**
* \brief Checks whether the type is a complex type
*
* Provides the member constant \a value which is equal to \a true,
* if the type is a complex type, i.e. \a std::complex<T>
*/
template<typename T>
struct is_complex : std::false_type
{
};

/**
* \brief \brief Checks whether the type is a complex number type,
* specialization for complex types
*/
template<typename T>
struct is_complex<std::complex<T>> : std::true_type
{
};


} /* namespace qpp */

#endif /* TRAITS_H_ */

