/* 
 * Quantum++
 *
 * Copyright (c) 2013 - 2015 Vlad Gheorghiu (vgheorgh@gmail.com)
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
* \file constants.h
* \brief Constants
*/

#ifndef CONSTANTS_H_
#define CONSTANTS_H_

// constants

namespace qpp
{

/**
* \brief User-defined literal for complex \f$i=\sqrt{-1}\f$ (integer overload)
*
* Example: \code auto z = 4_i; // type of z is std::complex<double> \endcode
*/
inline constexpr cplx operator""_i(unsigned long long int x) noexcept
{
    return {0., static_cast<double>(x)};
}

/**
* \brief User-defined literal for complex \f$i=\sqrt{-1}\f$ (real overload)
*
* Example: \code auto z = 4.5_i; // type of z is std::complex<double> \endcode
*/
inline constexpr cplx operator""_i(long double x) noexcept
{
    return {0., static_cast<double>(x)};
}

/**
* \brief  Used in qpp::disp() for setting to zero  numbers
* that have their absolute value smaller than qpp::chop
*/
constexpr double chop = 1e-10;

/**
* \brief Used to decide whether a number or expression in double precision
* is zero or not
*
* Example: \code if(std::abs(x) < qpp::eps) // x is zero \endcode
*/
constexpr double eps = 1e-12;

/**
* \brief Maximum number of allowed qu(d)its (subsystems)
*
* Used internally to allocate arrays on the stack (for speed reasons)
*/
constexpr idx maxn = 64; // maximum number of qudits a state may have

/**
* \brief \f$ \pi \f$
*/
constexpr double pi = 3.141592653589793238462643383279502884;
/**
* \brief Base of natural logarithm, \f$e\f$
*/
constexpr double ee = 2.718281828459045235360287471352662497;

/**
* \brief Used to denote infinity in double precision
*/
constexpr double infty = std::numeric_limits<double>::infinity();

/**
* \brief D-th root of unity
*
* \param D Non-negative integer
* \return D-th root of unity \f$\exp(2\pi i/D)\f$
*/
inline cplx omega(idx D)
{
    if (D == 0)
        throw Exception::Type::OUT_OF_RANGE;
    return exp(2.0 * pi * 1_i / static_cast<double>(D));
}

} /* namespace qpp */

#endif	/* CONSTANTS_H_ */
