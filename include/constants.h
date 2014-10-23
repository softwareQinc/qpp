/* 
 * File:   constants.h
 * Author: vlad
 *
 * Created on December 12, 2013, 10:41 PM
 */

#ifndef CONSTANTS_H_
#define	CONSTANTS_H_

// constants

namespace qpp
{
/**
 * @brief User-defined literal for complex \f$i=\sqrt{-1}\f$ (integer overload)
 *
 * Example: @code auto z = 4_i; // type of z is std::complex<double> @endcode
 */
constexpr std::complex<double> operator "" _i(unsigned long long int x)
{
	return
	{	0., static_cast<double>(x)};
}

/**
 * @brief User-defined literal for complex \f$i=\sqrt{-1}\f$ (real overload)
 *
 * Example: @code auto z = 4.5_i; // type of z is std::complex<double> @endcode
 */
constexpr std::complex<double> operator "" _i(long double x)
{
	return
	{	0., static_cast<double>(x)};
}

namespace ct
{

/**
 * @brief  Used in @a qpp::disp() and @a qpp::displn() for setting to zero
 * numbers that have their absolute value smaller than @a qpp::ct::chop
 *
 *
 */
constexpr double chop = 1e-10;

/**
 * @brief Used to decide whether a number or expression in double precision
 * is zero or not
 *
 * Example: @code if(std::abs(x) < qpp::ct::eps) // x is zero @endcode
 */
constexpr double eps = 1e-12;

/**
 * @brief Maximum number of qubits
 *
 * Used internally to statically allocate arrays (for speed reasons)
 */
constexpr std::size_t maxn = 64; // definitely cannot simulate more qubits :)

/**
 * @brief \f$ \pi \f$
 */
constexpr double pi = 3.141592653589793238462643383279502884;
/**
 * @brief Base of natural logarithm, \f$e\f$
 */
constexpr double ee = 2.718281828459045235360287471352662497;

/**
 * @brief D-th root of unity
 *
 * @param D Non-negative integer
 * @return D-th root of unity \f$\exp(2\pi i/D)\f$
 */
std::complex<double> omega(std::size_t D)
{
	return exp(2.0 * pi * 1_i / static_cast<double>(D));
}

} /* namespace ct */
} /* namespace qpp */

#endif	/* CONSTANTS_H_ */

