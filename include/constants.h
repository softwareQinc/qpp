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
// user defined literal for complex imaginary i (square root of -1)
constexpr std::complex<double> operator "" _i(unsigned long long int x)
{
	return
	{	0., static_cast<double>(x)};
}

constexpr std::complex<double> operator "" _i(long double x)
{
	return
	{	0., static_cast<double>(x)};
}

namespace ct
{

// used in disp(ln) function for setting to zero numbers for which |z|<chop
constexpr double chop = 1e-10;

// used to decide whether a number or expression in double precision
// is zero or not; example: if(std::abs(x)<eps) then x is 0
constexpr double eps = 1e-12;

// max number of systems, used to statically allocate arrays such as Cdims etc
constexpr std::size_t maxn = 64; // definitely cannot simulate more qubits :)

// math constants
//const std::complex<double> ii = { 0, 1 }; // Imaginary i (square root of -1)
constexpr double pi = 3.141592653589793238462643383279502884; // pi
constexpr double ee = 2.718281828459045235360287471352662497; // base of natural log

// D-th root of unity
std::complex<double> omega(std::size_t D) // D-th root of unity
{
	return exp(2.0 * pi * 1_i / static_cast<double>(D));
}

} /* namespace ct */
} /* namespace qpp */

#endif	/* CONSTANTS_H_ */

