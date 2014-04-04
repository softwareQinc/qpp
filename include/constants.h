/* 
 * File:   constants.h
 * Author: vlad
 *
 * Created on December 12, 2013, 10:41 PM
 */

#ifndef CONSTANTS_H_
#define	CONSTANTS_H_

#include "types.h"

// constants

namespace qpp
{
namespace ct
{

// used for setting to zero everything that is smaller in absolute magnitude
const double chop = 1e-10;

// math constants
const types::cplx ii = { 0, 1 }; // Imaginary i (square root of -1)
const double pi = 3.141592653589793238462643383279502884; // pi
const double ee = 2.718281828459045235360287471352662497; // base of natural log

// D-th root of unity
inline types::cplx omega(size_t D) // D-th root of unity
{
	return exp(2.0 * pi * ii / static_cast<double>(D));
}

}
}

#endif	/* CONSTANTS_H_ */

