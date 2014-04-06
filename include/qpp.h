/* 
 * File:   qpp.h
 * Author: vlad
 *
 * Created on December 12, 2013, 10:42 PM
 */

// IMPORTANT: contains implementation/initialization of external variables
// ALWAYS include it in main.cpp and run _init() from main.cpp
#ifndef QPP_H_
#define	QPP_H_

#include <cstdlib>
#include "types.h"
#include "constants.h"
#include "gates.h"
#include "stat.h"
#include "functions.h"
#include "random.h"
#include "entropies.h"
#include "io.h"
#include "timer.h"
#include "exception.h"

// internal.h and matlab.h should not be included by default

namespace qpp
{

// make visible the gates declared as external in gates.h
namespace gt
{
types::cmat Id2, H, X, Y, Z, S, T;
types::cmat CNOTab, CNOTba, CZ, CS, SWAP;
types::cmat TOF, FRED;
types::cmat x0, x1, y0, y1, z0, z1;
types::cmat b00, b01, b10, b11;
}

// make the random_device generator visible
std::random_device stat::_rd;

// make the mt19937 generator visible and seed it
std::mt19937 stat::_rng(_rd());

// initialize the rest of the library, call this at the start of main.cpp
int _init()
{
	// seed the standard C random number generator (needed by Eigen)
	std::srand(static_cast<unsigned int>(stat::_rd()));

	// initialize the non-constant gates
	gt::_init_gates();

	return 0;
}

}

#endif	/* QPP_H_ */

