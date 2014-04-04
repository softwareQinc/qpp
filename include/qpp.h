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
#include "util.h"
#include "constants.h"
#include "gates.h"
#include "stat.h"
#include "functional.h"
#include "random.h"
#include "entropy.h"
#include "io.h"
#include "timer.h"
#include "exception.h"

// internal.h and matlab.h should not be included by default

namespace qpp
{

// make visible the gates declared as external in gates.h
namespace gt
{
types::cmat H, Id2, X, Y, Z, S, T;
types::cmat CNOT, CP;
types::cmat TOF(8, 8);
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

