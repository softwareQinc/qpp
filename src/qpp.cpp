/* 
 * File:   qpp.cpp
 * Author: vlad
 *
 * Created on December 12, 2013, 10:43 PM
 */

#include <cstdlib>
#include <ctime>
#include "qpp.h"

namespace qpp
{

// Random number initialization
std::random_device stat::_rd; // make the random_device generator visible
std::mt19937 stat::_rng(_rd()); // make the mt19937 generator visible

// Gate initialization
namespace gt
{
// various matrices (declared extern in "gates.h") ; make the gates visible
types::cmat H, Id2, X, Y, Z, S, T;
types::cmat CNOT, CP;
types::cmat TOF(8, 8);
}


// initialize the library
int _init()
{
	// initialize the gates
	gt::_init_gates();

// seed the standard random number generator (needed by Eigen)
	std::srand(static_cast<unsigned int>(stat::_rd()));

	return 0;
}

}
