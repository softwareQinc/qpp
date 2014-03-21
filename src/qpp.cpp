/* 
 * File:   qpp.cpp
 * Author: vlad
 *
 * Created on December 12, 2013, 10:43 PM
 */

#include <cstdlib>
#include <ctime>
#include "qpp.h"
#include "gates.h"
#include "stat.h"

namespace qpp
{

std::random_device stat::_rd; // make the random_device generator visible
std::mt19937 stat::_rng(_rd()); // make the mt19937 generator visible

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
