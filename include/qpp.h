/* 
 * File:   qpp.h
 * Author: vlad
 *
 * Created on December 12, 2013, 10:42 PM
 */

// IMPORTANT: instantiation of global singletons Gates and RandomDevices
// Any additional singletons should be instantiated here
// Includes all necessary headers
// ALWAYS include it in main.cpp
#ifndef QPP_H_
#define	QPP_H_

#include <algorithm>
#include <cmath>
#include <chrono>
#include <complex>
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <numeric>
#include <ostream>
#include <random>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/SVD>

// do not change the order of these!
#include "constants.h"
#include "types.h"
#include "classes/exception.h"
#include "internal.h"
#include "functions.h"
#include "classes/states.h"
#include "classes/gates.h"
#include "classes/randevs.h"
#include "classes/stat.h"
#include "channels.h"
#include "entanglement.h"
#include "entropies.h"
#include "io.h"
#include "random.h"
#include "classes/qudit.h"
#include "classes/timer.h"

namespace qpp
{

// initialize the random devices -- Singleton
RandomDevices& rdevs = RandomDevices::getInstance();

// initialize the gates -- Singleton
const Gates& gt = Gates::getInstance();

// initialize the states -- Singleton
const States& st = States::getInstance();

} /* namespace qpp */

#endif	/* QPP_H_ */

