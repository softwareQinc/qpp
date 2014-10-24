/* 
 * File:   qpp.h
 * Author: vlad
 *
 * Created on December 12, 2013, 10:42 PM
 */

// IMPORTANT: instantiation of global singletons RandomDevices, Gates and States
// Any additional singletons should be instantiated here
// Includes all necessary headers (except "matlab.h")
// ALWAYS include it in main.cpp

#ifndef QPP_H_
#define	QPP_H_

#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <cstring>
#include <exception>
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
#include <sstream>
#include <type_traits>
#include <utility>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/SVD>

// do not change the order of these!
#include "constants.h"
#include "types.h"
#include "classes/exception.h"
#include "classes/singleton.h"
#include "classes/states.h"
#include "classes/randevs.h"
#include "internal.h"
#include "functions.h"
#include "classes/gates.h"
#include "classes/stat.h"
#include "entropies.h"
#include "entanglement.h"

#include "channels.h"
#include "io.h"
#include "random.h"
#include "classes/qudit.h"
#include "classes/timer.h"

namespace qpp
{

/**
 * \brief qpp::RandomDevices Singleton
 *
 * Initializes the random devices, see the class \a qpp::RandomDevices
 */
RandomDevices& rdevs = RandomDevices::get_instance();

/**
 * \brief qpp::Gates const Singleton
 *
 * Initializes the gates, see the class \a qpp::Gates
 */
const Gates& gt = Gates::get_instance();

/**
 * \brief qpp::States const Singleton
 *
 * Initializes the states, see the class \a qpp::States
 */const States& st = States::get_instance();

} /* namespace qpp */

#endif	/* QPP_H_ */

