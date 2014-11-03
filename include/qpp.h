/* 
 * File:   qpp.h
 * Author: vlad
 *
 * Created on December 12, 2013, 10:42 PM
 */

// IMPORTANT: instantiation of global singletons
// RandomDevices, Gates, States and Init
// Any additional singletons should be instantiated here
// Includes all necessary headers (except "matlab.h")
// ALWAYS include it in main.cpp
/**
 * \mainpage quantum++ - A C++11 quantum computing library
 * \version 0.1
 * \author Vlad Gheorghiu, vgheorgh@gmail.com
 * \date \f$\mathrm{\today}\f$
 *
 * This is the main page of the documentation. More coming soon.
 */

#ifndef INCLUDE_QPP_H_
#define	INCLUDE_QPP_H_

#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <exception>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <numeric>
#include <ostream>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
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
#include "internal/functions.h"
#include "classes/init.h"
#include "functions.h"
#include "classes/gates.h"
#include "entropies.h"
#include "entanglement.h"

#include "operations.h"
#include "io.h"
#include "random.h"
#include "classes/timer.h"

#include "experimental/test.h"
#include "experimental/classes/qudit.h"

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
 */
const States& st = States::get_instance();

/**
 * \brief qpp::Init const Singleton
 *
 * Additional initializations/cleanups
 */
const Init& init = Init::get_instance();

} /* namespace qpp */

#endif	/* INCLUDE_QPP_H_ */

