/*
 * Quantum++
 *
 * Copyright (c) 2013 - 2015 Vlad Gheorghiu (vgheorgh@gmail.com)
 *
 * This file is part of Quantum++.
 *
 * Quantum++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Quantum++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Quantum++.  If not, see <http://www.gnu.org/licenses/>.
 */

// IMPORTANT: instantiation of global singletons
// Init, Codes, Gates, States and RandomDevices
//
// Any additional singletons should be instantiated here
// Includes all necessary headers (except "matlab.h")
// ALWAYS include it in main.cpp

/**
* \file qpp.h
* \brief Quantum++ main header file, includes all other necessary headers
*/

#ifndef QPP_H_
#define QPP_H_

// silence bogus warning -Wunused-variable for Singletons
#if (__GNUC__)
#pragma GCC diagnostic ignored "-Wunused-variable"
#endif

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
#include <initializer_list>
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

// do not change the order in this group, inter-dependencies
#include "types.h"
#include "classes/exception.h"
#include "constants.h"
#include "traits.h"
#include "classes/idisplay.h"
#include "internal/util.h"
#include "internal/classes/iomanip.h"
#include "input_output.h"

// do not change the order in this group, inter-dependencies
#include "internal/classes/singleton.h"
#include "classes/init.h"
#include "functions.h"
#include "classes/codes.h"
#include "classes/gates.h"
#include "classes/states.h"
#include "classes/random_devices.h"

// do not change the order in this group, inter-dependencies
#include "statistics.h"
#include "operations.h"
#include "entropies.h"
#include "entanglement.h"

// the ones below can be in any order, no inter-dependencies
#include "random.h"
#include "classes/timer.h"
#include "instruments.h"
#include "number_theory.h"


/**
* \namespace qpp
* \brief Quantum++ main namespace
*/
namespace qpp
{

/**
* \brief qpp::Init const Singleton
*
* Additional initializations/cleanups, see the class qpp::Init
*/
static const Init& init = Init::get_instance();

/**
* \brief qpp::Codes const Singleton
*
* Initializes the codes, see the class qpp::Codes
*/
static const Codes& codes = Codes::get_instance();

/**
* \brief qpp::Gates const Singleton
*
* Initializes the gates, see the class qpp::Gates
*/
static const Gates& gt = Gates::get_instance();

/**
* \brief qpp::States const Singleton
*
* Initializes the states, see the class qpp::States
*/
static const States& st = States::get_instance();

/**
* \brief qpp::RandomDevices Singleton
*
* Initializes the random devices, see the class qpp::RandomDevices
*
* \note Has thread storage duration, due to mutability of its public member
* std::mt19937 and possible data races
*/
#ifndef _NO_THREAD_LOCAL
thread_local static RandomDevices& rdevs =
        RandomDevices::get_thread_local_instance();
#else
static RandomDevices& rdevs =
    RandomDevices::get_instance();
#endif // _NO_THREAD_LOCAL

} /* namespace qpp */

#endif  /* QPP_H_ */
