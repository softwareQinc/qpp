/*
 * Quantum++
 *
 * Copyright (c) 2013 - 2016 Vlad Gheorghiu (vgheorgh@gmail.com)
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

// silence bogus warning -Wunused-variable for singletons
#if (__GNUC__ && !__clang__)
#define QPP_UNUSED_  __attribute__ ((unused))
#else
#define QPP_UNUSED_
#endif

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <fstream>
#include <functional>
#include <initializer_list>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <memory>
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

// pre-processor macros, make them visible to the whole library
#include "macros.h"

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
static const Init& init QPP_UNUSED_ = Init::get_instance();

/**
* \brief qpp::Codes const Singleton
*
* Initializes the codes, see the class qpp::Codes
*/
static const Codes& codes QPP_UNUSED_ = Codes::get_instance();

/**
* \brief qpp::Gates const Singleton
*
* Initializes the gates, see the class qpp::Gates
*/
static const Gates& gt QPP_UNUSED_ = Gates::get_instance();

/**
* \brief qpp::States const Singleton
*
* Initializes the states, see the class qpp::States
*/
static const States& st QPP_UNUSED_ = States::get_instance();

/**
* \brief qpp::RandomDevices Singleton
*
* Initializes the random devices, see the class qpp::RandomDevices
*
* \note Has thread storage duration, due to mutability of its public member
* std::mt19937 and possible data races
*/
#ifdef NO_THREAD_LOCAL_
static RandomDevices& rdevs QPP_UNUSED_= RandomDevices::get_instance();
#else
thread_local static RandomDevices& rdevs QPP_UNUSED_ =
        RandomDevices::get_thread_local_instance();
#endif // NO_THREAD_LOCAL_

} /* namespace qpp */

#endif  /* QPP_H_ */
