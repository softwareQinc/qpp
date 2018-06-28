/*
 * This file is part of Quantum++.
 *
 * MIT License
 *
 * Copyright (c) 2013 - 2018 Vlad Gheorghiu (vgheorgh@gmail.com)
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
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
#define QPP_UNUSED_ __attribute__((unused))
#else
#define QPP_UNUSED_
#endif

// standard C++ library headers
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

// Eigen headers
#include <Eigen/Dense>
#include <Eigen/SVD>

// Quantum++ headers

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
#include "classes/reversible.h"

/**
* \namespace qpp
* \brief Quantum++ main namespace
*/
namespace qpp {
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
static RandomDevices& rdevs QPP_UNUSED_ = RandomDevices::get_instance();
#else
thread_local static RandomDevices& rdevs QPP_UNUSED_ =
    RandomDevices::get_thread_local_instance();
#endif // NO_THREAD_LOCAL_

} /* namespace qpp */

#endif /* QPP_H_ */
