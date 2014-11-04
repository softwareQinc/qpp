/* 
 * Quantum++
 *
 * Copyright (c) 2013 - 2014 Vlad Gheorghiu (vgheorgh@gmail.com)
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

// do not change the order in this group, inter-dependencies
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

// the ones below can be in any order, no inter-dependencies
#include "instruments.h"
#include "io.h"
#include "operations.h"
#include "random.h"
#include "classes/timer.h"

// experimental/testing
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

