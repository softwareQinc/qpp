/*
 * This file is part of Quantum++.
 *
 * Copyright (c) 2017 - 2025 softwareQ Inc. All rights reserved.
 *
 * MIT License
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
// Includes all necessary headers (except "matlab.hpp")
// ALWAYS include it in main.cpp

/**
 * \file qpp/qpp.hpp
 * \brief Quantum++ main header file, includes all other required headers
 */

#ifndef QPP_QPP_HPP_
#define QPP_QPP_HPP_

// Ignore warnings for unknown C++17 attributes (we use such "custom"
// attributes internally, the compiler is supposed to ignore them according to
// the C++17 standard)

// Intel
#if defined(__INTEL_COMPILER)
#pragma warning(disable : 3924)

// Clang
#elif defined(__clang__)
#pragma clang diagnostic ignored "-Wunknown-attributes"

// GCC
#elif defined(__GNUC__) && !defined(__clang__) && !defined(__INTEL_COMPILER)
#pragma GCC diagnostic ignored "-Wattributes"

// MSVC
#elif defined(_MSC_VER) && !defined(__INTEL_COMPILER)
#pragma warning(disable : 5030)

#endif

// SunOS Solaris/OpenIndiana issues
#if defined(__sun)
// CTRL defined as a macro
#ifdef CTRL
#undef CTRL
#endif // CTRL
#endif // __sun

// Quantum++ library headers

#include "qpp/constants.hpp"
#include "qpp/entanglement.hpp"
#include "qpp/entropies.hpp"
#include "qpp/functions.hpp"
#include "qpp/input_output.hpp"
#include "qpp/instruments.hpp"
#include "qpp/number_theory.hpp"
#include "qpp/operations.hpp"
#include "qpp/options.hpp"
#include "qpp/random.hpp"
#include "qpp/statistics.hpp"
#include "qpp/traits.hpp"
#include "qpp/types.hpp"

#include "qpp/classes/codes.hpp"
#include "qpp/classes/exception.hpp"
#include "qpp/classes/gates.hpp"
#include "qpp/classes/idisplay.hpp"
#include "qpp/classes/ijson.hpp"
#include "qpp/classes/init.hpp"
#include "qpp/classes/layouts.hpp"
#include "qpp/classes/noise.hpp"
#include "qpp/classes/random_devices.hpp"
#include "qpp/classes/reversible.hpp"
#include "qpp/classes/states.hpp"
#include "qpp/classes/timer.hpp"

#include "qpp/classes/qbase_engine.hpp"
#include "qpp/classes/qdummy_engine.hpp"
#include "qpp/classes/qengine.hpp"
#include "qpp/classes/qengine_traits.hpp"
#include "qpp/classes/qnoisy_engine.hpp"

#include "qpp/classes/qcircuit.hpp"
#include "qpp/classes/qcircuit_traits.hpp"

#include "qpp/qasm/qasm.hpp"

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
static const Init& init [[maybe_unused]] = Init::get_no_thread_local_instance();

/**
 * \brief qpp::Codes const Singleton
 *
 * Initializes the codes, see the class qpp::Codes
 */
static const Codes& codes [[maybe_unused]] =
    Codes::get_no_thread_local_instance();

/**
 * \brief qpp::Gates const Singleton
 *
 * Initializes the gates, see the class qpp::Gates
 */
static const Gates& gt [[maybe_unused]] = Gates::get_no_thread_local_instance();

/**
 * \brief qpp::States const Singleton
 *
 * Initializes the states, see the class qpp::States
 */
static const States& st [[maybe_unused]] =
    States::get_no_thread_local_instance();

/**
 * \brief qpp::RandomDevices Singleton
 *
 * Initializes the random devices, see the class qpp::RandomDevices
 *
 * \note If the compiler supports thread_local, has thread_local storage
 * duration, due to mutability of its public member std::mt19937 and
 * possible data races
 */

#ifndef NO_THREAD_LOCAL_
thread_local
#endif
    static RandomDevices& rdevs [[maybe_unused]] =
        RandomDevices::get_instance();

/**
 * \namespace qpp::obsolete
 * \brief Obsolete/deprecated code, may be removed without notice in future
 * releases
 */
namespace obsolete {} /* namespace obsolete */

} /* namespace qpp */

#endif /* QPP_QPP_HPP_ */
