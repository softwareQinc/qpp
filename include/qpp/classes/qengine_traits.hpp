/*
 * This file is part of Quantum++.
 *
 * Copyright (c) 2017 - 2024 softwareQ Inc. All rights reserved.
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

/**
 * \file qpp/classes/qengine_traits.hpp
 * \brief Quantum engine traits (run-time)
 */

#ifndef QPP_CLASSES_QENGINE_TRAITS_HPP_
#define QPP_CLASSES_QENGINE_TRAITS_HPP_

#include <string>

namespace qpp {
/**
 * \class qpp::IQEngineTraits
 * \brief Traits for quantum engines, resolved at run-time
 * \note All engines must implement this trait (i.e., inherit from it)
 */
struct IQEngineTraits {
    /**
     * \brief Engine name
     *
     * \return Engine name
     */
    virtual std::string traits_get_name() const = 0;

    /**
     * \brief Determines if the engine is noisy
     *
     * \return True if the engine simulates noisy execution, false if not
     */

    virtual bool traits_is_noisy() const = 0;
    /**
     * \brief Determines if the engine operates on pure states
     *
     * \return True if the engine operates on pure states, false otherwise
     */
    virtual bool traits_is_pure() const = 0;

    /**
     * \brief Default virtual destructor
     */
    virtual ~IQEngineTraits() = default;
};
} /* namespace qpp */

#endif /* QPP_CLASSES_QENGINE_TRAITS_HPP_ */
