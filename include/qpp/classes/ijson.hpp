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
 * \file classes/ijson.hpp
 * \brief Basic JSON serialization
 */

#ifndef QPP_CLASSES_IJSON_HPP_
#define QPP_CLASSES_IJSON_HPP_

#include <string>

namespace qpp {
/**
 * \class qpp::IJSON
 * \brief Abstract class (interface) that mandates the definition of
 * very basic JSON serialization support
 */
class IJSON {
  public:
    /**
     * \brief Default virtual destructor
     */
    virtual ~IJSON() = default;

    /**
     * \brief JSON representation of the derived instance, must be overridden by
     * all derived classes
     *
     * \param enclosed_in_curly_brackets If true, encloses the result in curly
     * brackets
     */
    virtual std::string
    to_JSON(bool enclosed_in_curly_brackets = true) const = 0;

}; /* class IJSON */

} /* namespace qpp */

#endif /* QPP_CLASSES_IJSON_HPP_ */
