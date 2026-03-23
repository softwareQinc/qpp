/*
 * This file is part of pyqpp.
 *
 * Copyright (c) 2017 - 2026 softwareQ Inc. All rights reserved.
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
 * @file <pyqpp/classes/codes_bind.hpp>
 * @brief Bindings for <qpp/classes/codes.hpp>
 */

#ifndef PYQPP_CLASSES_CODES_BIND_HPP_
#define PYQPP_CLASSES_CODES_BIND_HPP_

#include "pyqpp/pyqpp_common.hpp"

inline void init_classes_codes(py::module_& m) {
    using namespace qpp;

    auto codes = m.def_submodule("codes");
    // Wrap the Type enum inside the codes submodule
    py::enum_<Codes::Type>(codes, "Type")
        .value("FIVE_QUBIT", Codes::Type::FIVE_QUBIT)
        .value("STEANE_SEVEN_QUBIT", Codes::Type::STEANE_SEVEN_QUBIT)
        .value("SHOR_NINE_QUBIT", Codes::Type::SHOR_NINE_QUBIT)
        .export_values();

    // Wrap the codeword static method
    // Note: Since codeword returns a 'ket', pybind11 handles the Eigen
    // conversion
    codes.def(
        "codeword",
        [](Codes::Type type, idx i) {
            return qpp::cmat(Codes::codeword(type, i));
        },
        "Returns the i-th codeword of the specified quantum code type",
        py::arg("type"), py::arg("i"));
}

#endif /* PYQPP_CLASSES_CODES_BIND_HPP_ */
