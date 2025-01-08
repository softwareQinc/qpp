/*
 * This file is part of pyqpp.
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

#ifndef PYQPP_QASM_QASM_BIND_HPP_
#define PYQPP_QASM_QASM_BIND_HPP_

#include "pyqpp/pyqpp_common.hpp"

/* OpenQASM interfacing */
inline void init_qasm_qasm(py::module_& m) {
    using namespace qpp;

    auto py_qasm = m.def_submodule("qasm");
    py_qasm.def("read_from_file", &qpp::qasm::read_from_file,
                "Get QCircuit representation of OpenQASM circuit");

    py_qasm.def("read_from_string", &qpp::qasm::read_from_string,
                "Get QCircuit representation of OpenQASM circuit");
}

#endif /* PYQPP_QASM_QASM_BIND_HPP_ */
