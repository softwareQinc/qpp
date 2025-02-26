/*
 * This file is part of pyqpp.
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

#ifndef PYQPP_TYPES_BIND_HPP_
#define PYQPP_TYPES_BIND_HPP_

#include "pyqpp/pyqpp_common.hpp"

/* Types from types.hpp */
inline void init_types(py::module_& m) {
    using namespace qpp;

    // supports only complex
    using py_dirac_t = dirac_t<cplx>;

    /* qpp::dirac_t */
    auto pydirac_t = py::class_<py_dirac_t>(m, "dirac_t");
    pydirac_t.def(py::self == py::self);
    pydirac_t.def(py::self != py::self);
    pydirac_t.def("__copy__",
                  [](const py_dirac_t& self) { return py_dirac_t(self); });
    pydirac_t.def("__deepcopy__", [](const py_dirac_t& self, py::dict) {
        return py_dirac_t(self);
    });
    pydirac_t.def("__repr__", [](const py_dirac_t& self) {
        std::ostringstream oss;
        oss << disp(self);
        return oss.str();
    });
}

#endif /* PYQPP_TYPES_BIND_HPP_ */
