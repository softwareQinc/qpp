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
    auto py_dirac = py::class_<py_dirac_t>(m, "dirac_t");
    py_dirac.def(py::self == py::self);
    py_dirac.def(py::self != py::self);
    py_dirac.def("__copy__",
                 [](const py_dirac_t& self) { return py_dirac_t(self); });
    py_dirac.def("__deepcopy__", [](const py_dirac_t& self, py::dict) {
        return py_dirac_t(self);
    });
    py_dirac.def("__repr__", [](const py_dirac_t& self) {
        std::ostringstream oss;
        oss << disp(self);
        return oss.str();
    });

    using py_proxy_to_engine_dits_t = internal::LabelledVectorProxy<idx, false>;
    using py_const_proxy_to_engine_dits_t =
        internal::LabelledVectorProxy<idx, true>;

    auto py_proxy_to_engine_dits =
        py::class_<py_proxy_to_engine_dits_t>(m, "proxy_to_engine_dits");
    py_proxy_to_engine_dits.def(
        "__getitem__",
        [](const py_proxy_to_engine_dits_t& self, idx i) { return self[i]; });
    py_proxy_to_engine_dits.def(
        "__setitem__",
        [](py_proxy_to_engine_dits_t& self, idx i) { return self[i]; });

    auto py_const_proxy_to_engine_dits =
        py::class_<py_const_proxy_to_engine_dits_t>(
            m, "const_proxy_to_engine_dits");
    py_const_proxy_to_engine_dits.def(
        "__getitem__",
        [](py_const_proxy_to_engine_dits_t& self, idx i) { return self[i]; });
}

#endif /* PYQPP_TYPES_BIND_HPP_ */
