/*
 * This file is part of pyqpp.
 *
 * Copyright (c) 2019 - 2023 softwareQ Inc. All rights reserved.
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

/* Types from types.hpp */
inline void init_types(py::module_& m) {
    using namespace qpp;

    // supports only complex
    using py_io_braket = io_braket<cplx>;

    /* qpp::io_braket */
    auto pyio_braket =
        py::class_<py_io_braket>(m, "io_braket")
            .def(py::self == py::self)
            .def(py::self != py::self)
            .def("__copy__",
                 [](const py_io_braket& self) { return io_braket(self); })
            .def("__deepcopy__", [](const py_io_braket& self,
                                    py::dict) { return io_braket(self); })
            .def("__repr__", [](const py_io_braket& iob) {
                std::ostringstream oss;
                oss << disp(iob, false, "\n", " * ");
                return oss.str();
            });
}

#endif /* PYQPP_TYPES_BIND_HPP_ */
