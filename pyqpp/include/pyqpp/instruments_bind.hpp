/*
 * This file is part of pyqpp.
 *
 * Copyright (c) 2019 - 2024 softwareQ Inc. All rights reserved.
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

#ifndef PYQPP_INSTRUMENTS_BIND_HPP_
#define PYQPP_INSTRUMENTS_BIND_HPP_

#include "pyqpp/pyqpp_common.h"

/* Some free functions (non-exhaustive list) from instruments.hpp */
inline void init_instruments(py::module_& m) {
    using namespace qpp;

    m.def(
        "sample",
        [](const cmat& A, const std::vector<idx>& target,
           const std::vector<idx>& dims) {
            return qpp::sample(A, target, dims);
        },
        "Samples from a quantum state in the computational basis (Z-basis)",
        py::arg("A"), py::arg("target"), py::arg("dims"));
    m.def(
        "sample",
        [](const cmat& A, const std::vector<idx>& target, idx d = 2) {
            return qpp::sample(A, target, d);
        },
        "Samples from a quantum state in the computational basis (Z-basis)",
        py::arg("A"), py::arg("target"), py::arg("d") = 2);
    m.def(
        "sample",
        [](idx num_samples, const cmat& A, const std::vector<idx>& target,
           const std::vector<idx>& dims) {
            std::map<std::string, idx> result;
            auto stats = qpp::sample(num_samples, A, target, dims);
            for (auto&& elem : stats) {
                std::stringstream ss;
                ss << qpp::disp(
                    elem.first,
                    IOManipContainerOpts{}.set_sep("").set_left("").set_right(
                        ""));
                result[ss.str()] = elem.second;
            }
            return result;
        },
        "Samples repeatedly from a quantum state in the computational basis "
        "(Z-basis)",
        py::arg("num_samples"), py::arg("A"), py::arg("target"),
        py::arg("dims"));
    m.def(
        "sample",
        [](idx num_samples, const cmat& A, const std::vector<idx>& target,
           idx d = 2) {
            std::map<std::string, idx> result;
            auto stats = qpp::sample(num_samples, A, target, d);
            for (auto&& elem : stats) {
                std::stringstream ss;
                ss << qpp::disp(
                    elem.first,
                    IOManipContainerOpts{}.set_sep("").set_left("").set_right(
                        ""));
                result[ss.str()] = elem.second;
            }
            return result;
        },
        "Samples repeatedly from a quantum state in the computational basis "
        "(Z-basis)",
        py::arg("num_samples"), py::arg("A"), py::arg("target"),
        py::arg("d") = 2);
}

#endif /* PYQPP_INSTRUMENTS_BIND_HPP_ */
