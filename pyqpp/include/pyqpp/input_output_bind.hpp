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

/**
 * \file <pyqpp/input_output_bind.hpp>
 * \brief Bindings for <qpp/input_output.hpp>
 */

#ifndef PYQPP_INPUT_OUTPUT_BIND_HPP_
#define PYQPP_INPUT_OUTPUT_BIND_HPP_

#include "pyqpp/pyqpp_common.hpp"

inline void init_input_output(py::module_& m) {
    using namespace qpp;

    // --- qpp::save() bindings ---
    m.def(
        "save",
        [](const cmat& A, const std::string& filename) {
            std::ofstream fout(filename);
            if (!fout.is_open()) {
                throw std::runtime_error("Could not open file for writing: " +
                                         filename);
            }
            qpp::save(A, fout);
        },
        py::arg("A"), py::arg("filename"),
        "Saves a complex matrix to a text file");

    m.def(
        "save",
        [](const rmat& A, const std::string& filename) {
            std::ofstream fout(filename);
            if (!fout.is_open()) {
                throw std::runtime_error("Could not open file for writing: " +
                                         filename);
            }
            qpp::save(A, fout);
        },
        py::arg("A"), py::arg("filename"),
        "Saves a real matrix to a text file");

    // --- qpp::load() bindings ---
    m.def(
        "load_cmat",
        [](const std::string& filename) {
            std::ifstream fin(filename);
            if (!fin.is_open()) {
                throw std::runtime_error("Could not open file for reading: " +
                                         filename);
            }
            return qpp::load<cmat>(fin);
        },
        py::arg("filename"), "Loads a complex matrix from a text file");

    m.def(
        "load_rmat",
        [](const std::string& filename) {
            std::ifstream fin(filename);
            if (!fin.is_open()) {
                throw std::runtime_error("Could not open file for reading: " +
                                         filename);
            }
            return qpp::load<rmat>(fin);
        },
        py::arg("filename"), "Loads a real matrix from a text file");
}

#endif /* PYQPP_INPUT_OUTPUT_BIND_HPP_ */
