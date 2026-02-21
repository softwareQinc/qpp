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
 * @file <pyqpp/entanglement_bind.hpp>
 * @file <pyqpp/entanglement_bind.hpp>
 * @brief Bindings for <qpp/entanglement.hpp>
 * @brief Bindings for <qpp/entanglement.hpp>
 */

#ifndef PYQPP_ENTANGLEMENT_BIND_HPP_
#define PYQPP_ENTANGLEMENT_BIND_HPP_

#include "pyqpp/pyqpp_common.hpp"

inline void init_entanglement(py::module_& m) {
    using namespace qpp;

    // --- Schmidt Decomposition and Coefficients ---

    m.def(
        "schmidtcoeffs",
        [](const cmat& A, const std::vector<idx>& dims) {
            return qpp::schmidtcoeffs(A, dims);
        },
        "Schmidt coefficients of the bi-partite pure state A.", py::arg("A"),
        py::arg("dims"));

    m.def(
        "schmidtcoeffs",
        [](const cmat& A, idx d) { return qpp::schmidtcoeffs(A, d); },
        "Schmidt coefficients of the bi-partite pure state A (uniform "
        "dimensions).",
        py::arg("A"), py::arg("d") = 2);

    m.def(
        "schmidtA",
        [](const cmat& A, const std::vector<idx>& dims) {
            return qpp::schmidtA(A, dims);
        },
        "Schmidt basis on Alice side.", py::arg("A"), py::arg("dims"));

    m.def(
        "schmidtB",
        [](const cmat& A, const std::vector<idx>& dims) {
            return qpp::schmidtB(A, dims);
        },
        "Schmidt basis on Bob side.", py::arg("A"), py::arg("dims"));

    m.def(
        "schmidtprobs",
        [](const cmat& A, const std::vector<idx>& dims) {
            return qpp::schmidtprobs(A, dims);
        },
        "Schmidt probabilities (squares of coefficients) of the bi-partite "
        "pure state A.",
        py::arg("A"), py::arg("dims"));

    m.def(
        "schmidt",
        [](const cmat& A, const std::vector<idx>& dims) {
            return qpp::schmidt(A, dims);
        },
        "Full Schmidt decomposition: returns tuple(U, V, coeffs, probs).",
        py::arg("A"), py::arg("dims"));

    // --- Entanglement Measures ---

    m.def(
        "entanglement",
        [](const cmat& A, const std::vector<idx>& dims) {
            return qpp::entanglement(A, dims);
        },
        "von-Neumann entanglement entropy of the bi-partite pure state A.",
        py::arg("A"), py::arg("dims"));

    m.def(
        "gconcurrence", [](const cmat& A) { return qpp::gconcurrence(A); },
        "G-concurrence of the bi-partite pure state A.", py::arg("A"));

    m.def(
        "negativity",
        [](const cmat& A, const std::vector<idx>& dims) {
            return qpp::negativity(A, dims);
        },
        "Negativity of the bi-partite mixed state A.", py::arg("A"),
        py::arg("dims"));

    m.def(
        "lognegativity",
        [](const cmat& A, const std::vector<idx>& dims) {
            return qpp::lognegativity(A, dims);
        },
        "Logarithmic negativity of the bi-partite mixed state A.", py::arg("A"),
        py::arg("dims"));

    m.def(
        "concurrence", [](const cmat& A) { return qpp::concurrence(A); },
        "Wootters concurrence of the bi-partite qubit mixed state A.",
        py::arg("A"));
}
#endif /* PYQPP_ENTANGLEMENT_BIND_HPP_ */
