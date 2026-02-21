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
 * @file <pyqpp/entropies_bind.hpp>
 * @brief Bindings for <qpp/entropies.hpp>
 */

#ifndef PYQPP_ENTROPIES_BIND_HPP_
#define PYQPP_ENTROPIES_BIND_HPP_

#include "pyqpp/pyqpp_common.hpp"

inline void init_entropies(py::module_& m) {
    using namespace qpp;

    // --- Entropy Bindings ---

    // von-Neumann entropy (Templated for MatrixBase)
    m.def(
        "entropy", [](const cmat& A) { return qpp::entropy(A); },
        "Computes the von-Neumann entropy of a density matrix (base 2).",
        py::arg("A"));

    // Shannon entropy (Non-templated, but kept as lambda for consistency)
    m.def(
        "entropy",
        [](const std::vector<realT>& prob) { return qpp::entropy(prob); },
        "Computes the Shannon entropy of a probability distribution (base 2).",
        py::arg("prob"));

    // Renyi entropy (Templated for MatrixBase)
    m.def(
        "renyi",
        [](const cmat& A, realT alpha) { return qpp::renyi(A, alpha); },
        "Computes the Renyi-alpha entropy of a density matrix.", py::arg("A"),
        py::arg("alpha"));

    // Renyi entropy (Vector version)
    m.def(
        "renyi",
        [](const std::vector<realT>& prob, realT alpha) {
            return qpp::renyi(prob, alpha);
        },
        "Computes the Renyi-alpha entropy of a probability distribution.",
        py::arg("prob"), py::arg("alpha"));

    // Tsallis entropy (Templated for MatrixBase)
    m.def(
        "tsallis", [](const cmat& A, realT q) { return qpp::tsallis(A, q); },
        "Computes the Tsallis-q entropy of a density matrix.", py::arg("A"),
        py::arg("q"));

    // Tsallis entropy (Vector version)
    m.def(
        "tsallis",
        [](const std::vector<realT>& prob, realT q) {
            return qpp::tsallis(prob, q);
        },
        "Computes the Tsallis-q entropy of a probability distribution.",
        py::arg("prob"), py::arg("q"));

    // --- Mutual Information Bindings ---

    // Quantum Mutual Information (Templated, varied dimensions)
    m.def(
        "qmutualinfo",
        [](const cmat& A, const std::vector<idx>& subsysA,
           const std::vector<idx>& subsysB, const std::vector<idx>& dims) {
            return qpp::qmutualinfo(A, subsysA, subsysB, dims);
        },
        "Computes the quantum mutual information between two subsystems A and "
        "B with given dimensions.",
        py::arg("A"), py::arg("subsysA"), py::arg("subsysB"), py::arg("dims"));

    // Quantum Mutual Information (Templated, uniform dimensions)
    m.def(
        "qmutualinfo",
        [](const cmat& A, const std::vector<idx>& subsysA,
           const std::vector<idx>& subsysB,
           idx d) { return qpp::qmutualinfo(A, subsysA, subsysB, d); },
        "Computes the quantum mutual information between two subsystems A and "
        "B assuming uniform dimensions d.",
        py::arg("A"), py::arg("subsysA"), py::arg("subsysB"), py::arg("d") = 2);
}

#endif /* PYQPP_ENTROPIES_BIND_HPP_ */
