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
 * \file <pyqpp/random_bind.hpp>
 * \brief Bindings for <qpp/random.hpp>
 */

#ifndef PYQPP_RANDOM_BIND_HPP_
#define PYQPP_RANDOM_BIND_HPP_

#include "pyqpp/pyqpp_common.hpp"

inline void init_random(py::module_& m) {
    using namespace qpp;

    // --- Arithmetic Randomization ---

    m.def(
        "rand", [](realT a, realT b) { return qpp::rand(a, b); },
        "Generates a random real number in [a, b).", py::arg("a"),
        py::arg("b"));

    m.def(
        "randn",
        [](realT mean, realT sigma) { return qpp::randn(mean, sigma); },
        "Generates a random real number normally distributed in N(mean, "
        "sigma).",
        py::arg("mean") = 0, py::arg("sigma") = 1);

    m.def("randidx", &qpp::randidx,
          "Generates a random index uniformly distributed in [a, b].",
          py::arg("a") = std::numeric_limits<idx>::min(),
          py::arg("b") = std::numeric_limits<idx>::max());

    // --- Matrix Randomization ---

    m.def(
        "rand",
        [](idx rows, idx cols, realT a, realT b) {
            return qpp::rand<cmat>(rows, cols, a, b);
        },
        "Generates a random complex matrix with entries uniformly distributed "
        "in [a, b).",
        py::arg("rows"), py::arg("cols"), py::arg("a") = 0, py::arg("b") = 1);

    m.def(
        "randn",
        [](idx rows, idx cols, realT mean, realT sigma) {
            return qpp::randn<cmat>(rows, cols, mean, sigma);
        },
        "Generates a random complex matrix with entries normally distributed "
        "in N(mean, sigma).",
        py::arg("rows"), py::arg("cols"), py::arg("mean") = 0,
        py::arg("sigma") = 1);

    // --- Quantum State and Operator Randomization ---

    m.def("randU", &qpp::randU,
          "Generates a random unitary matrix (Haar measure).",
          py::arg("D") = 2);

    m.def("randV", &qpp::randV, "Generates a random isometry matrix.",
          py::arg("Din"), py::arg("Dout"));

    m.def("randH", &qpp::randH, "Generates a random Hermitian matrix.",
          py::arg("D") = 2);

    m.def("randket", &qpp::randket,
          "Generates a random normalized pure state (ket).", py::arg("D") = 2);

    m.def("randrho", &qpp::randrho, "Generates a random density matrix.",
          py::arg("D") = 2);

    m.def("randkraus", py::overload_cast<idx, idx, idx>(&qpp::randkraus),
          "Generates a set of random Kraus operators (rectangular).",
          py::arg("N"), py::arg("Din"), py::arg("Dout"));

    m.def("randkraus", py::overload_cast<idx, idx>(&qpp::randkraus),
          "Generates a set of random Kraus operators (square).", py::arg("N"),
          py::arg("D") = 2);

    // --- Combinatorial and Statistical Randomization ---

    m.def("randperm", &qpp::randperm,
          "Generates a random permutation of [0, ..., N-1].", py::arg("N"));

    m.def("randprob", &qpp::randprob,
          "Generates a random probability vector uniformly over the simplex.",
          py::arg("N"));

    m.def("bernoulli", &qpp::bernoulli,
          "Generates a random boolean from a Bernoulli-p distribution.",
          py::arg("p") = 0.5);
}

#endif /* PYQPP_RANDOM_BIND_HPP_ */
