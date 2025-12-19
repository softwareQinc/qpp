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
 * \file <pyqpp/statistics_bind.hpp>
 * \brief Bindings for <qpp/statistics.hpp>
 */

#ifndef PYQPP_STATISTICS_BIND_HPP_
#define PYQPP_STATISTICS_BIND_HPP_

#include "pyqpp/pyqpp_common.hpp"

/* Bindings for <qpp/statistics.hpp> */
inline void init_statistics(py::module_& m) {
    using namespace qpp;
    // --- Probability Distributions ---

    m.def("uniform", &qpp::uniform,
          "Generates a uniform probability distribution vector of size N.",
          py::arg("N"));

    m.def("marginalX", &qpp::marginalX,
          "Computes the marginal distribution of X from a joint probability "
          "matrix.",
          py::arg("probXY"));

    m.def("marginalY", &qpp::marginalY,
          "Computes the marginal distribution of Y from a joint probability "
          "matrix.",
          py::arg("probXY"));

    // --- Statistical Moments ---

    m.def(
        "avg",
        [](const std::vector<realT>& prob, const std::vector<realT>& X) {
            return qpp::avg(prob, X);
        },
        "Computes the average E[X] = sum(p_i * x_i).", py::arg("prob"),
        py::arg("X"));

    m.def(
        "var",
        [](const std::vector<realT>& prob, const std::vector<realT>& X) {
            return qpp::var(prob, X);
        },
        "Computes the variance of random variable X.", py::arg("prob"),
        py::arg("X"));

    m.def(
        "sigma",
        [](const std::vector<realT>& prob, const std::vector<realT>& X) {
            return qpp::sigma(prob, X);
        },
        "Computes the standard deviation of random variable X.",
        py::arg("prob"), py::arg("X"));

    // --- Multivariate Statistics ---

    m.def(
        "cov",
        [](const rmat& probXY, const std::vector<realT>& X,
           const std::vector<realT>& Y) { return qpp::cov(probXY, X, Y); },
        "Computes the covariance of random variables X and Y.",
        py::arg("probXY"), py::arg("X"), py::arg("Y"));

    m.def(
        "cor",
        [](const rmat& probXY, const std::vector<realT>& X,
           const std::vector<realT>& Y) { return qpp::cor(probXY, X, Y); },
        "Computes the correlation coefficient of random variables X and Y.",
        py::arg("probXY"), py::arg("X"), py::arg("Y"));
}

#endif /* PYQPP_STATISTICS_BIND_HPP_ */
