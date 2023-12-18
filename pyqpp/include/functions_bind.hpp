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

#ifndef PYQPP_FUNCTIONS_BIND_HPP_
#define PYQPP_FUNCTIONS_BIND_HPP_

/* Some free functions (non-exhaustive list) from functions.hpp */
inline void init_functions(py::module_& m) {
    using namespace qpp;

    /* Template methods must be explicitly instantiated, some examples below */
    m.def(
        "adjoint", [](const cmat& A) { return qpp::adjoint(A); }, "Adjoint",
        py::arg("A"));
    m.def(
        "conjugate", [](const cmat& A) { return qpp::conjugate(A); },
        "Complex conjugate", py::arg("A"));
    m.def(
        "det", [](const cmat& A) { return qpp::det(A); }, "Determinant",
        py::arg("A"));
    m.def(
        "inverse", [](const cmat& A) { return qpp::inverse(A); }, "Inverse",
        py::arg("A"));
    m.def(
        "kron", [](const cmat& A, const cmat& B) { return qpp::kron(A, B); },
        "Kronecker product", py::arg("A"), py::arg("B"));
    m.def("kron", static_cast<cmat (*)(const std::vector<cmat>&)>(&qpp::kron),
          "Kronecker product of a list of elements", py::arg("As"));
    m.def(
        "logdet", [](const cmat& A) { return qpp::logdet(A); },
        "Logarithm of the determinant", py::arg("A"));
    m.def(
        "norm", [](const cmat& A) { return qpp::norm(A); }, "Frobenius norm",
        py::arg("A"));
    m.def(
        "prod", [](const cmat& A) { return qpp::prod(A); },
        "Element-wise product", py::arg("As"));
    m.def(
        "prod", [](const std::vector<cmat>& As) { return qpp::prod(As); },
        "Products of the elements of the list", py::arg("As"));
    m.def(
        "sum", [](const cmat& A) { return qpp::sum(A); }, "Element-wise sum",
        py::arg("A"));
    m.def(
        "sum", [](const std::vector<cmat>& As) { return qpp::sum(As); },
        "Sum of the elements of the list", py::arg("A"));
    m.def(
        "trace", [](const cmat& A) { return qpp::trace(A); }, "trace",
        py::arg("A"));
    m.def(
        "transpose", [](const cmat& A) { return qpp::transpose(A); },
        "Transpose", py::arg("A"));
    m.def(
        "dirac",
        [](const cmat& A, std::vector<idx> dims_rows,
           std::vector<idx> dims_cols) {
            return qpp::dirac(A, dims_rows, dims_cols);
        },
        "Dirac notation", py::arg("A"), py::arg("dims_rows"),
        py::arg("dims_cols"));
    m.def(
        "dirac", [](const cmat& A, idx d = 2) { return qpp::dirac(A, d); },
        "Dirac notation", py::arg("A"), py::arg("d") = 2);
}

#endif /* PYQPP_FUNCTIONS_BIND_HPP_ */
