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
 * \file <pyqpp/functions_bind.hpp>
 * \brief Bindings for <qpp/functions.hpp>
 */

#ifndef PYQPP_FUNCTIONS_BIND_HPP_
#define PYQPP_FUNCTIONS_BIND_HPP_

#include "pyqpp/pyqpp_common.hpp"

inline void init_functions(py::module_& m) {
    using namespace qpp;

    /* Template methods must be explicitly instantiated, some examples below */
    // --- Eigen Wrappers ---
    m.def(
        "transpose", [](const cmat& A) { return qpp::transpose(A); },
        "Transpose");
    m.def(
        "conjugate", [](const cmat& A) { return qpp::conjugate(A); },
        "Complex conjugate");
    m.def(
        "adjoint", [](const cmat& A) { return qpp::adjoint(A); },
        "Adjoint (Hermitian conjugate)");
    m.def("inverse", [](const cmat& A) { return qpp::inverse(A); }, "Inverse");
    m.def("trace", [](const cmat& A) { return qpp::trace(A); }, "Trace");
    m.def("det", [](const cmat& A) { return qpp::det(A); }, "Determinant");
    m.def(
        "logdet", [](const cmat& A) { return qpp::logdet(A); },
        "Logarithm of the determinant");
    m.def(
        "sum", [](const cmat& A) { return qpp::sum(A); },
        "Element-wise sum of matrix");
    m.def(
        "prod", [](const cmat& A) { return qpp::prod(A); },
        "Element-wise product of matrix");
    m.def("norm", [](const cmat& A) { return qpp::norm(A); }, "Frobenius norm");
    m.def(
        "normalize", [](const cmat& A) { return qpp::normalize(A); },
        "Normalize state vector or density matrix");

    // --- Decompositions ---
    m.def(
        "eig", [](const cmat& A) { return qpp::eig(A); },
        "Full eigen decomposition");
    m.def("evals", [](const cmat& A) { return qpp::evals(A); }, "Eigenvalues");
    m.def(
        "evects", [](const cmat& A) { return qpp::evects(A); }, "Eigenvectors");
    m.def(
        "heig", [](const cmat& A) { return qpp::heig(A); },
        "Full eigen decomposition of Hermitian matrix");
    m.def(
        "hevals", [](const cmat& A) { return qpp::hevals(A); },
        "Hermitian eigenvalues");
    m.def(
        "hevects", [](const cmat& A) { return qpp::hevects(A); },
        "Hermitian eigenvectors");
    m.def(
        "svd", [](const cmat& A) { return qpp::svd(A); },
        "Full singular value decomposition");
    m.def(
        "svals", [](const cmat& A) { return qpp::svals(A); },
        "Singular values");
    m.def(
        "svdU", [](const cmat& A) { return qpp::svdU(A); },
        "Left singular vectors");
    m.def(
        "svdV", [](const cmat& A) { return qpp::svdV(A); },
        "Right singular vectors");

    // --- Matrix Functional Calculus ---
    m.def(
        "sqrtm", [](const cmat& A) { return qpp::sqrtm(A); },
        "Matrix square root");
    m.def(
        "absm", [](const cmat& A) { return qpp::absm(A); },
        "Matrix absolute value");
    m.def(
        "expm", [](const cmat& A) { return qpp::expm(A); },
        "Matrix exponential");
    m.def(
        "logm", [](const cmat& A) { return qpp::logm(A); }, "Matrix logarithm");
    m.def("sinm", [](const cmat& A) { return qpp::sinm(A); }, "Matrix sine");
    m.def("cosm", [](const cmat& A) { return qpp::cosm(A); }, "Matrix cosine");
    m.def(
        "spectralpowm",
        [](const cmat& A, cplx z) { return qpp::spectralpowm(A, z); },
        "Matrix power via spectral decomposition");
    m.def(
        "powm", [](const cmat& A, idx n) { return qpp::powm(A, n); },
        "Fast matrix power (square-and-multiply)");
    m.def(
        "schatten", [](const cmat& A, realT p) { return qpp::schatten(A, p); },
        "Schatten matrix norm");
    m.def(
        "funm",
        [](const cmat& A, cplx (*f)(const cplx&)) { return qpp::funm(A, f); },
        "Applies the scalar function f to the eigenvalues of matrix A",
        py::arg("A"), py::arg("f"));

    // --- Kronecker and Direct Sum ---
    m.def(
        "kron", [](const std::vector<cmat>& As) { return qpp::kron(As); },
        "Kronecker product of multiple matrices");
    m.def(
        "kronpow", [](const cmat& A, idx n) { return qpp::kronpow(A, n); },
        "Kronecker power");
    m.def(
        "dirsum", [](const std::vector<cmat>& As) { return qpp::dirsum(As); },
        "Direct sum of multiple matrices");
    m.def(
        "dirsumpow", [](const cmat& A, idx n) { return qpp::dirsumpow(A, n); },
        "Direct sum power");

    // --- Other Matrix Operations ---
    m.def(
        "reshape",
        [](const cmat& A, idx rows, idx cols) {
            return qpp::reshape(A, rows, cols);
        },
        "Reshape matrix");
    m.def(
        "comm", [](const cmat& A, const cmat& B) { return qpp::comm(A, B); },
        "Commutator");
    m.def(
        "anticomm",
        [](const cmat& A, const cmat& B) { return qpp::anticomm(A, B); },
        "Anticommutator");
    m.def(
        "prj", [](const cmat& A) { return qpp::prj(A); },
        "Projector onto state vector");
    m.def(
        "grams", [](const std::vector<cmat>& As) { return qpp::grams(As); },
        "Gram-Schmidt orthogonalization (from list)");
    m.def(
        "grams", [](const cmat& A) { return qpp::grams(A); },
        "Gram-Schmidt orthogonalization (from columns)");

    // --- Indices and Multi-indices ---
    m.def("n2multiidx", &qpp::n2multiidx<idx>,
          "Integer to multi-index conversion");
    m.def("multiidx2n", &qpp::multiidx2n<idx>,
          "Multi-index to integer conversion");

    // --- Projectors ---
    m.def("mprj",
          py::overload_cast<const std::vector<idx>&, const std::vector<idx>&>(
              &qpp::mprj),
          "Projector onto multi-partite qudit state");
    m.def("mprj", py::overload_cast<const std::vector<idx>&, idx>(&qpp::mprj),
          py::arg("mask"), py::arg("d") = 2,
          "Projector onto multi-partite qudit state (uniform d)");

    // --- STL-like Container Functions ---
    m.def(
        "abssq", [](const std::vector<cplx>& v) { return qpp::abssq(v); },
        "Absolute values squared of vector");
    m.def(
        "abssq", [](const cmat& A) { return qpp::abssq(A); },
        "Absolute values squared of matrix elements");
    m.def(
        "sum",
        [](const std::vector<cplx>& v) { return qpp::sum(v.begin(), v.end()); },
        "Sum of vector elements");
    m.def(
        "prod", [](const std::vector<cplx>& v) { return qpp::prod(v); },
        "Product of vector elements");

    // --- Miscellaneous ---
    m.def(
        "rho2pure", [](const cmat& A) { return qpp::rho2pure(A); },
        "Finds pure state representation of rank-1 matrix");
    m.def("complement", &qpp::complement, "Complement of a subsystem");
    m.def(
        "hash_eigen",
        [](const cmat& A, std::size_t seed) {
            return qpp::hash_eigen(A, seed);
        },
        py::arg("A"), py::arg("seed") = 0, "Hash an Eigen expression");

    // --- State Generation ---
    m.def("mket",
          py::overload_cast<const std::vector<idx>&, const std::vector<idx>&>(
              &qpp::mket),
          "Creates a multi-partite ket from indices and dimensions",
          py::arg("states"), py::arg("dims"));

    m.def("mket", py::overload_cast<const std::vector<idx>&, idx>(&qpp::mket),
          "Creates a multi-partite ket from indices with uniform dimension d",
          py::arg("states"), py::arg("d") = 2);

    m.def(
        "zket2dits",
        [](const cmat& psi, const std::vector<idx>& dims, realT precision) {
            return qpp::zket2dits(psi, dims, precision);
        },
        "Extracts dits from a computational basis state given specific "
        "subsystem dimensions.\n"
        "Returns None if psi is not a computational basis state.",
        py::arg("psi"), py::arg("dims"), py::arg("precision") = 1e-12);
    m.def(
        "zket2dits",
        [](const cmat& psi, idx d, realT precision) {
            return qpp::zket2dits(psi, d, precision);
        },
        "Extracts dits from a computational basis state assuming uniform "
        "subsystem dimension d.\n"
        "Returns None if psi is not a computational basis state.",
        py::arg("psi"), py::arg("d") = 2, py::arg("precision") = 1e-12);

    // --- Bloch Sphere Operations ---
    m.def(
        "rho2bloch", [](const cmat& A) { return qpp::rho2bloch(A); },
        "Computes the 3-dimensional real Bloch vector corresponding to the "
        "qubit density matrix A",
        py::arg("A"));
    m.def("bloch2rho", &qpp::bloch2rho,
          "Computes the density matrix corresponding to the 3-dimensional real "
          "Bloch vector r",
          py::arg("r"));

    // --- Formatting ---
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
