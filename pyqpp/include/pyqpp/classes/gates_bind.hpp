/*
 * This file is part of pyqpp.
 *
 * Copyright (c) 2017 - 2024 softwareQ Inc. All rights reserved.
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

#ifndef PYQPP_CLASSES_GATES_BIND_HPP_
#define PYQPP_CLASSES_GATES_BIND_HPP_

#include "pyqpp/pyqpp_common.hpp"

/* qpp::Gates */
inline void init_classes_gates(py::module_& m) {
    using namespace qpp;

    auto gates = m.def_submodule("gates");
    gates.attr("Id2") = qpp::gt.Id2;

    gates.attr("X") = qpp::gt.X;
    gates.attr("Y") = qpp::gt.Y;
    gates.attr("Z") = qpp::gt.Z;
    gates.attr("H") = qpp::gt.H;
    gates.attr("S") = qpp::gt.S;
    gates.attr("T") = qpp::gt.T;

    gates.attr("CNOT") = qpp::gt.CNOT;
    gates.attr("CZ") = qpp::gt.CZ;
    gates.attr("CNOTba") = qpp::gt.CNOTba;
    gates.attr("SWAP") = qpp::gt.SWAP;

    gates.attr("RXX") = qpp::gt.RXX;
    gates.attr("RYY") = qpp::gt.RYY;

    gates.attr("TOF") = qpp::gt.TOF;
    gates.attr("FRED") = qpp::gt.FRED;

    gates.def(
        "Fd", [](idx D) { return qpp::gt.Fd(D); },
        "Quantum Fourier transform gate for qudits", py::arg("D") = 2);
    gates.def(
        "get_name", [](const cmat& U) { return qpp::gt.get_name(U); },
        "Get the name of the most common qubit gates", py::arg("U"));
    gates.def(
        "Id", [](idx D) { return qpp::gt.Id(D); }, "Identity gate",
        py::arg("D") = 2);
    gates.def(
        "MODMUL", [](idx a, idx N, idx n) { return qpp::gt.MODMUL(a, N, n); },
        "Modular multiplication gate for qubits", py::arg("a"), py::arg("N"),
        py::arg("n"));
    gates.def(
        "Rn",
        [](realT theta, const std::array<realT, 3>& n) {
            return qpp::gt.Rn(theta, n);
        },
        "Qubit rotation of theta about the 3-dimensional real (unit) vector n",
        py::arg("theta"), py::arg("n"));
    gates.def(
        "RX", [](realT theta) { return qpp::gt.RX(theta); },
        "Qubit rotation of theta about the X axis", py::arg("theta"));
    gates.def(
        "RY", [](realT theta) { return qpp::gt.RY(theta); },
        "Qubit rotation of theta about the Y axis", py::arg("theta"));
    gates.def(
        "RZ", [](realT theta) { return qpp::gt.RZ(theta); },
        "Qubit rotation of theta about the Z axis", py::arg("theta"));
    gates.def(
        "SWAPd", [](idx D) { return qpp::gt.SWAPd(D); }, "SWAP gate for qudits",
        py::arg("D") = 2);
    gates.def(
        "Xd", [](idx D) { return qpp::gt.Xd(D); },
        "Generalized X gate for qudits", py::arg("D") = 2);
    gates.def(
        "Zd", [](idx D) { return qpp::gt.Zd(D); },
        "Generalized Z gate for qudits", py::arg("D") = 2);
}

#endif /* PYQPP_CLASSES_GATES_BIND_HPP_ */
