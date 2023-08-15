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

#ifndef PYQPP_CLASSES_STATES_BIND_HPP_
#define PYQPP_CLASSES_STATES_BIND_HPP_

/* qpp::States */
inline void init_classes_states(py::module_& m) {
    using namespace qpp;

    auto states = m.def_submodule("states");
    states.attr("x0") = qpp::st.x0;
    states.attr("x1") = qpp::st.x1;
    states.attr("y0") = qpp::st.y0;
    states.attr("y1") = qpp::st.y1;
    states.attr("z0") = qpp::st.z0;
    states.attr("z1") = qpp::st.z1;

    states.attr("b00") = qpp::st.b00;
    states.attr("b01") = qpp::st.b01;
    states.attr("b10") = qpp::st.b10;
    states.attr("b11") = qpp::st.b11;

    states.attr("pb00") = qpp::st.pb00;
    states.attr("pb01") = qpp::st.pb01;
    states.attr("pb10") = qpp::st.pb10;
    states.attr("pb11") = qpp::st.pb11;

    states.attr("GHZ") = qpp::st.GHZ;
    states.attr("W") = qpp::st.W;

    states.attr("pGHZ") = qpp::st.pGHZ;
    states.attr("pW") = qpp::st.pW;

    states.def(
        "j", [](idx j, idx d) { return qpp::st.j(j, d); },
        "$|j\\rangle$ computational basis state of a single qudit",
        py::arg("j"), py::arg("d") = 2);
    states.def(
        "jn", [](idx j, idx n, idx d) { return qpp::st.jn(j, n, d); },
        "$|j\\rangle^{\\otimes n}$ state of n qudits", py::arg("j"),
        py::arg("n") = 1, py::arg("d") = 2);
    states.def(
        "mes", [](idx d) { return qpp::st.mes(d); },
        "Maximally entangled state of 2 qudits", py::arg("d") = 2);
    states.def(
        "minus", [](idx n) { return qpp::st.minus(n); },
        "Minus state of n qubits", py::arg("n") = 1);
    states.def(
        "one", [](idx n, idx d) { return qpp::st.one(n, d); },
        "One state of n qudits", py::arg("n") = 1, py::arg("d") = 2);
    states.def(
        "plus", [](idx n) { return qpp::st.plus(n); }, "Plus state of n qubits",
        py::arg("n") = 1);
    states.def(
        "zero", [](idx n, idx d) { return qpp::st.zero(n, d); },
        "Zero state of n qudits", py::arg("n") = 1, py::arg("d") = 2);
}

#endif /* PYQPP_CLASSES_STATES_BIND_HPP_ */
