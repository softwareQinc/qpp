/*
 * This file is part of pyqpp.
 *
 * Copyright (c) 2019 - 2021 softwareQ Inc. All rights reserved.
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

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <sstream>

#include "qpp.h"

namespace py = pybind11;
using QCircuit = qpp::QCircuit;
using QEngine = qpp::QEngine;
using idx = qpp::idx;
using cmat = qpp::cmat;
using ket = qpp::ket;

PYBIND11_MODULE(pyqpp, m) {
    m.doc() = "Python wrapper for qpp (https://github.com/softwareQinc/qpp)";

    auto pyQCircuit = py::class_<QCircuit>(m, "QCircuit")
        .def(py::init<idx, idx, idx, std::string>(), py::arg("nq") = 1,
             py::arg("nc") = 0, py::arg("d") = 2, py::arg("name") = "")
        .def(py::init<const QCircuit&>())
        .def("get_nq", &QCircuit::get_nq, "Number of qudits")
        .def("get_nc", &QCircuit::get_nc, "Number of classical dits")
        .def("get_d", &QCircuit::get_d, "Qudit dimension")
        .def("get_name", &QCircuit::get_name, "Description name")
        .def("get_measured",
             py::overload_cast<idx>(&QCircuit::get_measured, py::const_),
             "Whether qudit i was already measured", py::arg("i"))
        .def("get_measured",
             py::overload_cast<>(&QCircuit::get_measured, py::const_),
             "Already measured qudit indexes")
        .def("get_non_measured", &QCircuit::get_non_measured,
             "Non-measured qudit indexes")
        .def("get_gate_count",
             py::overload_cast<const std::string&>(&QCircuit::get_gate_count, py::const_),
             "Gate count", py::arg("name"))
        .def("get_gate_count",
             py::overload_cast<>(&QCircuit::get_gate_count, py::const_),
             "Total gate count")
        .def("get_gate_depth",
             py::overload_cast<const std::string&>(&QCircuit::get_gate_depth, py::const_),
             "Gate depth", py::arg("name"))
        .def("get_gate_depth",
             py::overload_cast<>(&QCircuit::get_gate_depth, py::const_),
             "Total gate depth")
        .def("get_measurement_depth",
             py::overload_cast<const std::string&>(&QCircuit::get_measurement_depth, py::const_),
             "Measurement depth", py::arg("name"))
        .def("get_measurement_depth",
             py::overload_cast<>(&QCircuit::get_measurement_depth, py::const_),
             "Total measurement depth")
        .def("get_depth", &QCircuit::get_depth,
             "Quantum circuit description total depth")
        .def("get_measurement_count",
             py::overload_cast<const std::string&>(&QCircuit::get_measurement_count, py::const_),
             "Measurement count", py::arg("name"))
        .def("get_measurement_count",
             py::overload_cast<>(&QCircuit::get_measurement_count, py::const_),
             "Total measurement count")
        .def("get_step_count", &QCircuit::get_step_count,
             "Total (gates + measurements) count")
        .def("get_nop_count", &QCircuit::get_nop_count, "No-op count")
        .def("get_resources", &QCircuit::get_resources,
             "Quantum circuit resources")
        .def("set_name", &QCircuit::set_name, "Sets name")
        .def("add_qudit", py::overload_cast<idx, idx>(&QCircuit::add_qudit),
             "Adds n additional qudits before qudit pos", py::arg("n"),
             py::arg("pos"))
        .def("add_qudit", py::overload_cast<idx>(&QCircuit::add_qudit),
             "Adds n additional qudits after the last qudit", py::arg("n") = 1)
        .def("add_dit", py::overload_cast<idx, idx>(&QCircuit::add_dit),
             "Adds n additional classical dits before qudit pos", py::arg("n"),
             py::arg("pos"))
        .def("add_dit", py::overload_cast<idx>(&QCircuit::add_dit),
             "Adds n additional classical dits after the last qudit",
             py::arg("n") = 1)
        .def("gate",
             py::overload_cast<const cmat&, idx, std::string>(&QCircuit::gate),
             "Applies the single qudit gate U on single qudit i", py::arg("U"),
             py::arg("i"), py::arg("name") = "")
        .def("gate",
             py::overload_cast<const cmat&, idx, idx, std::string>(&QCircuit::gate),
             "Applies the two qudit gate U on qudits i and j", py::arg("U"),
             py::arg("i"), py::arg("j"), py::arg("name") = "")
        .def("gate",
             py::overload_cast<const cmat&, idx, idx, idx, std::string>(&QCircuit::gate),
             "Applies the three qudit gate U on qudits i, j and k", py::arg("U"),
             py::arg("i"), py::arg("j"), py::arg("k"), py::arg("name") = "")
        .def("gate_fan",
             py::overload_cast<const cmat&, const std::vector<idx>&, std::string>(&QCircuit::gate_fan),
             "Applies the single qudit gate U on every qudit listed in target",
             py::arg("U"), py::arg("target"), py::arg("name") = "")
        .def("gate_fan",
             py::overload_cast<const cmat&, std::string>(&QCircuit::gate_fan),
             "Applies the single qudit gate U on all of the remaining non-measured qudits",
             py::arg("U"), py::arg("name") = "")
        .def("gate_joint", &QCircuit::gate_joint,
             "Jointly applies the multiple-qudit gate U on the qudit indexes specified by target",
             py::arg("U"), py::arg("target"), py::arg("name") = "")
        .def("QFT",
             py::overload_cast<const std::vector<idx>&, bool>(&QCircuit::QFT),
             "Applies the quantum Fourier transform on the qudit indexes specified by target",
             py::arg("target"), py::arg("swap") = true)
        .def("QFT", py::overload_cast<bool>(&QCircuit::QFT),
             "Applies the quantum Fourier transform on all of remaining non-measured qudits",
             py::arg("swap") = true)
        .def("TFQ",
             py::overload_cast<const std::vector<idx>&, bool>(&QCircuit::TFQ),
             "Applies the inverse quantum Fourier transform on the qudit indexes specified by target",
             py::arg("target"), py::arg("swap") = true)
        .def("TFQ", py::overload_cast<bool>(&QCircuit::TFQ),
             "Applies the inverse quantum Fourier transform on all of remaining non-measured qudits",
             py::arg("swap") = true)
        .def("CTRL",
             py::overload_cast<const cmat&, idx, idx, idx, std::string>(&QCircuit::CTRL),
             "Applies the single qudit controlled gate U", py::arg("U"),
             py::arg("ctrl"), py::arg("target"), py::arg("shift") = 0,
             py::arg("name") = "")
        .def("CTRL",
             py::overload_cast<const cmat&, idx, const std::vector<idx>&, idx, std::string>(&QCircuit::CTRL),
             "Applies the single qudit controlled gate U", py::arg("U"),
             py::arg("ctrl"), py::arg("target"), py::arg("shift") = 0,
             py::arg("name") = "")
        .def("CTRL",
             py::overload_cast<const cmat&, const std::vector<idx>&, idx, const std::vector<idx>&, std::string>(&QCircuit::CTRL),
             "Applies the single qudit controlled gate U", py::arg("U"),
             py::arg("ctrl"), py::arg("target"), py::arg("shift") = std::vector<idx>(),
             py::arg("name") = "")
        .def("CTRL",
             py::overload_cast<const cmat&, const std::vector<idx>&, const std::vector<idx>&, const std::vector<idx>&, std::string>(&QCircuit::CTRL),
             "Applies the single qudit controlled gate U", py::arg("U"),
             py::arg("ctrl"), py::arg("target"), py::arg("shift") = std::vector<idx>(),
             py::arg("name") = "")
        .def("CTRL_joint", &QCircuit::CTRL_joint,
             "Jointly applies the multiple-qudit controlled gate U",
             py::arg("U"), py::arg("ctrl"), py::arg("target"),
             py::arg("shift") = std::vector<idx>(), py::arg("name") = "")
        .def("cCTRL",
             py::overload_cast<const cmat&, idx, idx, idx, std::string>(&QCircuit::cCTRL),
             "Applies the single qudit controlled gate U with classical control dit",
             py::arg("U"), py::arg("ctrl_dit"), py::arg("target"),
             py::arg("shift") = 0, py::arg("name") = "")
        .def("cCTRL",
             py::overload_cast<const cmat&, idx, const std::vector<idx>&, idx, std::string>(&QCircuit::cCTRL),
             "Applies the single qudit controlled gate U with classical control dit",
             py::arg("U"), py::arg("ctrl_dit"), py::arg("target"),
             py::arg("shift") = 0, py::arg("name") = "")
        .def("cCTRL",
             py::overload_cast<const cmat&, const std::vector<idx>&, idx, const std::vector<idx>&, std::string>(&QCircuit::cCTRL),
             "Applies the single qudit controlled gate U with classical control dit",
             py::arg("U"), py::arg("ctrl_dits"), py::arg("target"),
             py::arg("shift") = std::vector<idx>(), py::arg("name") = "")
        .def("cCTRL",
             py::overload_cast<const cmat&, const std::vector<idx>&, const std::vector<idx>&, const std::vector<idx>&, std::string>(&QCircuit::cCTRL),
             "Applies the single qudit controlled gate U with classical control dit",
             py::arg("U"), py::arg("ctrl_dits"), py::arg("target"),
             py::arg("shift") = std::vector<idx>(), py::arg("name") = "")
        .def("cCTRL_joint", &QCircuit::cCTRL_joint,
             "Jointly applies the multiple-qudit controlled gate with classical control dits",
             py::arg("U"), py::arg("ctrl_dits"), py::arg("target"),
             py::arg("shift") = std::vector<idx>(), py::arg("name") = "")
        .def("measureZ",
             py::overload_cast<idx, idx, bool, std::string>(&QCircuit::measureZ),
             "Z measurement of single qudit", py::arg("target"),
             py::arg("c_reg"), py::arg("destructive") = true,
             py::arg("name") = "")
        .def("measureZ",
             py::overload_cast<const std::vector<idx>&, idx, bool, std::string>(&QCircuit::measureZ),
             "Z measurement of multiple qudits", py::arg("target"),
             py::arg("c_reg"), py::arg("destructive") = true,
             py::arg("name") = "")
        .def("measureV",
             py::overload_cast<const cmat&, idx, idx, bool, std::string>(&QCircuit::measureV),
             "Measurement of single qudit in the orthonormal basis specified by the columns of matrix V",
             py::arg("V"), py::arg("target"), py::arg("c_reg"),
             py::arg("destructive") = true, py::arg("name") = "")
        .def("measureV",
             py::overload_cast<const cmat&, const std::vector<idx>&, idx, bool, std::string>(&QCircuit::measureV),
             "Measurement of multiple qudits in the orthonormal basis specified by the columns of matrix V",
             py::arg("V"), py::arg("target"), py::arg("c_reg"),
             py::arg("destructive") = true, py::arg("name") = "")
        .def("discard", py::overload_cast<idx, std::string>(&QCircuit::discard),
             "Discards single qudit by measuring it destructively in the Z-basis",
             py::arg("target"), py::arg("name") = "")
        .def("discard",
             py::overload_cast<const std::vector<idx>&, std::string>(&QCircuit::discard),
             "Discards multiple qudits by measuring them destructively in the Z-basis",
             py::arg("target"), py::arg("name") = "")
        .def("nop", &QCircuit::nop, "No operation (no-op)")
        .def("reset", py::overload_cast<idx, std::string>(&QCircuit::reset),
             "Reset single qudit", py::arg("target"), py::arg("name") = "")
        .def("reset",
             py::overload_cast<const std::vector<idx>&, std::string>(&QCircuit::reset),
             "Reset multiple qudits", py::arg("target"), py::arg("name") = "")
        .def("replicate", &QCircuit::replicate,
             "Replicates the circuit, in place")
        .def("match_circuit_right", &QCircuit::match_circuit_right,
             "Matches a quantum circuit description to the current one, placed at the right (end) of the current one",
             py::arg("other"), py::arg("target"), py::arg("pos_dit") = -1)
        .def("match_circuit_left", &QCircuit::match_circuit_left,
             "Matches a quantum circuit description to the current one, placed at the left (beginning) of the current one",
             py::arg("other"), py::arg("target"), py::arg("pos_dit") = -1)
        .def("add_circuit", &QCircuit::add_circuit,
             "Appends (glues) a quantum circuit description to the current one",
             py::arg("other"), py::arg("pos_qudit"), py::arg("pos_dit") = -1)
        .def("kron", &QCircuit::kron,
             "Kronecker product with another quantum circuit description, in place")
        .def("adjoint", &QCircuit::adjoint,
             "Adjoint quantum circuit description, in place")
        .def("is_clean_qudit", &QCircuit::is_clean_qudit,
             "Whether a qudit in the circuit was used before or not",
             py::arg("i"))
        .def("is_clean_dit", &QCircuit::is_clean_dit,
             "Whether a classical dit in the circuit was used before or not",
             py::arg("i"))
        .def("get_clean_qudits", &QCircuit::get_clean_qudits,
             "Vector of clean qudits")
        .def("get_clean_dits", &QCircuit::get_clean_dits,
             "Vector of clean classical dits")
        .def("remove_clean_qudit", &QCircuit::remove_clean_qudit,
             "Removes clean qudit and relabels the rest of the qudits accordingly",
             py::arg("target"))
        .def("remove_clean_dit", &QCircuit::remove_clean_dit,
             "Removes clean classical dit and relabels the rest of the classical dits accordingly",
             py::arg("target"))
        .def("remove_clean_qudits", &QCircuit::remove_clean_qudits,
             "Removes clean qudits and relabels the rest of the qudits accordingly",
             py::arg("target"))
        .def("remove_clean_dits", &QCircuit::remove_clean_dits,
             "Removes clean classical dits and relabels the rest of the classical dits accordingly",
             py::arg("target"))
        .def("compress", &QCircuit::compress,
             "Removes all clean qudits and relabels the rest of the qudits accordingly",
             py::arg("compress_dits") = false)
        .def("to_JSON", &QCircuit::to_JSON,
             "Displays the quantum circuit description in JSON format",
             py::arg("enclosed_in_curly_brackets") = true)
        .def(py::self == py::self)
        .def(py::self != py::self)
        .def("__repr__", [](const QCircuit& qc){
             std::ostringstream oss;
             oss << qc;
             return oss.str();
        });

    py::class_<QCircuit::Resources>(pyQCircuit, "Resources")
        .def_readonly("nq", &QCircuit::Resources::nq)
        .def_readonly("nc", &QCircuit::Resources::nc)
        .def_readonly("d", &QCircuit::Resources::d)
        .def_readonly("name", &QCircuit::Resources::name)
        .def_readonly("step_count", &QCircuit::Resources::step_count)
        .def_readonly("gate_count", &QCircuit::Resources::gate_count)
        .def_readonly("gate_depth", &QCircuit::Resources::gate_depth)
        .def_readonly("measurement_count",
                      &QCircuit::Resources::measurement_count)
        .def_readonly("measurement_depth",
                      &QCircuit::Resources::measurement_depth)
        .def_readonly("total_depth", &QCircuit::Resources::total_depth)
        .def("__repr__", [](const QCircuit::Resources& r){
             std::ostringstream oss;
             oss << r;
             return oss.str();
        });

    py::class_<QEngine>(m, "QEngine")
        .def(py::init<const QCircuit&>())
        .def(py::init<const QEngine&>())
        .def("get_psi", &QEngine::get_psi, "Underlying quantum state")
        .def("get_dits", &QEngine::get_dits, "Underlying classical dits")
        .def("get_dit", &QEngine::get_dit,
             "Underlying classical dit at position i", py::arg("i"))
        .def("get_probs", &QEngine::get_probs,
             "Underlying measurement outcome probabilities")
        .def("get_measured",
             py::overload_cast<idx>(&QEngine::get_measured, py::const_),
             "Whether qudit i was already measured", py::arg("i"))
        .def("get_measured",
             py::overload_cast<>(&QEngine::get_measured, py::const_),
             "Vector of already measured qudit indexes")
        .def("get_non_measured", &QEngine::get_non_measured,
             "Non-measured qudit indexes")
        /* Always return QCircuit copy */
        .def("get_circuit", [](const QEngine& qe) -> QCircuit {
                return qe.get_circuit();
            },
            "Underlying quantum circuit description")
        .def("get_stats", &QEngine::get_stats,
             "Measurement statistics for multiple runs")
        .def("is_noisy", &QEngine::is_noisy, "Whether the engine is noisy")
        .def("set_dit", &QEngine::set_dit,
             "Sets the classical dit at position i", py::arg("i"),
             py::arg("value"))
        .def("set_dits", &QEngine::set_dits, "Set the classical dits")
        .def("set_psi", &QEngine::set_psi, "Sets the underlying quantum state")
        .def("reset_stats", &QEngine::reset_stats,
             "Resets the collected measurement statistics hash table")
        .def("reset", &QEngine::reset, "Resets the engine",
             py::arg("reset_stats") = true)
        .def("execute", py::overload_cast<idx, bool>(&QEngine::execute),
             "Executes the entire quantum circuit description",
             py::arg("reps") = 1, py::arg("clear_stats") = true)
        .def("to_JSON", &QEngine::to_JSON, "State of the engine in JSON format",
             py::arg("enclosed_in_curly_brackets") = true)
        .def("__repr__", [](const QEngine& qe){
             std::ostringstream oss;
             oss << qe;
             return oss.str();
        });

    auto py_qasm = m.def_submodule("qasm");
    py_qasm.def("read_from_file", &qpp::qasm::read_from_file,
                "Get QCircuit representation of OpenQASM circuit");

    auto gates = m.def_submodule("gates");
    gates.attr("Id2") = qpp::gt.Id2;
    gates.attr("H") = qpp::gt.H;
    gates.attr("X") = qpp::gt.X;
    gates.attr("Y") = qpp::gt.Y;
    gates.attr("Z") = qpp::gt.Z;
    gates.attr("S") = qpp::gt.S;
    gates.attr("T") = qpp::gt.T;
    gates.attr("CNOT") = qpp::gt.CNOT;
    gates.attr("CZ") = qpp::gt.CZ;
    gates.attr("CNOTba") = qpp::gt.CNOTba;
    gates.attr("SWAP") = qpp::gt.SWAP;
    gates.attr("TOF") = qpp::gt.TOF;
    gates.attr("FRED") = qpp::gt.FRED;
    gates.def("Rn", [](double theta, const std::array<double, 3>& n) -> cmat {
            return qpp::gt.Rn(theta, n);
        },
        "Qubit rotation of theta about the 3-dimensional real (unit) vector n",
        py::arg("theta"), py::arg("n"));
    gates.def("RX", [](double theta) -> cmat { return qpp::gt.RX(theta); },
              "Qubit rotation of theta about the X axis",
              py::arg("theta"));
    gates.def("RY", [](double theta) -> cmat { return qpp::gt.RY(theta); },
              "Qubit rotation of theta about the Y axis",
              py::arg("theta"));
    gates.def("RZ", [](double theta) -> cmat { return qpp::gt.RZ(theta); },
              "Qubit rotation of theta about the Z axis",
              py::arg("theta"));
    gates.def("Zd", [](idx D) -> cmat { return qpp::gt.Zd(D); },
              "Generalized Z gate for qudits",
              py::arg("D") = 2);
    gates.def("SWAPd", [](idx D) -> cmat { return qpp::gt.SWAPd(D); },
              "SWAP gate for qudits",
              py::arg("D") = 2);
    gates.def("Fd", [](idx D) -> cmat { return qpp::gt.Fd(D); },
              "Quantum Fourier transform gate for qudits",
              py::arg("D") = 2);
    gates.def("MODMUL", [](idx a, idx N, idx n) -> cmat {
            return qpp::gt.MODMUL(a, N, n);
        },
        "Modular multiplication gate for qubits",
        py::arg("a"), py::arg("N"), py::arg("n"));
    gates.def("Xd", [](idx D) -> cmat { return qpp::gt.Xd(D); },
              "Generalized X gate for qudits",
              py::arg("D") = 2);
    gates.def("Id", [](idx D) -> cmat { return qpp::gt.Id(D); },
              "Identity gate", py::arg("D") = 2);
    gates.def("get_name", [](const cmat& U) -> std::string {
            return qpp::gt.get_name(U);
        },
        "Get the name of the most common qubit gates", py::arg("U"));

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
    states.def("mes", [](idx d) -> ket { return qpp::st.mes(d); },
        "Maximally entangled state of 2 qudits", py::arg("d") = 2);
    states.def("zero", [](idx n, idx d) -> ket { return qpp::st.zero(n, d); },
        "Zero state of n qudits", py::arg("n") = 1, py::arg("d") = 2);
    states.def("one", [](idx n, idx d) -> ket { return qpp::st.one(n, d); },
        "One state of n qudits", py::arg("n") = 1, py::arg("d") = 2);
    states.def("jn", [](idx j, idx n, idx d) -> ket {
            return qpp::st.jn(j, n, d);
        },
        "$|j\\rangle^{\\otimes n}$ state of n qudits", py::arg("j"),
        py::arg("n") = 1, py::arg("d") = 2);
    states.def("j", [](idx j, idx d) -> ket { return qpp::st.j(j, d); },
        "$|j\\rangle$ computational basis state of a single qudit",
        py::arg("j"), py::arg("d") = 2);
    states.def("plus", [](idx n) -> ket { return qpp::st.plus(n); },
        "Plus state of n qubits", py::arg("n") = 1);
    states.def("minus", [](idx n) -> ket { return qpp::st.minus(n); },
        "Minus state of n qubits", py::arg("n") = 1);

    m.def("transpose", [](const cmat& A) -> cmat { return qpp::transpose(A); },
          "Transpose");
    m.def("conjugate", [](const cmat& A) -> cmat { return qpp::conjugate(A); },
          "Complex conjugate");
    m.def("adjoint", [](const cmat& A) -> cmat { return qpp::adjoint(A); },
          "Adjoint");
    m.def("inverse", [](const cmat& A) -> cmat { return qpp::inverse(A); },
          "Inverse");
    m.def("trace", [](const cmat& A) -> qpp::cplx { return qpp::trace(A); },
          "trace");
    m.def("det", [](const cmat& A) -> qpp::cplx { return qpp::det(A); },
          "Determinant");
    m.def("logdet", [](const cmat& A) -> qpp::cplx { return qpp::logdet(A); },
          "Logarithm of the determinant");
    m.def("sum", [](const cmat& A) -> qpp::cplx { return qpp::sum(A); },
          "Element-wise sum");
    m.def("prod", [](const cmat& A) -> qpp::cplx { return qpp::prod(A); },
          "Element-wise product");
    m.def("norm", [](const cmat& A) -> double { return qpp::norm(A); },
          "Frobenius norm");

    m.def("randU", &qpp::randU, "Generates a random unitary matrix",
          py::arg("D") = 2);

    m.attr("pi") = qpp::pi;
    m.attr("ee") = qpp::ee;
    m.def("omega", &qpp::omega, "D-th root of unity", py::arg("D"));
}
