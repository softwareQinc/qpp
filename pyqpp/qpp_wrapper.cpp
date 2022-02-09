/*
 * This file is part of pyqpp.
 *
 * Copyright (c) 2019 - 2022 softwareQ Inc. All rights reserved.
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

/* Data types aliases */
using idx = qpp::idx;
using cmat = qpp::cmat;
using ket = qpp::ket;

/* Class aliases */
using Bit_circuit = qpp::Bit_circuit;
using Dynamic_bitset = qpp::Dynamic_bitset;
using QCircuit = qpp::QCircuit;
using QEngine = qpp::QEngine;

/* qpp::QEngine trampoline class for virtual methods */
class PyQEngine : public QEngine {
  public:
    using QEngine::QEngine;

    bool is_noisy() const override {
        PYBIND11_OVERRIDE(
            bool,    /* Return type */
            QEngine, /* Parent class */
            is_noisy /* Name of function in C++ (must match Python name) */
        );
    }

    QEngine& execute(idx reps = 1, bool clear_stats = true) override {
        PYBIND11_OVERRIDE(
            QEngine&, /* Return type */
            QEngine,  /* Parent class */
            execute,  /* Name of function in C++ (must match Python name) */
            reps,     /* Argument(s) */
            clear_stats);
    }
};

/* qpp::QNoisyEngine instantiator */
template <typename NoiseModel, typename... CtorTypeList>
void declare_noisy_engine(py::module& m, const std::string& type) {
    py::class_<NoiseModel>(m, type.c_str())
        .def(py::init<CtorTypeList...>())
        .def("get_d", &NoiseModel::get_d, "Qudit dimension")
        .def("get_Ks", &NoiseModel::get_Ks, "Vector of noise operators")
        .def("get_probs", &NoiseModel::get_probs,
             "Vector of probabilities corresponding to each noise operator")
        .def("get_last_idx", &NoiseModel::get_last_idx,
             "Index of the last occurring noise element")
        .def("get_last_p", &NoiseModel::get_last_p,
             "Probability of the last occurring noise element")
        .def("get_last_K", &NoiseModel::get_last_K,
             "Last occurring noise element");

    std::string pyname = "QNoisyEngine_" + type;
    py::class_<qpp::QNoisyEngine<NoiseModel>, QEngine>(m, pyname.c_str())
        .def(py::init<const QCircuit&, const NoiseModel&>(),
             py::keep_alive<1, 2>())
        .def(
            "get_noise_results",
            &qpp::QNoisyEngine<NoiseModel>::get_noise_results,
            "Vector of noise results obtained before every step in the circuit")
        .def("__copy__",
             [](const qpp::QNoisyEngine<NoiseModel>& self) {
                 return qpp::QNoisyEngine<NoiseModel>(self);
             })
        .def("__deepcopy__",
             [](const qpp::QNoisyEngine<NoiseModel>& self, py::dict) {
                 return qpp::QNoisyEngine<NoiseModel>(self);
             });

    m.def(
        "QNoisyEngine",
        [](const QCircuit& qc, const NoiseModel& nm) {
            return qpp::QNoisyEngine(qc, nm);
        },
        py::keep_alive<0, 1>());
}

PYBIND11_MODULE(pyqpp, m) {
    m.doc() =
        "Python 3 wrapper for Quantum++ (https://github.com/softwareQinc/qpp)";

    /* qpp::Dynamic_bitset */
    auto pyDynamic_bitset =
        py::class_<Dynamic_bitset>(m, "Dynamic_bitset")
            .def(py::init<idx>(), py::arg("n"))
            .def(py::init<std::string, char, char>(), py::arg("str"),
                 py::arg("zero") = '0', py::arg("one") = '1')
            .def("all", &Dynamic_bitset::all, "True if all of the bits are set")
            .def("any", &Dynamic_bitset::any, "True if any of the bits is set")
            .def("count", &Dynamic_bitset::count,
                 "Number of bits set to one in the bitset (Hamming weight)")
            .def("flip", py::overload_cast<idx>(&Dynamic_bitset::flip),
                 "Flips the bit at position pos", py::arg("pos"))
            .def("flip", py::overload_cast<>(&Dynamic_bitset::flip),
                 "Flips all bits")
            .def("get", &Dynamic_bitset::get,
                 "The value of the bit at position pos", py::arg("pos"))
            .def("none", &Dynamic_bitset::none,
                 "True if none of the bits are set")
            .def("__sub__", &Dynamic_bitset::operator-,
                 "Number of places the two bitsets differ (Hamming distance)")
            .def("rand", py::overload_cast<double>(&Dynamic_bitset::rand),
                 "Sets all bits according to a Bernoulli(p) distribution",
                 py::arg("p") = 0.5)
            .def("rand", py::overload_cast<idx, double>(&Dynamic_bitset::rand),
                 "Sets the bit at position pos according to a Bernoulli(p) "
                 "distribution",
                 py::arg("pos"), py::arg("p") = 0.5)
            .def("reset", py::overload_cast<idx>(&Dynamic_bitset::reset),
                 "Sets the bit at position pos to false", py::arg("pos"))
            .def("reset", py::overload_cast<>(&Dynamic_bitset::reset),
                 "Sets all bits to false")
            .def("set", py::overload_cast<idx, bool>(&Dynamic_bitset::set),
                 "Sets the bit at position pos", py::arg("pos"),
                 py::arg("value") = true)
            .def("set", py::overload_cast<>(&Dynamic_bitset::set),
                 "Sets all bits to true")
            .def("size", &Dynamic_bitset::size,
                 "Number of bits stored in the bitset")
            .def("storage_size", &Dynamic_bitset::storage_size,
                 "Size of the underlying storage space (in units of "
                 "qpp::Dynamic_bitset::value_type, unsigned int by default)")
            .def("to_string", &Dynamic_bitset::to_string,
                 "String representation", py::arg("zero") = '0',
                 py::arg("one") = '1')
            .def(py::self == py::self)
            .def(py::self != py::self)
            .def("__repr__",
                 [](const Dynamic_bitset& dbs) {
                     std::ostringstream oss;
                     oss << dbs;
                     return oss.str();
                 })
            .def(
                "__copy__",
                [](const Dynamic_bitset& self) { return Dynamic_bitset(self); })
            .def("__deepcopy__", [](const Dynamic_bitset& self, py::dict) {
                return Dynamic_bitset(self);
            });

    /* qpp::Bit_circuit */
    auto pyBit_circuit =
        py::class_<Bit_circuit, Dynamic_bitset>(m, "Bit_circuit")
            .def(py::init<idx>(), py::arg("n"))
            .def(py::init<std::string, char, char>(), py::arg("str"),
                 py::arg("zero") = '0', py::arg("one") = '1')
            .def(py::init<const Dynamic_bitset&>(), py::keep_alive<1, 2>())
            .def("CNOT", &Bit_circuit::CNOT, "Controlled-NOT gate",
                 py::arg("ctrl"), py::arg("target"))
            .def("FRED", &Bit_circuit::FRED, " Fredkin gate (Controlled-SWAP)",
                 py::arg("i"), py::arg("j"), py::arg("k"))
            .def("get_gate_count",
                 py::overload_cast<const std::string&>(
                     &Bit_circuit::get_gate_count, py::const_),
                 "Bit circuit gate count from name. Possible names are NOT "
                 "(X), CNOT, SWAP, TOF, FRED",
                 py::arg("name"))
            .def("get_gate_count",
                 py::overload_cast<>(&Bit_circuit::get_gate_count, py::const_),
                 "Total gate count")
            .def("get_gate_depth",
                 py::overload_cast<const std::string&>(
                     &Bit_circuit::get_gate_depth, py::const_),
                 "Bit circuit gate depth from name. Possible names are NOT "
                 "(X), CNOT, SWAP, TOF, FRED",
                 py::arg("name"))
            .def("get_gate_depth",
                 py::overload_cast<>(&Bit_circuit::get_gate_depth, py::const_),
                 "Total gate depth")
            .def("NOT", &Bit_circuit::NOT, "NOT gate (bit flip)", py::arg("i"))
            .def("reset", &Bit_circuit::reset,
                 "Resets the circuit to all-zero, clears all gates")
            .def("SWAP", &Bit_circuit::SWAP, "Swap gate", py::arg("i"),
                 py::arg("j"))
            .def("to_JSON", &Bit_circuit::to_JSON,
                 "Displays the bit circuit in JSON format",
                 py::arg("enclosed_in_curly_brackets") = true)
            .def("to_string", &Bit_circuit::to_string, "String representation",
                 py::arg("zero") = '0', py::arg("one") = '1')
            .def("TOF", &Bit_circuit::TOF, "Toffoli gate", py::arg("i"),
                 py::arg("j"), py::arg("k"))
            .def("X", &Bit_circuit::X, "NOT gate (bit flip)", py::arg("i"))
            .def(py::self == py::self)
            .def(py::self != py::self)
            .def("__repr__",
                 [](const Bit_circuit& bc) {
                     std::ostringstream oss;
                     oss << bc;
                     return oss.str();
                 })
            .def("__copy__",
                 [](const Bit_circuit& self) { return Bit_circuit(self); })
            .def("__deepcopy__", [](const Bit_circuit& self, py::dict) {
                return Bit_circuit(self);
            });

    /* qpp::QCircuit */
    auto pyQCircuit =
        py::class_<QCircuit>(m, "QCircuit")
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
                 "Vector of already measured qudit indexes")
            .def("get_measured_nd",
                 py::overload_cast<idx>(&QCircuit::get_measured_nd, py::const_),
                 "Whether qudit i was already measured non-destructively",
                 py::arg("i"))
            .def("get_measured_nd",
                 py::overload_cast<>(&QCircuit::get_measured_nd, py::const_),
                 "Vector of already measured non-destructively qudit indexes")
            .def("get_non_measured", &QCircuit::get_non_measured,
                 "Non-measured qudit indexes")
            .def("get_gate_count",
                 py::overload_cast<const cmat&>(&QCircuit::get_gate_count,
                                                py::const_),
                 "Gate count", py::arg("U"))
            .def("get_gate_count",
                 py::overload_cast<>(&QCircuit::get_gate_count, py::const_),
                 "Total gate count")
            .def("get_gate_depth",
                 py::overload_cast<const cmat&>(&QCircuit::get_gate_depth,
                                                py::const_),
                 "Gate depth", py::arg("U"))
            .def("get_gate_depth",
                 py::overload_cast<>(&QCircuit::get_gate_depth, py::const_),
                 "Total gate depth")
            .def("get_measurement_depth",
                 py::overload_cast<const cmat&>(
                     &QCircuit::get_measurement_depth, py::const_),
                 "Measurement depth", py::arg("V"))
            .def("get_measurement_depth",
                 py::overload_cast<>(&QCircuit::get_measurement_depth,
                                     py::const_),
                 "Total measurement depth")
            .def("get_depth", &QCircuit::get_depth,
                 "Quantum circuit description total depth")
            .def("get_measurement_count",
                 py::overload_cast<const cmat&>(
                     &QCircuit::get_measurement_count, py::const_),
                 "Measurement count", py::arg("V"))
            .def("get_measurement_count",
                 py::overload_cast<>(&QCircuit::get_measurement_count,
                                     py::const_),
                 "Total measurement count")
            .def("get_step_count", &QCircuit::get_step_count,
                 "Total (gates + measurements) count")
            .def("get_nop_count", &QCircuit::get_nop_count, "No-op count")
            .def("get_resources", &QCircuit::get_resources,
                 "Quantum circuit resources")
            .def("set_name", &QCircuit::set_name, "Sets name", py::arg("name"))
            .def("add_qudit", py::overload_cast<idx, idx>(&QCircuit::add_qudit),
                 "Adds n additional qudits before qudit pos", py::arg("n"),
                 py::arg("pos"))
            .def("add_qudit", py::overload_cast<idx>(&QCircuit::add_qudit),
                 "Adds n additional qudits after the last qudit",
                 py::arg("n") = 1)
            .def("add_dit", py::overload_cast<idx, idx>(&QCircuit::add_dit),
                 "Adds n additional classical dits before qudit pos",
                 py::arg("n"), py::arg("pos"))
            .def("add_dit", py::overload_cast<idx>(&QCircuit::add_dit),
                 "Adds n additional classical dits after the last qudit",
                 py::arg("n") = 1)
            .def("gate",
                 py::overload_cast<const cmat&, idx, std::string>(
                     &QCircuit::gate),
                 "Applies the single qudit gate U on single qudit i",
                 py::arg("U"), py::arg("i"), py::arg("name") = "")
            .def("gate",
                 py::overload_cast<const cmat&, idx, idx, std::string>(
                     &QCircuit::gate),
                 "Applies the two qudit gate U on qudits i and j", py::arg("U"),
                 py::arg("i"), py::arg("j"), py::arg("name") = "")
            .def("gate",
                 py::overload_cast<const cmat&, idx, idx, idx, std::string>(
                     &QCircuit::gate),
                 "Applies the three qudit gate U on qudits i, j and k",
                 py::arg("U"), py::arg("i"), py::arg("j"), py::arg("k"),
                 py::arg("name") = "")
            .def("gate_fan",
                 py::overload_cast<const cmat&, const std::vector<idx>&,
                                   std::string>(&QCircuit::gate_fan),
                 "Applies the single qudit gate U on every qudit listed in "
                 "target",
                 py::arg("U"), py::arg("target"), py::arg("name") = "")
            .def("gate_fan",
                 py::overload_cast<const cmat&, std::string>(
                     &QCircuit::gate_fan),
                 "Applies the single qudit gate U on all of the remaining "
                 "non-measured qudits",
                 py::arg("U"), py::arg("name") = "")
            .def("gate_joint", &QCircuit::gate_joint,
                 "Jointly applies the multiple-qudit gate U on the qudit "
                 "indexes specified by target",
                 py::arg("U"), py::arg("target"), py::arg("name") = "")
            .def("QFT",
                 py::overload_cast<const std::vector<idx>&, bool>(
                     &QCircuit::QFT),
                 "Applies the quantum Fourier transform on the qudit indexes "
                 "specified by target",
                 py::arg("target"), py::arg("swap") = true)
            .def("QFT", py::overload_cast<bool>(&QCircuit::QFT),
                 "Applies the quantum Fourier transform on all of remaining "
                 "non-measured qudits",
                 py::arg("swap") = true)
            .def("TFQ",
                 py::overload_cast<const std::vector<idx>&, bool>(
                     &QCircuit::TFQ),
                 "Applies the inverse quantum Fourier transform on the qudit "
                 "indexes specified by target",
                 py::arg("target"), py::arg("swap") = true)
            .def("TFQ", py::overload_cast<bool>(&QCircuit::TFQ),
                 "Applies the inverse quantum Fourier transform on all of "
                 "remaining non-measured qudits",
                 py::arg("swap") = true)
            .def("CTRL",
                 py::overload_cast<const cmat&, idx, idx, idx, std::string>(
                     &QCircuit::CTRL),
                 "Applies the single qudit controlled gate U", py::arg("U"),
                 py::arg("ctrl"), py::arg("target"), py::arg("shift") = 0,
                 py::arg("name") = "")
            .def("CTRL",
                 py::overload_cast<const cmat&, idx, const std::vector<idx>&,
                                   idx, std::string>(&QCircuit::CTRL),
                 "Applies the single qudit controlled gate U", py::arg("U"),
                 py::arg("ctrl"), py::arg("target"), py::arg("shift") = 0,
                 py::arg("name") = "")
            .def("CTRL",
                 py::overload_cast<const cmat&, const std::vector<idx>&, idx,
                                   const std::vector<idx>&, std::string>(
                     &QCircuit::CTRL),
                 "Applies the single qudit controlled gate U", py::arg("U"),
                 py::arg("ctrl"), py::arg("target"),
                 py::arg("shift") = std::vector<idx>(), py::arg("name") = "")
            .def("CTRL",
                 py::overload_cast<const cmat&, const std::vector<idx>&,
                                   const std::vector<idx>&,
                                   const std::vector<idx>&, std::string>(
                     &QCircuit::CTRL),
                 "Applies the single qudit controlled gate U", py::arg("U"),
                 py::arg("ctrl"), py::arg("target"),
                 py::arg("shift") = std::vector<idx>(), py::arg("name") = "")
            .def("CTRL_joint", &QCircuit::CTRL_joint,
                 "Jointly applies the multiple-qudit controlled gate U",
                 py::arg("U"), py::arg("ctrl"), py::arg("target"),
                 py::arg("shift") = std::vector<idx>(), py::arg("name") = "")
            .def("cCTRL",
                 py::overload_cast<const cmat&, idx, idx, idx, std::string>(
                     &QCircuit::cCTRL),
                 "Applies the single qudit controlled gate U with classical "
                 "control dit",
                 py::arg("U"), py::arg("ctrl_dit"), py::arg("target"),
                 py::arg("shift") = 0, py::arg("name") = "")
            .def("cCTRL",
                 py::overload_cast<const cmat&, idx, const std::vector<idx>&,
                                   idx, std::string>(&QCircuit::cCTRL),
                 "Applies the single qudit controlled gate U with classical "
                 "control dit",
                 py::arg("U"), py::arg("ctrl_dit"), py::arg("target"),
                 py::arg("shift") = 0, py::arg("name") = "")
            .def("cCTRL",
                 py::overload_cast<const cmat&, const std::vector<idx>&, idx,
                                   const std::vector<idx>&, std::string>(
                     &QCircuit::cCTRL),
                 "Applies the single qudit controlled gate U with classical "
                 "control dit",
                 py::arg("U"), py::arg("ctrl_dits"), py::arg("target"),
                 py::arg("shift") = std::vector<idx>(), py::arg("name") = "")
            .def("cCTRL",
                 py::overload_cast<const cmat&, const std::vector<idx>&,
                                   const std::vector<idx>&,
                                   const std::vector<idx>&, std::string>(
                     &QCircuit::cCTRL),
                 "Applies the single qudit controlled gate U with classical "
                 "control dit",
                 py::arg("U"), py::arg("ctrl_dits"), py::arg("target"),
                 py::arg("shift") = std::vector<idx>(), py::arg("name") = "")
            .def("cCTRL_joint", &QCircuit::cCTRL_joint,
                 "Jointly applies the multiple-qudit controlled gate with "
                 "classical control dits",
                 py::arg("U"), py::arg("ctrl_dits"), py::arg("target"),
                 py::arg("shift") = std::vector<idx>(), py::arg("name") = "")
            .def("measureZ",
                 py::overload_cast<idx, idx, bool, std::string>(
                     &QCircuit::measureZ),
                 "Z measurement of single qudit", py::arg("target"),
                 py::arg("c_reg"), py::arg("destructive") = true,
                 py::arg("name") = "")
            .def("measureZ",
                 py::overload_cast<const std::vector<idx>&, idx, bool,
                                   std::string>(&QCircuit::measureZ),
                 "Z measurement of multiple qudits", py::arg("target"),
                 py::arg("c_reg"), py::arg("destructive") = true,
                 py::arg("name") = "")
            .def("measureV",
                 py::overload_cast<const cmat&, idx, idx, bool, std::string>(
                     &QCircuit::measureV),
                 "Measurement of single qudit in the orthonormal basis "
                 "specified by the columns of matrix V",
                 py::arg("V"), py::arg("target"), py::arg("c_reg"),
                 py::arg("destructive") = true, py::arg("name") = "")
            .def("measureV",
                 py::overload_cast<const cmat&, const std::vector<idx>&, idx,
                                   bool, std::string>(&QCircuit::measureV),
                 "Measurement of multiple qudits in the orthonormal basis "
                 "specified by the columns of matrix V",
                 py::arg("V"), py::arg("target"), py::arg("c_reg"),
                 py::arg("destructive") = true, py::arg("name") = "")
            .def("discard",
                 py::overload_cast<idx, std::string>(&QCircuit::discard),
                 "Discards single qudit by measuring it destructively in the "
                 "Z-basis",
                 py::arg("target"), py::arg("name") = "")
            .def("discard",
                 py::overload_cast<const std::vector<idx>&, std::string>(
                     &QCircuit::discard),
                 "Discards multiple qudits by measuring them destructively in "
                 "the Z-basis",
                 py::arg("target"), py::arg("name") = "")
            .def("nop", &QCircuit::nop, "No operation (no-op)")
            .def("reset", py::overload_cast<idx, std::string>(&QCircuit::reset),
                 "Reset single qudit", py::arg("target"), py::arg("name") = "")
            .def("reset",
                 py::overload_cast<const std::vector<idx>&, std::string>(
                     &QCircuit::reset),
                 "Reset multiple qudits", py::arg("target"),
                 py::arg("name") = "")
            .def("replicate", &QCircuit::replicate,
                 "Replicates the quantum circuit description, in place",
                 py::arg("n"))
            .def("match_circuit_right", &QCircuit::match_circuit_right,
                 "Matches a quantum circuit description to the current one, "
                 "placed at the right (end) of the current one",
                 py::arg("other"), py::arg("target"), py::arg("pos_dit") = -1)
            .def("match_circuit_left", &QCircuit::match_circuit_left,
                 "Matches a quantum circuit description to the current one, "
                 "placed at the left (beginning) of the current one",
                 py::arg("other"), py::arg("target"), py::arg("pos_dit") = -1)
            .def("add_circuit", &QCircuit::add_circuit,
                 "Appends (glues) a quantum circuit description to the current "
                 "one",
                 py::arg("other"), py::arg("pos_qudit"),
                 py::arg("pos_dit") = -1)
            .def("kron", &QCircuit::kron,
                 "Kronecker product with another quantum circuit description, "
                 "in place",
                 py::arg("qc"))
            .def("adjoint", &QCircuit::adjoint,
                 "Adjoint quantum circuit description, in place")
            .def("is_clean_qudit", &QCircuit::is_clean_qudit,
                 "Whether qudit i in the quantum circuit description was used "
                 "before or not",
                 py::arg("i"))
            .def("is_clean_dit", &QCircuit::is_clean_dit,
                 "Whether classical dit i in the quantum circuit description "
                 "was used before or not",
                 py::arg("i"))
            .def("is_measurement_dit", &QCircuit::is_measurement_dit,
                 "Whether classical dit i in the quantum circuit description "
                 "was used to store the result of a measurement (either "
                 "destructive or non-destructive)",
                 py::arg("i"))
            .def("get_clean_qudits", &QCircuit::get_clean_qudits,
                 "Vector of clean qudits")
            .def("get_clean_dits", &QCircuit::get_clean_dits,
                 "Vector of clean classical dits")
            .def("get_dirty_qudits", &QCircuit::get_dirty_qudits,
                 "Vector of dirty qudits")
            .def("get_dirty_dits", &QCircuit::get_dirty_dits,
                 "Vector of dirty classical dits")
            .def("get_measurement_dits", &QCircuit::get_measurement_dits,
                 "Vector of classical dits that were used to store results of "
                 "measurements (either destructive or non-destructive)")
            .def("remove_clean_qudit", &QCircuit::remove_clean_qudit,
                 "Removes clean qudit and relabels the rest of the qudits "
                 "accordingly",
                 py::arg("target"))
            .def("remove_clean_dit", &QCircuit::remove_clean_dit,
                 "Removes clean classical dit and relabels the rest of the "
                 "classical dits accordingly",
                 py::arg("target"))
            .def("remove_clean_qudits", &QCircuit::remove_clean_qudits,
                 "Removes clean qudits and relabels the rest of the qudits "
                 "accordingly",
                 py::arg("target"))
            .def("remove_clean_dits", &QCircuit::remove_clean_dits,
                 "Removes clean classical dits and relabels the rest of the "
                 "classical dits accordingly",
                 py::arg("target"))
            .def("compress", &QCircuit::compress,
                 "Removes all clean qudits and relabels the rest of the qudits "
                 "accordingly",
                 py::arg("compress_dits") = false)
            .def("to_JSON", &QCircuit::to_JSON,
                 "Displays the quantum circuit description in JSON format",
                 py::arg("enclosed_in_curly_brackets") = true)
            .def(py::self == py::self)
            .def(py::self != py::self)
            .def("__repr__",
                 [](const QCircuit& qc) {
                     std::ostringstream oss;
                     oss << qc;
                     return oss.str();
                 })
            .def("__copy__",
                 [](const QCircuit& self) { return QCircuit(self); })
            .def("__deepcopy__",
                 [](const QCircuit& self, py::dict) { return QCircuit(self); });

    /* qpp::QCircuit related free functions */
    m.def("add_circuit", &qpp::add_circuit,
          "Appends (glues) the second quantum circuit description to the first "
          "one",
          py::arg("qc1"), py::arg("qc2"), py::arg("pos_qudit"),
          py::arg("pos_dit") = -1);
    m.def("adjoint", static_cast<QCircuit (*)(QCircuit)>(&qpp::adjoint),
          "Adjoint quantum circuit description", py::arg("qc"));
    m.def("kron",
          static_cast<QCircuit (*)(QCircuit, const QCircuit&)>(&qpp::kron),
          "Kronecker product between two quantum circuit descriptions",
          py::arg("qc1"), py::arg("qc2"));
    m.def("match_circuit_left", &qpp::match_circuit_left,
          "Matches the second quantum circuit description to the left "
          "(beginning) of the first one",
          py::arg("qc1"), py::arg("qc2"), py::arg("target"),
          py::arg("pos_dit") = -1);
    m.def("match_circuit_right", &qpp::match_circuit_right,
          "Matches the second quantum circuit description to the right (end) "
          "of the first one",
          py::arg("qc1"), py::arg("qc2"), py::arg("target"),
          py::arg("pos_dit") = -1);
    m.def("random_circuit_count", &qpp::random_circuit_count,
          "Random quantum circuit description generator for fixed gate count",
          py::arg("nq"), py::arg("num_steps"), py::arg("p_two"),
          py::arg("one_qudit_gate_set") = std::vector<cmat>{},
          py::arg("two_qudit_gate_set") = std::vector<cmat>{}, py::arg("d") = 2,
          py::arg("one_qudit_gate_names") = std::vector<std::string>{},
          py::arg("two_qudit_gate_names") = std::vector<std::string>{});
    m.def("random_circuit_depth", &qpp::random_circuit_depth,
          "Random quantum circuit description generator for fixed gate depth",
          py::arg("nq"), py::arg("num_steps"), py::arg("p_two"),
          py::arg("gate_depth") = cmat{},
          py::arg("one_qudit_gate_set") = std::vector<cmat>{},
          py::arg("two_qudit_gate_set") = std::vector<cmat>{}, py::arg("d") = 2,
          py::arg("one_qudit_gate_names") = std::vector<std::string>{},
          py::arg("two_qudit_gate_names") = std::vector<std::string>{});
    m.def("replicate", &qpp::replicate,
          "Replicates a quantum circuit description", py::arg("qc"),
          py::arg("n"));

    /* qpp::QCircuit::Resources */
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
        .def("__repr__", [](const QCircuit::Resources& r) {
            std::ostringstream oss;
            oss << r;
            return oss.str();
        });

    /* qpp::QEngine */
    py::class_<QEngine, PyQEngine>(m, "QEngine")
        .def(py::init<const QCircuit&>(), py::keep_alive<1, 2>())
        .def("get_psi", &QEngine::get_psi, "Underlying quantum state")
        .def("get_dits", &QEngine::get_dits, "Underlying classical dits")
        .def("get_dit", &QEngine::get_dit,
             "Underlying classical dit at position i", py::arg("i"))
        .def("get_probs", &QEngine::get_probs,
             "Underlying measurement outcome probabilities")
        .def("get_measured",
             py::overload_cast<idx>(&QEngine::get_measured, py::const_),
             "Whether qudit i was already measured destructively", py::arg("i"))
        .def("get_measured",
             py::overload_cast<>(&QEngine::get_measured, py::const_),
             "Vector of already measured qudit indexes")
        .def("get_non_measured", &QEngine::get_non_measured,
             "Non-measured qudit indexes")
        .def(
            "get_circuit", [](const QEngine& qe) { return qe.get_circuit(); },
            "Underlying quantum circuit description")
        .def("get_stats", &QEngine::get_stats,
             "Measurement statistics for multiple runs")
        .def("is_noisy", &QEngine::is_noisy, "Whether the engine is noisy")
        .def("set_dit", &QEngine::set_dit,
             "Sets the classical dit at position i", py::arg("i"),
             py::arg("value"))
        .def("set_dits", &QEngine::set_dits, "Set the classical dits",
             py::arg("dits"))
        .def("set_psi", &QEngine::set_psi, "Sets the underlying quantum state",
             py::arg("psi"))
        .def("reset_stats", &QEngine::reset_stats,
             "Resets the collected measurement statistics hash table")
        .def("reset", &QEngine::reset, "Resets the engine",
             py::arg("reset_stats") = true)
        .def("execute", py::overload_cast<idx, bool>(&QEngine::execute),
             "Executes the entire quantum circuit description",
             py::arg("reps") = 1, py::arg("clear_stats") = true)
        .def("to_JSON", &QEngine::to_JSON, "State of the engine in JSON format",
             py::arg("enclosed_in_curly_brackets") = true)
        .def("__repr__",
             [](const QEngine& qe) {
                 std::ostringstream oss;
                 oss << qe;
                 return oss.str();
             })
        .def("__copy__", [](const QEngine& self) { return QEngine(self); })
        .def("__deepcopy__",
             [](const QEngine& self, py::dict) { return QEngine(self); });

    /* qpp::QNoisyEngine instantiations with different noise models */
    declare_noisy_engine<qpp::QubitDepolarizingNoise, double>(
        m, "QubitDepolarizingNoise");
    declare_noisy_engine<qpp::QuditDepolarizingNoise, double, idx>(
        m, "QuditDepolarizingNoise");
    declare_noisy_engine<qpp::QubitPhaseFlipNoise, double>(
        m, "QubitPhaseFlipNoise");
    declare_noisy_engine<qpp::QubitBitFlipNoise, double>(m,
                                                         "QubitBitFlipNoise");
    declare_noisy_engine<qpp::QubitBitPhaseFlipNoise, double>(
        m, "QubitBitPhaseFlipNoise");

    /* OpenQASM interfacing */
    auto py_qasm = m.def_submodule("qasm");
    py_qasm.def("read_from_file", &qpp::qasm::read_from_file,
                "Get QCircuit representation of OpenQASM circuit");

    /* qpp::Gates */
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
    gates.def(
        "Rn",
        [](double theta, const std::array<double, 3>& n) {
            return qpp::gt.Rn(theta, n);
        },
        "Qubit rotation of theta about the 3-dimensional real (unit) vector n",
        py::arg("theta"), py::arg("n"));
    gates.def(
        "RX", [](double theta) { return qpp::gt.RX(theta); },
        "Qubit rotation of theta about the X axis", py::arg("theta"));
    gates.def(
        "RY", [](double theta) { return qpp::gt.RY(theta); },
        "Qubit rotation of theta about the Y axis", py::arg("theta"));
    gates.def(
        "RZ", [](double theta) { return qpp::gt.RZ(theta); },
        "Qubit rotation of theta about the Z axis", py::arg("theta"));
    gates.def(
        "Zd", [](idx D) { return qpp::gt.Zd(D); },
        "Generalized Z gate for qudits", py::arg("D") = 2);
    gates.def(
        "SWAPd", [](idx D) { return qpp::gt.SWAPd(D); }, "SWAP gate for qudits",
        py::arg("D") = 2);
    gates.def(
        "Fd", [](idx D) { return qpp::gt.Fd(D); },
        "Quantum Fourier transform gate for qudits", py::arg("D") = 2);
    gates.def(
        "MODMUL", [](idx a, idx N, idx n) { return qpp::gt.MODMUL(a, N, n); },
        "Modular multiplication gate for qubits", py::arg("a"), py::arg("N"),
        py::arg("n"));
    gates.def(
        "Xd", [](idx D) { return qpp::gt.Xd(D); },
        "Generalized X gate for qudits", py::arg("D") = 2);
    gates.def(
        "Id", [](idx D) { return qpp::gt.Id(D); }, "Identity gate",
        py::arg("D") = 2);
    gates.def(
        "get_name", [](const cmat& U) { return qpp::gt.get_name(U); },
        "Get the name of the most common qubit gates", py::arg("U"));

    /* qpp::States */
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
        "mes", [](idx d) { return qpp::st.mes(d); },
        "Maximally entangled state of 2 qudits", py::arg("d") = 2);
    states.def(
        "zero", [](idx n, idx d) { return qpp::st.zero(n, d); },
        "Zero state of n qudits", py::arg("n") = 1, py::arg("d") = 2);
    states.def(
        "one", [](idx n, idx d) { return qpp::st.one(n, d); },
        "One state of n qudits", py::arg("n") = 1, py::arg("d") = 2);
    states.def(
        "jn", [](idx j, idx n, idx d) { return qpp::st.jn(j, n, d); },
        "$|j\\rangle^{\\otimes n}$ state of n qudits", py::arg("j"),
        py::arg("n") = 1, py::arg("d") = 2);
    states.def(
        "j", [](idx j, idx d) { return qpp::st.j(j, d); },
        "$|j\\rangle$ computational basis state of a single qudit",
        py::arg("j"), py::arg("d") = 2);
    states.def(
        "plus", [](idx n) { return qpp::st.plus(n); }, "Plus state of n qubits",
        py::arg("n") = 1);
    states.def(
        "minus", [](idx n) { return qpp::st.minus(n); },
        "Minus state of n qubits", py::arg("n") = 1);

    /* Constants */
    m.attr("pi") = qpp::pi;
    m.attr("ee") = qpp::ee;

    /* Some free functions (non-exhaustive list */
    m.def("omega", &qpp::omega, "D-th root of unity", py::arg("D"));
    m.def("randU", &qpp::randU, "Generates a random unitary matrix",
          py::arg("D") = 2);
    /* template methods must be explicitly instantiated, some examples below */
    m.def(
        "transpose", [](const cmat& A) { return qpp::transpose(A); },
        "Transpose", py::arg("A"));
    m.def(
        "conjugate", [](const cmat& A) { return qpp::conjugate(A); },
        "Complex conjugate", py::arg("A"));
    m.def(
        "adjoint", [](const cmat& A) { return qpp::adjoint(A); }, "Adjoint",
        py::arg("A"));
    m.def(
        "inverse", [](const cmat& A) { return qpp::inverse(A); }, "Inverse",
        py::arg("A"));
    m.def(
        "trace", [](const cmat& A) { return qpp::trace(A); }, "trace",
        py::arg("A"));
    m.def(
        "det", [](const cmat& A) { return qpp::det(A); }, "Determinant",
        py::arg("A"));
    m.def(
        "logdet", [](const cmat& A) { return qpp::logdet(A); },
        "Logarithm of the determinant", py::arg("A"));
    m.def(
        "norm", [](const cmat& A) { return qpp::norm(A); }, "Frobenius norm",
        py::arg("A"));
    m.def(
        "kron", [](const cmat& A, const cmat& B) { return qpp::kron(A, B); },
        "Kronecker product", py::arg("A"), py::arg("B"));
    m.def("kron", static_cast<cmat (*)(const std::vector<cmat>&)>(&qpp::kron),
          "Kronecker product of a list of elements", py::arg("As"));
    m.def(
        "sum", [](const cmat& A) { return qpp::sum(A); }, "Element-wise sum",
        py::arg("A"));
    m.def(
        "sum", [](const std::vector<cmat>& As) { return qpp::sum(As); },
        "Sum of the elements of the list", py::arg("A"));
    m.def(
        "prod", [](const cmat& A) { return qpp::prod(A); },
        "Element-wise product", py::arg("As"));
    m.def(
        "prod", [](const std::vector<cmat>& As) { return qpp::prod(As); },
        "Products of the elements of the list", py::arg("As"));
}
