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

#ifndef PYQPP_CLASSES_QCIRCUIT_BIND_HPP_
#define PYQPP_CLASSES_QCIRCUIT_BIND_HPP_

#include "pyqpp/pyqpp_common.h"

/* qpp::QCircuit and related free functions */
inline void init_classes_qcircuit(py::module_& m) {
    using namespace qpp;

    auto pyQCircuit =
        py::class_<QCircuit>(m, "QCircuit")
            .def(py::init<idx, idx, idx, std::optional<std::string>>(),
                 py::arg("nq") = 1, py::arg("nc") = 0, py::arg("d") = 2,
                 py::arg("name") = std::nullopt)
            .def(py::init<const QCircuit&>())

            .def("compose_circuit", &QCircuit::compose_circuit,
                 "Composes (appends) a quantum circuit description to the end "
                 "of the current one",
                 py::arg("other"), py::arg("pos_qudit"),
                 py::arg("pos_dit") = std::nullopt)
            .def("compose_CTRL_circuit", &QCircuit::compose_CTRL_circuit,
                 "Composes (appends) a controlled quantum circuit description "
                 "to the end of the current one, with the current instance "
                 "acting as the control",
                 py::arg("ctrl"), py::arg("qc_target"), py::arg("pos_qudit"),
                 py::arg("shift") = std::nullopt,
                 py::arg("pos_dit") = std::nullopt)
            .def("couple_circuit_left", &QCircuit::couple_circuit_left,
                 "Couples (in place) a quantum circuit description to the "
                 "current one, placed at the left (beginning) of the current "
                 "one",
                 py::arg("other"), py::arg("target"),
                 py::arg("pos_dit") = std::nullopt)
            .def("couple_circuit_right", &QCircuit::couple_circuit_right,
                 "Couples (in place) a quantum circuit description to the "
                 "current one, placed at the right (end) of the current one",
                 py::arg("other"), py::arg("target"),
                 py::arg("pos_dit") = std::nullopt)
            .def("add_dit", py::overload_cast<idx>(&QCircuit::add_dit),
                 "Adds n additional classical dits after the last qudit",
                 py::arg("n") = 1)
            .def("add_dit", py::overload_cast<idx, idx>(&QCircuit::add_dit),
                 "Adds n additional classical dits before qudit pos",
                 py::arg("n"), py::arg("pos"))
            .def("add_qudit", py::overload_cast<idx>(&QCircuit::add_qudit),
                 "Adds n additional qudits after the last qudit",
                 py::arg("n") = 1)
            .def("add_qudit", py::overload_cast<idx, idx>(&QCircuit::add_qudit),
                 "Adds n additional qudits before qudit pos", py::arg("n"),
                 py::arg("pos"))
            .def("adjoint", &QCircuit::adjoint,
                 "Adjoint quantum circuit description, in place")
            .def(
                "cCTRL",
                py::overload_cast<const cmat&, idx, idx, std::optional<idx>,
                                  std::optional<std::string>>(&QCircuit::cCTRL),
                "Applies the single qubit controlled gate U with classical "
                "control dit ctrl and target qudit target, i.e., cCTRL-U",
                py::arg("U"), py::arg("ctrl_dit"), py::arg("target"),
                py::arg("shift") = std::nullopt, py::arg("name") = std::nullopt)
            .def(
                "cCTRL",
                py::overload_cast<const cmat&, idx, const std::vector<idx>&,
                                  std::optional<idx>,
                                  std::optional<std::string>>(&QCircuit::cCTRL),
                "Jointly applies the single qudit controlled gate U with "
                "classical control dit ctrl on the qudit indexes specified "
                "by target, i.e., cCTRL-U_{joint}",
                py::arg("U"), py::arg("ctrl_dit"), py::arg("target"),
                py::arg("shift") = std::nullopt, py::arg("name") = std::nullopt)
            .def(
                "cCTRL",
                py::overload_cast<const cmat&, const std::vector<idx>&, idx,
                                  std::optional<std::vector<idx>>,
                                  std::optional<std::string>>(&QCircuit::cCTRL),
                "Applies the single qudit controlled gate U with multiple "
                "classical control dits listed in ctrl on the target qudit "
                "target, i.e., cCTRL-cCTRL-...-cCTRL-U",
                py::arg("U"), py::arg("ctrl_dits"), py::arg("target"),
                py::arg("shift") = std::nullopt, py::arg("name") = std::nullopt)
            .def(
                "cCTRL",
                py::overload_cast<const cmat&, const std::vector<idx>&,
                                  const std::vector<idx>&,
                                  std::optional<std::vector<idx>>,
                                  std::optional<std::string>>(&QCircuit::cCTRL),
                "Jointly applies the multiple-qudit controlled gate U with "
                "multiple classical control dits listed in ctrl on the qudit "
                "indexes specified by target, i.e., "
                "cCTRL-cCTRL-...-cCTRL-U_{joint}",
                py::arg("U"), py::arg("ctrl_dits"), py::arg("target"),
                py::arg("shift") = std::nullopt, py::arg("name") = std::nullopt)
            .def("cCTRL_fan",
                 py::overload_cast<const cmat&, const std::vector<idx>&,
                                   const std::vector<idx>&,
                                   std::optional<std::vector<idx>>,
                                   std::optional<std::string>>(
                     &QCircuit::cCTRL_fan),
                 "Applies the single qudit controlled gate U with multiple "
                 "classical control dits listed in ctrl on every qudit "
                 "listed in target, i.e., cCTRL-cCTRL-...-cCTRL-U-U-...-U",
                 py::arg("U"), py::arg("ctrl_dits"), py::arg("target"),
                 py::arg("shift") = std::nullopt,
                 py::arg("name") = std::nullopt)
            .def("cCTRL_fan",
                 py::overload_cast<const cmat&, idx, const std::vector<idx>&,
                                   std::optional<idx>,
                                   std::optional<std::string>>(
                     &QCircuit::cCTRL_fan),
                 "Applies the single qudit controlled gate U with classical "
                 "control dit ctrl on every qudit listed in target, "
                 "i.e.,cCTRL-U-U-...-U",
                 py::arg("U"), py::arg("ctrl_dit"), py::arg("target"),
                 py::arg("shift") = std::nullopt,
                 py::arg("name") = std::nullopt)
            .def("compress", &QCircuit::compress,
                 "Removes all clean qudits and relabels the rest of the qudits "
                 "accordingly",
                 py::arg("compress_dits") = false)
            .def("CTRL",
                 py::overload_cast<const cmat&, idx, idx, std::optional<idx>,
                                   std::optional<std::string>>(&QCircuit::CTRL),
                 "Applies the single qudit controlled gate U with control "
                 "qudit ctrl and target qudit target, i.e., CTRL-U",
                 py::arg("U"), py::arg("ctrl"), py::arg("target"),
                 py::arg("shift") = std::nullopt,
                 py::arg("name") = std::nullopt)
            .def("CTRL",
                 py::overload_cast<const cmat&, idx, const std::vector<idx>&,
                                   std::optional<idx>,
                                   std::optional<std::string>>(&QCircuit::CTRL),
                 "Jointly applies the single qudit controlled gate U with "
                 "control qudit ctrl on the qudit indexes specified by target, "
                 "i.e., CTRL-U_{joint}",
                 py::arg("U"), py::arg("ctrl"), py::arg("target"),
                 py::arg("shift") = std::nullopt,
                 py::arg("name") = std::nullopt)
            .def("CTRL",
                 py::overload_cast<const cmat&, const std::vector<idx>&, idx,
                                   std::optional<std::vector<idx>>,
                                   std::optional<std::string>>(&QCircuit::CTRL),
                 "Applies the multiple-qudit controlled gate U with multiple "
                 "control qudits listed in ctrl on the target qudit specified "
                 "by target, i.e., CTRL-CTRL-...-CTRL-U",
                 py::arg("U"), py::arg("ctrl"), py::arg("target"),
                 py::arg("shift") = std::nullopt,
                 py::arg("name") = std::nullopt)
            .def("CTRL",
                 py::overload_cast<const cmat&, const std::vector<idx>&,
                                   const std::vector<idx>&,
                                   std::optional<std::vector<idx>>,
                                   std::optional<std::string>>(&QCircuit::CTRL),
                 "Jointly applies the multiple-qudit controlled gate U with "
                 "multiple control qudits listed in ctrl on the qudit indexes "
                 "specified by target, i.e., CTRL-CTRL-...-CTRL-U_{joint}",
                 py::arg("U"), py::arg("ctrl"), py::arg("target"),
                 py::arg("shift") = std::nullopt,
                 py::arg("name") = std::nullopt)
            .def("CTRL_fan",
                 py::overload_cast<const cmat&, const std::vector<idx>&,
                                   const std::vector<idx>&,
                                   std::optional<std::vector<idx>>,
                                   std::optional<std::string>>(
                     &QCircuit::CTRL_fan),
                 "Applies the single qudit controlled gate U with multiple "
                 "control qudits listed in ctrl on every qudit listed in "
                 "target, i.e., CTRL-CTRL-...-CTRL-U-U-...-U",
                 py::arg("U"), py::arg("ctrl"), py::arg("target"),
                 py::arg("shift") = std::nullopt,
                 py::arg("name") = std::nullopt)
            .def("CTRL_fan",
                 py::overload_cast<const cmat&, idx, const std::vector<idx>&,
                                   std::optional<idx>,
                                   std::optional<std::string>>(
                     &QCircuit::CTRL_fan),
                 "Applies the single qudit controlled gate U with control "
                 "qudit ctrl on every qudit listed in target, i.e., "
                 "CTRL-U-U-...-U",
                 py::arg("U"), py::arg("ctrl"), py::arg("target"),
                 py::arg("shift") = std::nullopt,
                 py::arg("name") = std::nullopt)
            .def("discard",
                 py::overload_cast<idx, std::optional<std::string>>(
                     &QCircuit::discard),
                 "Discards single qudit by measuring it destructively in the "
                 "Z-basis",
                 py::arg("target"), py::arg("name") = "discard")
            .def("discard",
                 py::overload_cast<const std::vector<idx>&,
                                   std::optional<std::string>>(
                     &QCircuit::discard),
                 "Discards multiple qudits by measuring them destructively in "
                 "the Z-basis",
                 py::arg("target"), py::arg("name") = "discard")
            .def(
                "gate",
                py::overload_cast<const cmat&, idx, std::optional<std::string>>(
                    &QCircuit::gate),
                "Applies the single qudit gate U on single qudit i",
                py::arg("U"), py::arg("i"), py::arg("name") = std::nullopt)
            .def("gate",
                 py::overload_cast<const cmat&, idx, idx,
                                   std::optional<std::string>>(&QCircuit::gate),
                 "Applies the two qudit gate U on qudits i and j", py::arg("U"),
                 py::arg("i"), py::arg("j"), py::arg("name") = std::nullopt)
            .def("gate",
                 py::overload_cast<const cmat&, idx, idx, idx,
                                   std::optional<std::string>>(&QCircuit::gate),
                 "Applies the three qudit gate U on qudits i, j and k",
                 py::arg("U"), py::arg("i"), py::arg("j"), py::arg("k"),
                 py::arg("name") = std::nullopt)
            .def("gate",
                 py::overload_cast<const cmat&, const std::vector<idx>&,
                                   std::optional<std::string>>(&QCircuit::gate),
                 "Jointly applies the multiple-qudit gate U on the qudit "
                 "indexes specified by target",
                 py::arg("U"), py::arg("target"),
                 py::arg("name") = std::nullopt)
            .def("gate_fan",
                 py::overload_cast<const cmat&, const std::vector<idx>&,
                                   std::optional<std::string>>(
                     &QCircuit::gate_fan),
                 "Applies the single qudit gate U on every qudit listed in "
                 "target",
                 py::arg("U"), py::arg("target"),
                 py::arg("name") = std::nullopt)
            .def("gate_fan",
                 py::overload_cast<const cmat&, std::optional<std::string>>(
                     &QCircuit::gate_fan),
                 "Applies the single qudit gate U on all of the remaining "
                 "non-measured qudits",
                 py::arg("U"), py::arg("name") = std::nullopt)
            .def("get_clean_dits", &QCircuit::get_clean_dits,
                 "Vector of clean classical dits")
            .def("get_clean_qudits", &QCircuit::get_clean_qudits,
                 "Vector of clean qudits")
            .def("get_d", &QCircuit::get_d, "Qudit dimension")
            .def("get_depth", &QCircuit::get_depth,
                 "Quantum circuit description total depth")
            .def("get_dirty_dits", &QCircuit::get_dirty_dits,
                 "Vector of dirty classical dits")
            .def("get_dirty_qudits", &QCircuit::get_dirty_qudits,
                 "Vector of dirty qudits")
            .def("get_gate_count", &QCircuit::get_gate_count,
                 "(Total) Gate count", py::arg("U") = std::nullopt)
            .def("get_gate_depth", &QCircuit::get_gate_depth,
                 "(Total) Gate depth", py::arg("U") = std::nullopt)
            .def("get_measured", &QCircuit::get_measured,
                 "Vector of already measured qudit indexes")
            .def("get_measured_nd", &QCircuit::get_measured_nd,
                 "Vector of already measured non-destructively qudit indexes")
            .def("get_measurement_count", &QCircuit::get_measurement_count,
                 "(Total) Measurement count", py::arg("V") = std::nullopt)
            .def("get_measurement_depth", &QCircuit::get_measurement_depth,
                 "(Total) Measurement depth", py::arg("V") = std::nullopt)
            .def("get_measurement_dits", &QCircuit::get_measurement_dits,
                 "Vector of classical dits that were used to store results of "
                 "measurements (either destructive or non-destructive)")
            .def("get_name", &QCircuit::get_name, "Description name")
            .def("get_nc", &QCircuit::get_nc, "Number of classical dits")
            .def("get_non_measured", &QCircuit::get_non_measured,
                 "Non-measured qudit indexes")
            .def("get_nop_count", &QCircuit::get_nop_count, "No-op count")
            .def("get_nq", &QCircuit::get_nq, "Number of qudits")
            .def("get_resources", &QCircuit::get_resources,
                 "Quantum circuit resources")
            .def("get_step_count", &QCircuit::get_step_count,
                 "Total (gates + measurements) count")
            .def("has_measurements", &QCircuit::has_measurements,
                 "True if the quantum circuit description contains any "
                 "measurements, false otherwise")
            .def_static(
                "is_cCTRL", &QCircuit::is_cCTRL,
                "True if the gate step is a classically-controlled gate, "
                "false otherwise",
                py::arg("gate_step"))
            .def("is_clean_dit", &QCircuit::is_clean_dit,
                 "Whether classical dit i in the quantum circuit description "
                 "was used before or not",
                 py::arg("i"))
            .def("is_clean_qudit", &QCircuit::is_clean_qudit,
                 "Whether qudit i in the quantum circuit description was used "
                 "before or not",
                 py::arg("i"))
            .def_static(
                "is_CTRL", &QCircuit::is_CTRL,
                "True if the gate step is a controlled gate, false otherwise",
                py::arg("gate_step"))
            .def("is_measurement_dit", &QCircuit::is_measurement_dit,
                 "Whether classical dit i in the quantum circuit description "
                 "was used to store the result of a measurement (either "
                 "destructive or non-destructive)",
                 py::arg("i"))
            .def("kron", &QCircuit::kron,
                 "Kronecker product with another quantum circuit description, "
                 "in place",
                 py::arg("qc"))
            .def("measure",
                 py::overload_cast<idx, idx, bool, std::optional<std::string>>(
                     &QCircuit::measure),
                 "Z measurement of single qudit", py::arg("target"),
                 py::arg("c_reg"), py::arg("destructive") = true,
                 py::arg("name") = "mZ")
            .def("measure",
                 py::overload_cast<const std::vector<idx>&, idx, bool,
                                   std::optional<std::string>>(
                     &QCircuit::measure),
                 "Z measurement of multiple qudits", py::arg("target"),
                 py::arg("c_reg") = 0, py::arg("destructive") = true,
                 py::arg("name") = "mZ")
            .def("measure_all", &QCircuit::measure_all,
                 "Z measurement of all qudits", py::arg("c_reg") = 0,
                 py::arg("destructive") = true, py::arg("name") = "mZ")
            .def("measureV",
                 py::overload_cast<const cmat&, idx, idx, bool,
                                   std::optional<std::string>>(
                     &QCircuit::measureV),
                 "Measurement of single qudit in the orthonormal basis "
                 "specified by the columns of matrix V",
                 py::arg("V"), py::arg("target"), py::arg("c_reg"),
                 py::arg("destructive") = true, py::arg("name") = std::nullopt)
            .def("measureV",
                 py::overload_cast<const cmat&, const std::vector<idx>&, idx,
                                   bool, std::optional<std::string>>(
                     &QCircuit::measureV),
                 "Measurement of multiple qudits in the orthonormal basis "
                 "specified by the columns of matrix V",
                 py::arg("V"), py::arg("target"), py::arg("c_reg"),
                 py::arg("destructive") = true, py::arg("name") = std::nullopt)
            .def("nop", &QCircuit::nop, "No operation (no-op)")
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
            .def("remove_clean_dit", &QCircuit::remove_clean_dit,
                 "Removes clean classical dit and relabels the rest of the "
                 "classical dits accordingly",
                 py::arg("target"))
            .def("remove_clean_dits", &QCircuit::remove_clean_dits,
                 "Removes clean classical dits and relabels the rest of the "
                 "classical dits accordingly",
                 py::arg("target"))
            .def("remove_clean_qudit", &QCircuit::remove_clean_qudit,
                 "Removes clean qudit and relabels the rest of the qudits "
                 "accordingly",
                 py::arg("target"))
            .def("remove_clean_qudits", &QCircuit::remove_clean_qudits,
                 "Removes clean qudits and relabels the rest of the qudits "
                 "accordingly",
                 py::arg("target"))
            .def("removes_qudits", &QCircuit::removes_qudits,
                 "True if the quantum circuit description contains any "
                 "measurements that remove qudits, false otherwise")
            .def("replicate", &QCircuit::replicate,
                 "Replicates the quantum circuit description, in place",
                 py::arg("n"))
            .def(
                "reset",
                py::overload_cast<const std::vector<idx>&,
                                  std::optional<std::string>>(&QCircuit::reset),
                "Reset multiple qudits", py::arg("target"),
                py::arg("name") = "reset")
            .def("reset",
                 py::overload_cast<idx, std::optional<std::string>>(
                     &QCircuit::reset),
                 "Reset single qudit", py::arg("target"),
                 py::arg("name") = "reset")
            .def("set_name", &QCircuit::set_name, "Sets name", py::arg("name"))
            .def("TFQ", py::overload_cast<bool>(&QCircuit::TFQ),
                 "Applies the inverse quantum Fourier transform on all of "
                 "remaining non-measured qudits",
                 py::arg("swap") = true)
            .def("TFQ",
                 py::overload_cast<const std::vector<idx>&, bool>(
                     &QCircuit::TFQ),
                 "Applies the inverse quantum Fourier transform on the qudit "
                 "indexes specified by target",
                 py::arg("target"), py::arg("swap") = true)
            .def("to_JSON", &QCircuit::to_JSON,
                 "Displays the quantum circuit description in JSON format",
                 py::arg("enclosed_in_curly_brackets") = true)
            .def("was_measured", &QCircuit::was_measured,
                 "Whether qudit i was already measured", py::arg("i"))
            .def("was_measured_nd", &QCircuit::was_measured_nd,
                 "Whether qudit i was already measured non-destructively",
                 py::arg("i"))

            .def(py::self == py::self)
            .def(py::self != py::self)
            .def("__repr__",
                 [](const QCircuit& self) {
                     std::ostringstream oss;
                     oss << self;
                     return oss.str();
                 })
            .def("__copy__",
                 [](const QCircuit& self) { return QCircuit(self); })
            .def("__deepcopy__",
                 [](const QCircuit& self, py::dict) { return QCircuit(self); });

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
        .def("__repr__", [](const QCircuit::Resources& self) {
            std::ostringstream oss;
            oss << self;
            return oss.str();
        });

    /* qpp::QCircuit related free functions */
    m.def("compose_circuit", &qpp::compose_circuit,
          "Composes (appends) the second quantum circuit description to the "
          "end of the first one; qc_ctrl controls the qc_target.",
          py::arg("qc1"), py::arg("qc2"), py::arg("pos_qudit"),
          py::arg("name") = std::nullopt, py::arg("pos_dit") = std::nullopt);
    m.def("compose_CTRL_circuit", &qpp::compose_CTRL_circuit,
          "Composes (appends) the qc_target controlled quantum circuit "
          "description to the end of the qc_ctrl quantum circuit description",
          py::arg("qc_ctrl"), py::arg("ctrl"), py::arg("qc_target"),
          py::arg("pos_qudit"), py::arg("shift") = std::nullopt,
          py::arg("pos_dit") = std::nullopt, py::arg("name") = std::nullopt);
    m.def("couple_circuit_left", &qpp::couple_circuit_left,
          "Couples (in place) the second quantum circuit description to the "
          "left (beginning) of the first one",
          py::arg("qc1"), py::arg("qc2"), py::arg("target"),
          py::arg("name") = std::nullopt, py::arg("pos_dit") = std::nullopt);
    m.def("couple_circuit_right", &qpp::couple_circuit_right,
          "Couples (in place) the second quantum circuit description to the "
          "right (end) of the first one",
          py::arg("qc1"), py::arg("qc2"), py::arg("target"),
          py::arg("name") = std::nullopt, py::arg("pos_dit") = std::nullopt);
    m.def("adjoint",
          static_cast<QCircuit (*)(QCircuit, std::optional<std::string>)>(
              &qpp::adjoint),
          "Adjoint quantum circuit description", py::arg("qc"),
          py::arg("name") = std::nullopt);
    m.def("kron",
          static_cast<QCircuit (*)(QCircuit, const QCircuit&,
                                   std::optional<std::string>)>(&qpp::kron),
          "Kronecker product between two quantum circuit descriptions",
          py::arg("qc1"), py::arg("qc2"), py::arg("name") = std::nullopt);
    m.def("qpe_circuit", &qpp::qpe_circuit,
          "Quantum phase estimation circuit with n bits of precision",
          py::arg("U"), py::arg("n"), py::arg("omit_measurements") = true,
          py::arg("d") = 2, py::arg("name") = "qpe");
    m.def("random_circuit_count", &qpp::random_circuit_count,
          "Random quantum circuit description generator for fixed gate count",
          py::arg("nq"), py::arg("d"), py::arg("gate_count"),
          py::arg("p_two") = std::nullopt,
          py::arg("with_respect_to_gate") = std::nullopt,
          py::arg("one_qudit_gate_set") = std::nullopt,
          py::arg("two_qudit_gate_set") = std::nullopt,
          py::arg("one_qudit_gate_names") = std::nullopt,
          py::arg("two_qudit_gate_names") = std::nullopt);
    m.def("random_circuit_depth", &qpp::random_circuit_depth,
          "Random quantum circuit description generator for fixed gate depth",
          py::arg("nq"), py::arg("d"), py::arg("gate_depth"),
          py::arg("p_two") = std::nullopt,
          py::arg("with_respect_to_gate") = std::nullopt,
          py::arg("one_qudit_gate_set") = std::nullopt,
          py::arg("two_qudit_gate_set") = std::nullopt,
          py::arg("one_qudit_gate_names") = std::nullopt,
          py::arg("two_qudit_gate_names") = std::nullopt);
    m.def("replicate", &qpp::replicate,
          "Replicates a quantum circuit description", py::arg("qc"),
          py::arg("n"), py::arg("name") = std::nullopt);
}

#endif /* PYQPP_CLASSES_QCIRCUIT_BIND_HPP_ */
