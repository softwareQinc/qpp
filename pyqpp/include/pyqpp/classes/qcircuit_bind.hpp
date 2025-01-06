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

#ifndef PYQPP_CLASSES_QCIRCUIT_BIND_HPP_
#define PYQPP_CLASSES_QCIRCUIT_BIND_HPP_

#include <pybind11/functional.h>

#include "pyqpp/pyqpp_common.hpp"

/* qpp::QCircuit and related free functions */
inline void init_classes_qcircuit(py::module_& m) {
    using namespace qpp;

    auto pyQCircuit = py::class_<QCircuit>(m, "QCircuit");
    pyQCircuit.def(py::init<idx, idx, idx, std::optional<std::string>>(),
                   py::arg("nq") = 1, py::arg("nc") = 0, py::arg("d") = 2,
                   py::arg("name") = std::nullopt);
    pyQCircuit.def(py::init<const QCircuit&>());

    pyQCircuit.def(
        "compose_circuit", &QCircuit::compose_circuit,
        "Composes (appends) a quantum circuit description to the end "
        "of the current one",
        py::arg("other"), py::arg("pos_qudit"),
        py::arg("pos_dit") = std::nullopt);
    pyQCircuit.def(
        "compose_CTRL_circuit", &QCircuit::compose_CTRL_circuit,
        "Composes (appends) a controlled quantum circuit description "
        "to the end of the current one, with the current instance "
        "acting as the control",
        py::arg("ctrl"), py::arg("qc_target"), py::arg("pos_qudit"),
        py::arg("shift") = std::nullopt, py::arg("pos_dit") = std::nullopt);
    pyQCircuit.def("couple_circuit_left", &QCircuit::couple_circuit_left,
                   "Couples (in place) a quantum circuit description to the "
                   "current one, placed at the left (beginning) of the current "
                   "one",
                   py::arg("other"), py::arg("target"),
                   py::arg("pos_dit") = std::nullopt);
    pyQCircuit.def("couple_circuit_right", &QCircuit::couple_circuit_right,
                   "Couples (in place) a quantum circuit description to the "
                   "current one, placed at the right (end) of the current one",
                   py::arg("other"), py::arg("target"),
                   py::arg("pos_dit") = std::nullopt);
    pyQCircuit.def("add_dit", py::overload_cast<idx>(&QCircuit::add_dit),
                   "Adds n additional classical dits after the last qudit",
                   py::arg("n") = 1);
    pyQCircuit.def("add_dit", py::overload_cast<idx, idx>(&QCircuit::add_dit),
                   "Adds n additional classical dits before qudit pos",
                   py::arg("n"), py::arg("pos"));
    pyQCircuit.def("cond_if",
                   py::overload_cast<std::function<bool(std::vector<idx>)>>(
                       &QCircuit::cond_if),
                   "Adds conditional if", py::arg("pred"));
    pyQCircuit.def("cond_while",
                   py::overload_cast<std::function<bool(std::vector<idx>)>>(
                       &QCircuit::cond_while),
                   "Adds conditional while", py::arg("pred"));
    pyQCircuit.def("cond_else", &QCircuit::cond_else, "Adds conditional else");
    pyQCircuit.def("cond_end", &QCircuit::cond_end, "Adds conditional end");
    pyQCircuit.def("add_qudit", py::overload_cast<idx>(&QCircuit::add_qudit),
                   "Adds n additional qudits after the last qudit",
                   py::arg("n") = 1);
    pyQCircuit.def("add_qudit",
                   py::overload_cast<idx, idx>(&QCircuit::add_qudit),
                   "Adds n additional qudits before qudit pos", py::arg("n"),
                   py::arg("pos"));
    pyQCircuit.def("adjoint", &QCircuit::adjoint,
                   "Adjoint quantum circuit description, in place");
    pyQCircuit.def(
        "cCTRL",
        py::overload_cast<const cmat&, idx, idx, std::optional<idx>,
                          std::optional<std::string>>(&QCircuit::cCTRL),
        "Applies the single qudit controlled gate U with classical "
        "control dit ctrl and target qudit target, i.e., cCTRL-U",
        py::arg("U"), py::arg("ctrl_dit"), py::arg("target"),
        py::arg("shift") = std::nullopt, py::arg("name") = std::nullopt);
    pyQCircuit.def(
        "cCTRL",
        py::overload_cast<const cmat&, idx, const std::vector<idx>&,
                          std::optional<idx>, std::optional<std::string>>(
            &QCircuit::cCTRL),
        "Jointly applies the single qudit controlled gate U with "
        "classical control dit ctrl on the qudit indexes specified "
        "by target, i.e., cCTRL-U_{joint}",
        py::arg("U"), py::arg("ctrl_dit"), py::arg("target"),
        py::arg("shift") = std::nullopt, py::arg("name") = std::nullopt);
    pyQCircuit.def(
        "cCTRL",
        py::overload_cast<const cmat&, const std::vector<idx>&, idx,
                          std::optional<std::vector<idx>>,
                          std::optional<std::string>>(&QCircuit::cCTRL),
        "Applies the single qudit controlled gate U with multiple "
        "classical control dits listed in ctrl on the target qudit "
        "target, i.e., cCTRL-cCTRL-...-cCTRL-U",
        py::arg("U"), py::arg("ctrl_dits"), py::arg("target"),
        py::arg("shift") = std::nullopt, py::arg("name") = std::nullopt);
    pyQCircuit.def(
        "cCTRL",
        py::overload_cast<
            const cmat&, const std::vector<idx>&, const std::vector<idx>&,
            std::optional<std::vector<idx>>, std::optional<std::string>>(
            &QCircuit::cCTRL),
        "Jointly applies the multiple-qudit controlled gate U with "
        "multiple classical control dits listed in ctrl on the qudit "
        "indexes specified by target, i.e., "
        "cCTRL-cCTRL-...-cCTRL-U_{joint}",
        py::arg("U"), py::arg("ctrl_dits"), py::arg("target"),
        py::arg("shift") = std::nullopt, py::arg("name") = std::nullopt);
    pyQCircuit.def(
        "cCTRL_fan",
        py::overload_cast<
            const cmat&, const std::vector<idx>&, const std::vector<idx>&,
            std::optional<std::vector<idx>>, std::optional<std::string>>(
            &QCircuit::cCTRL_fan),
        "Applies the single qudit controlled gate U with multiple "
        "classical control dits listed in ctrl on every qudit "
        "listed in target, i.e., cCTRL-cCTRL-...-cCTRL-U-U-...-U",
        py::arg("U"), py::arg("ctrl_dits"), py::arg("target"),
        py::arg("shift") = std::nullopt, py::arg("name") = std::nullopt);
    pyQCircuit.def(
        "cCTRL_fan",
        py::overload_cast<const cmat&, idx, const std::vector<idx>&,
                          std::optional<idx>, std::optional<std::string>>(
            &QCircuit::cCTRL_fan),
        "Applies the single qudit controlled gate U with classical "
        "control dit ctrl on every qudit listed in target, "
        "i.e.,cCTRL-U-U-...-U",
        py::arg("U"), py::arg("ctrl_dit"), py::arg("target"),
        py::arg("shift") = std::nullopt, py::arg("name") = std::nullopt);
    pyQCircuit.def(
        "compress", &QCircuit::compress,
        "Removes all clean qudits and relabels the rest of the qudits "
        "accordingly",
        py::arg("compress_dits") = false);
    pyQCircuit.def(
        "CTRL",
        py::overload_cast<const cmat&, idx, idx, std::optional<idx>,
                          std::optional<std::string>>(&QCircuit::CTRL),
        "Applies the single qudit controlled gate U with control "
        "qudit ctrl and target qudit target, i.e., CTRL-U",
        py::arg("U"), py::arg("ctrl"), py::arg("target"),
        py::arg("shift") = std::nullopt, py::arg("name") = std::nullopt);
    pyQCircuit.def(
        "CTRL",
        py::overload_cast<const cmat&, idx, const std::vector<idx>&,
                          std::optional<idx>, std::optional<std::string>>(
            &QCircuit::CTRL),
        "Jointly applies the single qudit controlled gate U with "
        "control qudit ctrl on the qudit indexes specified by target, "
        "i.e., CTRL-U_{joint}",
        py::arg("U"), py::arg("ctrl"), py::arg("target"),
        py::arg("shift") = std::nullopt, py::arg("name") = std::nullopt);
    pyQCircuit.def(
        "CTRL",
        py::overload_cast<const cmat&, const std::vector<idx>&, idx,
                          std::optional<std::vector<idx>>,
                          std::optional<std::string>>(&QCircuit::CTRL),
        "Applies the multiple-qudit controlled gate U with multiple "
        "control qudits listed in ctrl on the target qudit specified "
        "by target, i.e., CTRL-CTRL-...-CTRL-U",
        py::arg("U"), py::arg("ctrl"), py::arg("target"),
        py::arg("shift") = std::nullopt, py::arg("name") = std::nullopt);
    pyQCircuit.def(
        "CTRL",
        py::overload_cast<
            const cmat&, const std::vector<idx>&, const std::vector<idx>&,
            std::optional<std::vector<idx>>, std::optional<std::string>>(
            &QCircuit::CTRL),
        "Jointly applies the multiple-qudit controlled gate U with "
        "multiple control qudits listed in ctrl on the qudit indexes "
        "specified by target, i.e., CTRL-CTRL-...-CTRL-U_{joint}",
        py::arg("U"), py::arg("ctrl"), py::arg("target"),
        py::arg("shift") = std::nullopt, py::arg("name") = std::nullopt);
    pyQCircuit.def(
        "CTRL_fan",
        py::overload_cast<
            const cmat&, const std::vector<idx>&, const std::vector<idx>&,
            std::optional<std::vector<idx>>, std::optional<std::string>>(
            &QCircuit::CTRL_fan),
        "Applies the single qudit controlled gate U with multiple "
        "control qudits listed in ctrl on every qudit listed in "
        "target, i.e., CTRL-CTRL-...-CTRL-U-U-...-U",
        py::arg("U"), py::arg("ctrl"), py::arg("target"),
        py::arg("shift") = std::nullopt, py::arg("name") = std::nullopt);
    pyQCircuit.def(
        "CTRL_fan",
        py::overload_cast<const cmat&, idx, const std::vector<idx>&,
                          std::optional<idx>, std::optional<std::string>>(
            &QCircuit::CTRL_fan),
        "Applies the single qudit controlled gate U with control "
        "qudit ctrl on every qudit listed in target, i.e., "
        "CTRL-U-U-...-U",
        py::arg("U"), py::arg("ctrl"), py::arg("target"),
        py::arg("shift") = std::nullopt, py::arg("name") = std::nullopt);
    pyQCircuit.def(
        "discard",
        py::overload_cast<idx, std::optional<std::string>>(&QCircuit::discard),
        "Discards single qudit by measuring it destructively in the "
        "Z-basis",
        py::arg("target"), py::arg("name") = "discard");
    pyQCircuit.def(
        "discard",
        py::overload_cast<const std::vector<idx>&, std::optional<std::string>>(
            &QCircuit::discard),
        "Discards multiple qudits by measuring them destructively in "
        "the Z-basis",
        py::arg("target"), py::arg("name") = "discard");
    pyQCircuit.def(
        "gate",
        py::overload_cast<const cmat&, idx, std::optional<std::string>>(
            &QCircuit::gate),
        "Applies the single qudit gate U on single qudit i", py::arg("U"),
        py::arg("i"), py::arg("name") = std::nullopt);
    pyQCircuit.def(
        "gate",
        py::overload_cast<const cmat&, idx, idx, std::optional<std::string>>(
            &QCircuit::gate),
        "Applies the two qudit gate U on qudits i and j", py::arg("U"),
        py::arg("i"), py::arg("j"), py::arg("name") = std::nullopt);
    pyQCircuit.def(
        "gate",
        py::overload_cast<const cmat&, idx, idx, idx,
                          std::optional<std::string>>(&QCircuit::gate),
        "Applies the three qudit gate U on qudits i, j and k", py::arg("U"),
        py::arg("i"), py::arg("j"), py::arg("k"),
        py::arg("name") = std::nullopt);
    pyQCircuit.def(
        "gate",
        py::overload_cast<const cmat&, const std::vector<idx>&,
                          std::optional<std::string>>(&QCircuit::gate),
        "Jointly applies the multiple-qudit gate U on the qudit "
        "indexes specified by target",
        py::arg("U"), py::arg("target"), py::arg("name") = std::nullopt);
    pyQCircuit.def(
        "gate_fan",
        py::overload_cast<const cmat&, const std::vector<idx>&,
                          std::optional<std::string>>(&QCircuit::gate_fan),
        "Applies the single qudit gate U on every qudit listed in "
        "target",
        py::arg("U"), py::arg("target"), py::arg("name") = std::nullopt);
    pyQCircuit.def("gate_fan",
                   py::overload_cast<const cmat&, std::optional<std::string>>(
                       &QCircuit::gate_fan),
                   "Applies the single qudit gate U on all of the remaining "
                   "non-measured qudits",
                   py::arg("U"), py::arg("name") = std::nullopt);
    pyQCircuit.def("get_clean_dits", &QCircuit::get_clean_dits,
                   "Vector of clean classical dits");
    pyQCircuit.def("get_clean_qudits", &QCircuit::get_clean_qudits,
                   "Vector of clean qudits");
    pyQCircuit.def("get_d", &QCircuit::get_d, "Qudit dimension");
    pyQCircuit.def("get_depth", &QCircuit::get_depth,
                   "Quantum circuit description total depth");
    pyQCircuit.def("get_dirty_dits", &QCircuit::get_dirty_dits,
                   "Vector of dirty classical dits");
    pyQCircuit.def("get_dirty_qudits", &QCircuit::get_dirty_qudits,
                   "Vector of dirty qudits");
    pyQCircuit.def("get_gate_count", &QCircuit::get_gate_count,
                   "(Total) Gate count", py::arg("U") = std::nullopt);
    pyQCircuit.def("get_gate_depth", &QCircuit::get_gate_depth,
                   "(Total) Gate depth", py::arg("U") = std::nullopt);
    pyQCircuit.def("get_measured_d", &QCircuit::get_measured_d,
                   "Vector of already measured qudit indexes");
    pyQCircuit.def(
        "get_measured_nd", &QCircuit::get_measured_nd,
        "Vector of already measured non-destructively qudit indexes");
    pyQCircuit.def("get_measurement_count", &QCircuit::get_measurement_count,
                   "(Total) Measurement count", py::arg("V") = std::nullopt);
    pyQCircuit.def("get_measurement_depth", &QCircuit::get_measurement_depth,
                   "(Total) Measurement depth", py::arg("V") = std::nullopt);
    pyQCircuit.def(
        "get_measurement_dits", &QCircuit::get_measurement_dits,
        "Vector of classical dits that were used to store results of "
        "measurements (either destructive or non-destructive)");
    pyQCircuit.def("get_name", &QCircuit::get_name, "Description name");
    pyQCircuit.def("get_nc", &QCircuit::get_nc, "Number of classical dits");
    pyQCircuit.def("get_non_measured_d", &QCircuit::get_non_measured_d,
                   "Non-measured qudit indexes");
    pyQCircuit.def("get_nop_count", &QCircuit::get_nop_count, "No-op count");
    pyQCircuit.def("get_nq", &QCircuit::get_nq, "Number of qudits");
    pyQCircuit.def("get_resources", &QCircuit::get_resources,
                   "Quantum circuit resources");
    pyQCircuit.def("validate_conditionals", &QCircuit::validate_conditionals,
                   "True if valid conditionals, false otherwise");
    pyQCircuit.def("get_step_count", &QCircuit::get_step_count,
                   "Total (gates + measurements) count");
    pyQCircuit.def("has_measurements", &QCircuit::has_measurements,
                   "True if the quantum circuit description contains any "
                   "measurements, false otherwise");
    pyQCircuit.def("has_conditionals", &QCircuit::has_conditionals,
                   "True if the quantum circuit description contains any "
                   "conditionals, false otherwise");
    pyQCircuit.def_static(
        "is_cCTRL", &QCircuit::is_cCTRL,
        "True if the gate step is a classically-controlled gate, "
        "false otherwise",
        py::arg("gate_step"));
    pyQCircuit.def("is_clean_dit", &QCircuit::is_clean_dit,
                   "Whether classical dit i in the quantum circuit description "
                   "was used before or not",
                   py::arg("i"));
    pyQCircuit.def(
        "is_clean_qudit", &QCircuit::is_clean_qudit,
        "Whether qudit i in the quantum circuit description was used "
        "before or not",
        py::arg("i"));
    pyQCircuit.def_static(
        "is_CTRL", &QCircuit::is_CTRL,
        "True if the gate step is a controlled gate, false otherwise",
        py::arg("gate_step"));
    pyQCircuit.def("is_measurement_dit", &QCircuit::is_measurement_dit,
                   "Whether classical dit i in the quantum circuit description "
                   "was used to store the result of a measurement (either "
                   "destructive or non-destructive)",
                   py::arg("i"));
    pyQCircuit.def(
        "kron", &QCircuit::kron,
        "Kronecker product with another quantum circuit description, "
        "in place",
        py::arg("qc"));
    pyQCircuit.def(
        "measure",
        py::overload_cast<idx, idx, bool, std::optional<std::string>>(
            &QCircuit::measure),
        "Z measurement of single qudit", py::arg("target"), py::arg("c_reg"),
        py::arg("destructive") = true, py::arg("name") = "mZ");
    pyQCircuit.def(
        "measure",
        py::overload_cast<const std::vector<idx>&, idx, bool,
                          std::optional<std::string>>(&QCircuit::measure),
        "Z measurement of multiple qudits", py::arg("target"),
        py::arg("c_reg") = 0, py::arg("destructive") = true,
        py::arg("name") = "mZ");
    pyQCircuit.def("measure_all", &QCircuit::measure_all,
                   "Z measurement of all qudits", py::arg("c_reg") = 0,
                   py::arg("destructive") = true, py::arg("name") = "mZ");
    pyQCircuit.def(
        "measureV",
        py::overload_cast<const cmat&, idx, idx, bool,
                          std::optional<std::string>>(&QCircuit::measureV),
        "Measurement of single qudit in the orthonormal basis "
        "specified by the columns of matrix V",
        py::arg("V"), py::arg("target"), py::arg("c_reg"),
        py::arg("destructive") = true, py::arg("name") = std::nullopt);
    pyQCircuit.def(
        "measureV",
        py::overload_cast<const cmat&, const std::vector<idx>&, idx, bool,
                          std::optional<std::string>>(&QCircuit::measureV),
        "Measurement of multiple qudits in the orthonormal basis "
        "specified by the columns of matrix V",
        py::arg("V"), py::arg("target"), py::arg("c_reg"),
        py::arg("destructive") = true, py::arg("name") = std::nullopt);
    pyQCircuit.def("nop", &QCircuit::nop, "No operation (no-op)");
    pyQCircuit.def(
        "post_select",
        py::overload_cast<idx, idx, idx, bool, std::optional<std::string>>(
            &QCircuit::post_select),
        "Z post-selection of single qudit", py::arg("target"),
        py::arg("ps_val"), py::arg("c_reg"), py::arg("destructive") = true,
        py::arg("name") = std::nullopt);
    pyQCircuit.def(
        "post_select",
        py::overload_cast<const std::vector<idx>&, const std::vector<idx>&, idx,
                          bool, std::optional<std::string>>(
            &QCircuit::post_select),
        "Z post-selection of multiple qudits", py::arg("target"),
        py::arg("ps_vals"), py::arg("c_reg"), py::arg("destructive") = true,
        py::arg("name") = std::nullopt);
    pyQCircuit.def(
        "post_selectV",
        py::overload_cast<const cmat&, idx, idx, idx, bool,
                          std::optional<std::string>>(&QCircuit::post_selectV),
        "Post-selection of single qudit in the orthonormal basis "
        "specified by the columns of matrix V",
        py::arg("V"), py::arg("target"), py::arg("ps_val"), py::arg("c_reg"),
        py::arg("destructive") = true, py::arg("name") = std::nullopt);
    pyQCircuit.def(
        "post_selectV",
        py::overload_cast<const cmat&, const std::vector<idx>&, idx, idx, bool,
                          std::optional<std::string>>(&QCircuit::post_selectV),
        "Post-selection of multiple qudits in the orthonormal basis "
        "specified by the columns of matrix V",
        py::arg("V"), py::arg("target"), py::arg("ps_val"), py::arg("c_reg"),
        py::arg("destructive") = true, py::arg("name") = std::nullopt);
    pyQCircuit.def(
        "QFT", py::overload_cast<const std::vector<idx>&, bool>(&QCircuit::QFT),
        "Applies the quantum Fourier transform on the qudit indexes "
        "specified by target",
        py::arg("target"), py::arg("swap") = true);
    pyQCircuit.def("QFT", py::overload_cast<bool>(&QCircuit::QFT),
                   "Applies the quantum Fourier transform on all of remaining "
                   "non-measured qudits",
                   py::arg("swap") = true);
    pyQCircuit.def("remove_clean_dit", &QCircuit::remove_clean_dit,
                   "Removes clean classical dit and relabels the rest of the "
                   "classical dits accordingly",
                   py::arg("target"));
    pyQCircuit.def("remove_clean_dits", &QCircuit::remove_clean_dits,
                   "Removes clean classical dits and relabels the rest of the "
                   "classical dits accordingly",
                   py::arg("target"));
    pyQCircuit.def("remove_clean_qudit", &QCircuit::remove_clean_qudit,
                   "Removes clean qudit and relabels the rest of the qudits "
                   "accordingly",
                   py::arg("target"));
    pyQCircuit.def("remove_clean_qudits", &QCircuit::remove_clean_qudits,
                   "Removes clean qudits and relabels the rest of the qudits "
                   "accordingly",
                   py::arg("target"));
    pyQCircuit.def("removes_qudits", &QCircuit::removes_qudits,
                   "True if the quantum circuit description contains any "
                   "measurements that remove qudits, false otherwise");
    pyQCircuit.def("replicate", &QCircuit::replicate,
                   "Replicates the quantum circuit description, in place",
                   py::arg("n"));
    pyQCircuit.def(
        "reset",
        py::overload_cast<const std::vector<idx>&, std::optional<std::string>>(
            &QCircuit::reset),
        "Reset multiple qudits", py::arg("target"), py::arg("name") = "reset");
    pyQCircuit.def(
        "reset",
        py::overload_cast<idx, std::optional<std::string>>(&QCircuit::reset),
        "Reset single qudit", py::arg("target"), py::arg("name") = "reset");
    pyQCircuit.def("set_name", &QCircuit::set_name, "Sets name",
                   py::arg("name"));
    pyQCircuit.def("TFQ", py::overload_cast<bool>(&QCircuit::TFQ),
                   "Applies the inverse quantum Fourier transform on all of "
                   "remaining non-measured qudits",
                   py::arg("swap") = true);
    pyQCircuit.def(
        "TFQ", py::overload_cast<const std::vector<idx>&, bool>(&QCircuit::TFQ),
        "Applies the inverse quantum Fourier transform on the qudit "
        "indexes specified by target",
        py::arg("target"), py::arg("swap") = true);
    pyQCircuit.def("to_JSON", &QCircuit::to_JSON,
                   "Displays the quantum circuit description in JSON format",
                   py::arg("enclosed_in_curly_brackets") = true);
    pyQCircuit.def("was_measured_d", &QCircuit::was_measured_d,
                   "Whether qudit i was already measured", py::arg("i"));
    pyQCircuit.def("was_measured_nd", &QCircuit::was_measured_nd,
                   "Whether qudit i was already measured non-destructively",
                   py::arg("i"));

    pyQCircuit.def(py::self == py::self);
    pyQCircuit.def(py::self != py::self);
    pyQCircuit.def("__repr__", [](const QCircuit& self) {
        std::ostringstream oss;
        oss << self;
        return oss.str();
    });
    pyQCircuit.def("__copy__",
                   [](const QCircuit& self) { return QCircuit(self); });
    pyQCircuit.def("__deepcopy__", [](const QCircuit& self, py::dict) {
        return QCircuit(self);
    });

    /* qpp::internal::QCircuitResources */
    auto pyQCircuitResources =
        py::class_<internal::QCircuitResources>(pyQCircuit, "Resources");
    pyQCircuitResources.def_readonly("nq", &internal::QCircuitResources::nq);
    pyQCircuitResources.def_readonly("nc", &internal::QCircuitResources::nc);
    pyQCircuitResources.def_readonly("d", &internal::QCircuitResources::d);
    pyQCircuitResources.def_readonly("name",
                                     &internal::QCircuitResources::name);
    pyQCircuitResources.def_readonly("step_count",
                                     &internal::QCircuitResources::step_count);
    pyQCircuitResources.def_readonly("gate_count",
                                     &internal::QCircuitResources::gate_count);
    pyQCircuitResources.def_readonly("gate_depth",
                                     &internal::QCircuitResources::gate_depth);
    pyQCircuitResources.def_readonly(
        "measurement_count", &internal::QCircuitResources::measurement_count);
    pyQCircuitResources.def_readonly(
        "measurement_depth", &internal::QCircuitResources::measurement_depth);
    pyQCircuitResources.def_readonly("total_depth",
                                     &internal::QCircuitResources::total_depth);
    pyQCircuitResources.def("__repr__",
                            [](const internal::QCircuitResources& self) {
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
