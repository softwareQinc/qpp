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
 * \file <pyqpp/operations_bind.hpp>
 * \brief Bindings for <qpp/operations.hpp>
 */

#ifndef PYQPP_OPERATIONS_BIND_HPP_
#define PYQPP_OPERATIONS_BIND_HPP_

#include "pyqpp/pyqpp_common.hpp"

inline void init_operations(py::module_& m) {
    using namespace qpp;

    m.def(
        "apply",
        [](const cmat& state, const cmat& A, const std::vector<idx>& target,
           const std::vector<idx>& dims) {
            return qpp::apply(state, A, target, dims);
        },
        "Applies the gate A to the part target of the multi-partite state "
        "vector or density matrix state",
        py::arg("state"), py::arg("A"), py::arg("target"), py::arg("dims"));
    m.def(
        "apply",
        [](const cmat& state, const cmat& A, const std::vector<idx>& target,
           idx d) { return qpp::apply(state, A, target, d); },
        "Applies the gate A to the part target of the multi-partite state "
        "vector or density matrix state",
        py::arg("state"), py::arg("A"), py::arg("target"), py::arg("d") = 2);
    m.def(
        "apply",
        [](const cmat& A, const std::vector<cmat>& Ks) {
            return qpp::apply(A, Ks);
        },
        "Applies the channel specified by the set of Kraus operators Ks to the "
        "density matrix A",
        py::arg("A"), py::arg("Ks"));
    m.def(
        "apply",
        [](const cmat& A, const std::vector<cmat>& Ks,
           const std::vector<idx>& target, const std::vector<idx>& dims) {
            return qpp::apply(A, Ks, target, dims);
        },
        "Applies the channel specified by the set of Kraus operators Ks to the "
        "part target of the multi-partite density matrix A",
        py::arg("A"), py::arg("Ks"), py::arg("target"), py::arg("dims"));
    m.def(
        "apply",
        [](const cmat& A, const std::vector<cmat>& Ks,
           const std::vector<idx>& target,
           idx d) { return qpp::apply(A, Ks, target, d); },
        "Applies the channel specified by the set of Kraus operators Ks to the "
        "part target of the multi-partite density matrix A",
        py::arg("A"), py::arg("Ks"), py::arg("target"), py::arg("d") = 2);
    m.def(
        "applyCTRL",
        [](const cmat& state, const cmat& A, const std::vector<idx>& ctrl,
           const std::vector<idx>& target, const std::vector<idx>& dims,
           std::optional<std::vector<idx>> shift) {
            return qpp::applyCTRL(state, A, ctrl, target, dims, shift);
        },
        "Applies the controlled-gate A to the part target of the multi-partite "
        "state vector or density matrix state",
        py::arg("state"), py::arg("A"), py::arg("ctrl"), py::arg("target"),
        py::arg("dims"), py::arg("shift") = py::none());
    m.def(
        "applyCTRL",
        [](const cmat& state, const cmat& A, const std::vector<idx>& ctrl,
           const std::vector<idx>& target, idx d,
           std::optional<std::vector<idx>> shift) {
            return qpp::applyCTRL(state, A, ctrl, target, d, shift);
        },
        "Applies the controlled-gate A to the part target of the multi-partite "
        "state vector or density matrix state",
        py::arg("state"), py::arg("A"), py::arg("ctrl"), py::arg("target"),
        py::arg("d") = 2, py::arg("shift") = py::none());
    m.def(
        "applyCTRL_fan",
        [](const cmat& state, const cmat& A, const std::vector<idx>& ctrl,
           const std::vector<idx>& target, const std::vector<idx>& dims,
           std::optional<std::vector<idx>> shift) {
            return qpp::applyCTRL_fan(state, A, ctrl, target, dims, shift);
        },
        "Applies the single qudit controlled-gate A with multiple control "
        "qudits listed in ctrl on every qudit listed in target",
        py::arg("state"), py::arg("A"), py::arg("ctrl"), py::arg("target"),
        py::arg("dims"), py::arg("shift") = py::none());
    m.def(
        "applyCTRL_fan",
        [](const cmat& state, const cmat& A, const std::vector<idx>& ctrl,
           const std::vector<idx>& target, idx d,
           std::optional<std::vector<idx>> shift) {
            return qpp::applyCTRL_fan(state, A, ctrl, target, d, shift);
        },
        "Applies the single qudit controlled-gate A with multiple control "
        "qudits listed in ctrl on every qudit listed in target",
        py::arg("state"), py::arg("A"), py::arg("ctrl"), py::arg("target"),
        py::arg("d") = 2, py::arg("shift") = py::none());
    m.def(
        "applyQFT",
        [](const cmat& A, const std::vector<idx>& target, idx d, bool swap) {
            return qpp::applyQFT(A, target, d, swap);
        },
        "Applies the qudit quantum Fourier transform to the part target of the "
        "multi-partite state vector or density matrix A",
        py::arg("A"), py::arg("target"), py::arg("d") = 2,
        py::arg("swap") = true);
    m.def(
        "applyTFQ",
        [](const cmat& A, const std::vector<idx>& target, idx d, bool swap) {
            return qpp::applyTFQ(A, target, d, swap);
        },
        "Applies the inverse (adjoint) qudit quantum Fourier transform to the "
        "part target of the multi-partite state vector or density matrix A",
        py::arg("A"), py::arg("target"), py::arg("d") = 2,
        py::arg("swap") = true);
    m.def(
        "QFT",
        [](const cmat& A, idx d, bool swap) { return qpp::QFT(A, d, swap); },
        "Qudit quantum Fourier transform", py::arg("A"), py::arg("d") = 2,
        py::arg("swap") = true);
    m.def(
        "TFQ",
        [](const cmat& A, idx d, bool swap) { return qpp::TFQ(A, d, swap); },
        "Inverse (adjoint) qudit quantum Fourier transform", py::arg("A"),
        py::arg("d") = 2, py::arg("swap") = true);
    m.def(
        "kraus2choi",
        [](const std::vector<cmat>& Ks) { return qpp::kraus2choi(Ks); },
        "Choi matrix. Constructs the Choi matrix of the channel specified by "
        "the set of Kraus operators Ks",
        py::arg("Ks"));
    m.def(
        "choi2kraus",
        [](const cmat& A, idx Din, idx Dout) {
            return qpp::choi2kraus(A, Din, Dout);
        },
        "Orthogonal Kraus operators from Choi matrix. Extracts a set of "
        "orthogonal Kraus operators from the Choi matrix A",
        py::arg("A"), py::arg("Din"), py::arg("Dout"));
    m.def(
        "choi2kraus", [](const cmat& A) { return qpp::choi2kraus(A); },
        "Orthogonal Kraus operators from Choi matrix. Extracts a set of "
        "orthogonal Kraus operators from the Choi matrix A, assuming square "
        "Kraus operators",
        py::arg("A"));
    m.def(
        "choi2super",
        [](const cmat& A, idx Din, idx Dout) {
            return qpp::choi2super(A, Din, Dout);
        },
        "Converts Choi matrix to superoperator matrix", py::arg("A"),
        py::arg("Din"), py::arg("Dout"));
    m.def(
        "choi2super", [](const cmat& A) { return qpp::choi2super(A); },
        "Converts Choi matrix to superoperator matrix, assuming square Kraus "
        "operators",
        py::arg("A"));
    m.def(
        "super2choi", [](const cmat& A) { return qpp::super2choi(A); },
        "Converts superoperator matrix to Choi matrix", py::arg("A"));
    m.def(
        "kraus2super",
        [](const std::vector<cmat>& Ks) { return qpp::kraus2super(Ks); },
        "Superoperator matrix. Constructs the superoperator matrix of the "
        "channel specified by the set of Kraus operators Ks",
        py::arg("Ks"));
    m.def(
        "super2kraus", [](const cmat& A) { return qpp::super2kraus(A); },
        "Orthogonal Kraus operators from superoperator matrix. Extracts a set "
        "of orthogonal Kraus operators from the superoperator matrix A",
        py::arg("A"));
    m.def(
        "ptrace1",
        [](const cmat& A, const std::vector<idx>& dims) {
            return qpp::ptrace1(A, dims);
        },
        "Partial trace over the first subsystem of a bi-partite state vector "
        "or density matrix",
        py::arg("A"), py::arg("dims"));
    m.def(
        "ptrace1", [](const cmat& A, idx d) { return qpp::ptrace1(A, d); },
        "Partial trace over the first subsystem of a bi-partite state vector "
        "or density matrix, assuming uniform subsystem dimensions",
        py::arg("A"), py::arg("d") = 2);
    m.def(
        "ptrace2",
        [](const cmat& A, const std::vector<idx>& dims) {
            return qpp::ptrace2(A, dims);
        },
        "Partial trace over the second subsystem of a bi-partite state vector "
        "or density matrix",
        py::arg("A"), py::arg("dims"));
    m.def(
        "ptrace2", [](const cmat& A, idx d) { return qpp::ptrace2(A, d); },
        "Partial trace over the second subsystem of a bi-partite state vector "
        "or density matrix, assuming uniform subsystem dimensions",
        py::arg("A"), py::arg("d") = 2);
    m.def(
        "ptrace",
        [](const cmat& A, const std::vector<idx>& target,
           const std::vector<idx>& dims) {
            return qpp::ptrace(A, target, dims);
        },
        "Partial trace of the multi-partite state vector or density matrix "
        "over the list target of subsystems",
        py::arg("A"), py::arg("target"), py::arg("dims"));
    m.def(
        "ptrace",
        [](const cmat& A, const std::vector<idx>& target, idx d) {
            return qpp::ptrace(A, target, d);
        },
        "Partial trace of the multi-partite state vector or density matrix "
        "over the list target of subsystems, assuming uniform subsystem "
        "dimensions",
        py::arg("A"), py::arg("target"), py::arg("d") = 2);
    m.def(
        "ptranspose",
        [](const cmat& A, const std::vector<idx>& target,
           const std::vector<idx>& dims) {
            return qpp::ptranspose(A, target, dims);
        },
        "Partial transpose of the multi-partite state vector or density matrix "
        "over the list target of subsystems",
        py::arg("A"), py::arg("target"), py::arg("dims"));
    m.def(
        "ptranspose",
        [](const cmat& A, const std::vector<idx>& target, idx d) {
            return qpp::ptranspose(A, target, d);
        },
        "Partial transpose of the multi-partite state vector or density matrix "
        "over the list target of subsystems, assuming uniform subsystem "
        "dimensions",
        py::arg("A"), py::arg("target"), py::arg("d") = 2);
    m.def(
        "syspermute",
        [](const cmat& A, const std::vector<idx>& perm,
           const std::vector<idx>& dims) {
            return qpp::syspermute(A, perm, dims);
        },
        "Subsystem permutation. Permutes the subsystems of a state vector or "
        "density matrix. The subsystem perm[i] is moved to the location i.",
        py::arg("A"), py::arg("perm"), py::arg("dims"));
    m.def(
        "syspermute",
        [](const cmat& A, const std::vector<idx>& perm, idx d) {
            return qpp::syspermute(A, perm, d);
        },
        "Subsystem permutation. Permutes the subsystems of a state vector or "
        "density matrix, assuming uniform subsystem dimensions.",
        py::arg("A"), py::arg("perm"), py::arg("d") = 2);
    m.def(
        "qRAM",
        [](const ket& psi, const qram& data, idx DqRAM) {
            return qpp::qRAM(psi, data, DqRAM);
        },
        "Quantumly-accessible Random Access Memory (qRAM) over classical data. "
        "Implements the mapping sum_j alpha_j |j> -> sum_j alpha_j |j>|m_j>.",
        py::arg("psi"), py::arg("data"), py::arg("DqRAM"));
    m.def(
        "qRAM",
        [](const ket& psi, const qram& data) { return qpp::qRAM(psi, data); },
        "Quantumly-accessible Random Access Memory (qRAM) over classical data. "
        "The qRAM subsystem dimension is automatically set to 1 + maximum "
        "value stored in the data.",
        py::arg("psi"), py::arg("data"));
}

#endif /* PYQPP_OPERATIONS_BIND_HPP_ */
