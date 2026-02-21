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
 * @file <pyqpp/instruments_bind.hpp>
 * @file <pyqpp/instruments_bind.hpp>
 * @brief Bindings for <qpp/instruments.hpp>
 * @brief Bindings for <qpp/instruments.hpp>
 */

#ifndef PYQPP_INSTRUMENTS_BIND_HPP_
#define PYQPP_INSTRUMENTS_BIND_HPP_

#include "pyqpp/pyqpp_common.hpp"

inline void init_instruments(py::module_& m) {
    using namespace qpp;

    m.def(
        "ip",
        [](const ket& phi, const ket& psi, const std::vector<idx>& subsys,
           const std::vector<idx>& dims) {
            return qpp::ip(phi, psi, subsys, dims);
        },
        "Generalized inner product. Computes the projection of the "
        "multi-partite state psi onto the subsystem phi.",
        py::arg("phi"), py::arg("psi"), py::arg("subsys"), py::arg("dims"));
    m.def(
        "ip",
        [](const ket& phi, const ket& psi, const std::vector<idx>& subsys,
           idx d) { return qpp::ip(phi, psi, subsys, d); },
        "Generalized inner product. Computes the projection of the "
        "multi-partite state psi onto the subsystem phi, assuming uniform "
        "subsystem dimensions.",
        py::arg("phi"), py::arg("psi"), py::arg("subsys"), py::arg("d") = 2);
    m.def(
        "measure",
        [](const cmat& A, const std::vector<cmat>& Ks) {
            return qpp::measure(A, Ks);
        },
        "Measures the state vector or density operator A using the set of "
        "Kraus operators Ks. "
        "Returns a tuple of: 1. Result of the measurement (index), 2. Vector "
        "of outcome probabilities, "
        "and 3. Vector of post-measurement normalized states.",
        py::arg("A"), py::arg("Ks"));
    m.def(
        "measure",
        [](const cmat& A, const cmat& U) { return qpp::measure(A, U); },
        "Measures the state vector or density matrix A in the orthonormal "
        "basis specified by the unitary matrix U. "
        "Returns a tuple of: 1. Result of the measurement, 2. Vector of "
        "outcome probabilities, "
        "and 3. Vector of post-measurement normalized states.",
        py::arg("A"), py::arg("U"));
    m.def(
        "measure",
        [](const cmat& A, const std::vector<cmat>& Ks,
           const std::vector<idx>& target, const std::vector<idx>& dims,
           bool destructive) {
            return qpp::measure(A, Ks, target, dims, destructive);
        },
        "Measures the part 'target' of the multi-partite state vector or "
        "density matrix A "
        "using the set of Kraus operators Ks. If 'destructive' is True, the "
        "measured "
        "subsystems are traced away.",
        py::arg("A"), py::arg("Ks"), py::arg("target"), py::arg("dims"),
        py::arg("destructive") = true);
    m.def(
        "measure",
        [](const cmat& A, const std::vector<cmat>& Ks,
           const std::vector<idx>& target, idx d, bool destructive) {
            return qpp::measure(A, Ks, target, d, destructive);
        },
        "Measures the part 'target' of the multi-partite state vector or "
        "density matrix A "
        "using the set of Kraus operators Ks, assuming uniform subsystem "
        "dimensions. "
        "If 'destructive' is True, the measured subsystems are traced away.",
        py::arg("A"), py::arg("Ks"), py::arg("target"), py::arg("d") = 2,
        py::arg("destructive") = true);
    m.def(
        "measure",
        [](const cmat& A, const cmat& V, const std::vector<idx>& target,
           const std::vector<idx>& dims, bool destructive) {
            return qpp::measure(A, V, target, dims, destructive);
        },
        "Measures the part 'target' of the multi-partite state vector or "
        "density matrix A "
        "in the orthonormal basis specified by the columns of matrix V. "
        "If 'destructive' is True, the measured subsystems are traced away.",
        py::arg("A"), py::arg("V"), py::arg("target"), py::arg("dims"),
        py::arg("destructive") = true);
    m.def(
        "measure",
        [](const cmat& A, const cmat& V, const std::vector<idx>& target, idx d,
           bool destructive) {
            return qpp::measure(A, V, target, d, destructive);
        },
        "Measures the part 'target' of the multi-partite state vector or "
        "density matrix A "
        "in the orthonormal basis specified by the columns of matrix V, "
        "assuming "
        "uniform subsystem dimensions. If 'destructive' is True, the measured "
        "subsystems are traced away.",
        py::arg("A"), py::arg("V"), py::arg("target"), py::arg("d") = 2,
        py::arg("destructive") = true);
    m.def(
        "measure_seq",
        [](const cmat& A, std::vector<idx> target, std::vector<idx> dims,
           bool destructive) {
            return qpp::measure_seq(A, target, dims, destructive);
        },
        "Sequentially measures the subsystems specified by 'target' in the "
        "computational basis. "
        "Returns a tuple of: 1. Vector of measurement outcomes, 2. Vector of "
        "probabilities for each outcome, "
        "and 3. The resulting post-measurement normalized state.",
        py::arg("A"), py::arg("target"), py::arg("dims"),
        py::arg("destructive") = true);
    m.def(
        "measure_seq",
        [](const cmat& A, const std::vector<idx>& target, idx d,
           bool destructive) {
            return qpp::measure_seq(A, target, d, destructive);
        },
        "Sequentially measures the subsystems specified by 'target' in the "
        "computational basis, "
        "assuming uniform subsystem dimensions. Returns a tuple of: 1. Vector "
        "of measurement "
        "outcomes, 2. Outcome probabilities, and 3. The final post-measurement "
        "normalized state.",
        py::arg("A"), py::arg("target"), py::arg("d") = 2,
        py::arg("destructive") = true);
    m.def(
        "sample",
        [](const cmat& A, const std::vector<idx>& target,
           const std::vector<idx>& dims) {
            return qpp::sample(A, target, dims);
        },
        "Samples from a quantum state in the computational basis (Z-basis)",
        py::arg("A"), py::arg("target"), py::arg("dims"));
    m.def(
        "sample",
        [](const cmat& A, const std::vector<idx>& target, idx d = 2) {
            return qpp::sample(A, target, d);
        },
        "Samples from a quantum state in the computational basis (Z-basis)",
        py::arg("A"), py::arg("target"), py::arg("d") = 2);
    m.def(
        "sample",
        [](idx num_samples, const cmat& A, const std::vector<idx>& target,
           const std::vector<idx>& dims) {
            std::map<std::string, idx> result;
            auto stats = qpp::sample(num_samples, A, target, dims);
            for (auto&& elem : stats) {
                std::stringstream ss;
                ss << qpp::disp(
                    elem.first,
                    IOManipContainerOpts{}.set_sep("").set_left("").set_right(
                        ""));
                result[ss.str()] = elem.second;
            }
            return result;
        },
        "Samples repeatedly from a quantum state in the computational basis "
        "(Z-basis)",
        py::arg("num_samples"), py::arg("A"), py::arg("target"),
        py::arg("dims"));
    m.def(
        "sample",
        [](idx num_samples, const cmat& A, const std::vector<idx>& target,
           idx d = 2) {
            std::map<std::string, idx> result;
            auto stats = qpp::sample(num_samples, A, target, d);
            for (auto&& elem : stats) {
                std::stringstream ss;
                ss << qpp::disp(
                    elem.first,
                    IOManipContainerOpts{}.set_sep("").set_left("").set_right(
                        ""));
                result[ss.str()] = elem.second;
            }
            return result;
        },
        "Samples repeatedly from a quantum state in the computational basis "
        "(Z-basis)",
        py::arg("num_samples"), py::arg("A"), py::arg("target"),
        py::arg("d") = 2);
    m.def(
        "reset",
        [](const cmat& A, const std::vector<idx>& target,
           const std::vector<idx>& dims) {
            return qpp::reset(A, target, dims);
        },
        "Resets qudits in the multi-partite state A by performing a "
        "non-destructive "
        "computational basis measurement on the 'target' qudits and shifting "
        "them "
        "back to the |0...0> state.",
        py::arg("A"), py::arg("target"), py::arg("dims"));
    m.def(
        "reset",
        [](const cmat& A, const std::vector<idx>& target, idx d) {
            return qpp::reset(A, target, d);
        },
        "Resets qudits in the multi-partite state A, assuming uniform "
        "subsystem "
        "dimensions. Performs a non-destructive computational basis "
        "measurement "
        "and shifts target qudits back to the |0...0> state.",
        py::arg("A"), py::arg("target"), py::arg("d") = 2);
    m.def(
        "discard",
        [](const cmat& A, const std::vector<idx>& target,
           const std::vector<idx>& dims) {
            return qpp::discard(A, target, dims);
        },
        "Discards qudits from the multi-partite state A by performing a "
        "destructive "
        "measurement and discarding the results. Returns the reduced density "
        "matrix.",
        py::arg("A"), py::arg("target"), py::arg("dims"));
    m.def(
        "discard",
        [](const cmat& A, const std::vector<idx>& target,
           const std::vector<idx>& dims) {
            return qpp::discard(A, target, dims);
        },
        "Discards qudits from the multi-partite state A by performing a "
        "destructive "
        "measurement and discarding the results. Returns the reduced density "
        "matrix.",
        py::arg("A"), py::arg("target"), py::arg("dims"));
}

#endif /* PYQPP_INSTRUMENTS_BIND_HPP_ */
