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
 * \file <pyqpp/classes/reversible_bind.hpp>
 * \brief Bindings for <qpp/classes/reversible.hpp>
 */

#ifndef PYQPP_CLASSES_REVERSIBLE_BIND_HPP_
#define PYQPP_CLASSES_REVERSIBLE_BIND_HPP_

#include "pyqpp/pyqpp_common.hpp"

/* qpp::DynamicBitset and qpp::BitCircuit */
inline void init_classes_reversible(py::module_& m) {
    using namespace qpp;

    /* qpp::DynamicBitset */
    auto pyDynamicBitset = py::class_<DynamicBitset>(m, "DynamicBitset");

    pyDynamicBitset.def(py::init<idx>(), py::arg("n"));
    pyDynamicBitset.def(py::init<std::string, char, char>(), py::arg("str"),
                         py::arg("zero") = '0', py::arg("one") = '1');

    pyDynamicBitset.def("all", &DynamicBitset::all,
                         "True if all of the bits are set");
    pyDynamicBitset.def("any", &DynamicBitset::any,
                         "True if any of the bits is set");
    pyDynamicBitset.def(
        "count", &DynamicBitset::count,
        "Number of bits set to one in the bitset (Hamming weight)");
    pyDynamicBitset.def("flip", py::overload_cast<>(&DynamicBitset::flip),
                         "Flips all bits");
    pyDynamicBitset.def("flip", py::overload_cast<idx>(&DynamicBitset::flip),
                         "Flips the bit at position pos", py::arg("pos"));
    pyDynamicBitset.def("get", &DynamicBitset::get,
                         "The value of the bit at position pos",
                         py::arg("pos"));
    pyDynamicBitset.def("none", &DynamicBitset::none,
                         "True if none of the bits are set");
    pyDynamicBitset.def(
        "rand", py::overload_cast<realT>(&DynamicBitset::rand),
        "Sets all bits according to a Bernoulli(p) distribution",
        py::arg("p") = 0.5);
    pyDynamicBitset.def(
        "rand", py::overload_cast<idx, realT>(&DynamicBitset::rand),
        "Sets the bit at position pos according to a Bernoulli(p) "
        "distribution",
        py::arg("pos"), py::arg("p") = 0.5);
    pyDynamicBitset.def("reset", py::overload_cast<>(&DynamicBitset::reset),
                         "Sets all bits to false");
    pyDynamicBitset.def(
        "reset", py::overload_cast<idx>(&DynamicBitset::reset),
        "Sets the bit at position pos to false", py::arg("pos"));
    pyDynamicBitset.def("set", py::overload_cast<>(&DynamicBitset::set),
                         "Sets all bits to true");
    pyDynamicBitset.def("set",
                         py::overload_cast<idx, bool>(&DynamicBitset::set),
                         "Sets the bit at position pos", py::arg("pos"),
                         py::arg("value") = true);
    pyDynamicBitset.def("size", &DynamicBitset::size,
                         "Number of bits stored in the bitset");
    pyDynamicBitset.def(
        "storage_size", &DynamicBitset::storage_size,
        "Size of the underlying storage space (in units of "
        "qpp::DynamicBitset::value_type, unsigned int by default)");
    pyDynamicBitset.def("to_string", &DynamicBitset::to_string,
                         "String representation", py::arg("zero") = '0',
                         py::arg("one") = '1');

    pyDynamicBitset.def(py::self == py::self);
    pyDynamicBitset.def(py::self != py::self);
    pyDynamicBitset.def("__copy__", [](const DynamicBitset& self) {
        return DynamicBitset(self);
    });
    pyDynamicBitset.def("__deepcopy__",
                         [](const DynamicBitset& self, py::dict) {
                             return DynamicBitset(self);
                         });
    pyDynamicBitset.def("__repr__", [](const DynamicBitset& self) {
        std::ostringstream oss;
        oss << self;
        return oss.str();
    });
    pyDynamicBitset.def(
        "__sub__", &DynamicBitset::operator-,
        "Number of places the two bitsets differ (Hamming distance)");

    /* qpp::BitCircuit */
    auto pyBitCircuit =
        py::class_<BitCircuit, DynamicBitset>(m, "BitCircuit");

    pyBitCircuit.def(py::init<idx>(), py::arg("n"));
    pyBitCircuit.def(py::init<std::string, char, char>(), py::arg("str"),
                      py::arg("zero") = '0', py::arg("one") = '1');
    pyBitCircuit.def(py::init<const DynamicBitset&>(),
                      py::keep_alive<1, 2>());

    pyBitCircuit.def("CNOT", &BitCircuit::CNOT, "Controlled-NOT gate",
                      py::arg("ctrl"), py::arg("target"));
    pyBitCircuit.def("FRED", &BitCircuit::FRED,
                      " Fredkin gate (Controlled-SWAP)", py::arg("i"),
                      py::arg("j"), py::arg("k"));
    pyBitCircuit.def(
        "get_gate_count", &BitCircuit::get_gate_count,
        "(Total) Bit circuit gate count. Possible names are NOT (X), "
        "CNOT, SWAP, TOF, FRED",
        py::arg("name") = std::nullopt);
    pyBitCircuit.def(
        "get_gate_depth", &BitCircuit::get_gate_depth,
        "(Total) Bit circuit gate depth. Possible names are NOT (X), "
        "CNOT, SWAP, TOF, FRED",
        py::arg("name") = std::nullopt);
    pyBitCircuit.def("NOT", &BitCircuit::NOT, "NOT gate (bit flip)",
                      py::arg("i"));
    pyBitCircuit.def("reset", &BitCircuit::reset,
                      "Resets the circuit to all-zero, clears all gates");
    pyBitCircuit.def("SWAP", &BitCircuit::SWAP, "Swap gate", py::arg("i"),
                      py::arg("j"));
    pyBitCircuit.def("to_JSON", &BitCircuit::to_JSON,
                      "Displays the bit circuit in JSON format",
                      py::arg("enclosed_in_curly_brackets") = true);
    pyBitCircuit.def("to_string", &BitCircuit::to_string,
                      "String representation", py::arg("zero") = '0',
                      py::arg("one") = '1');
    pyBitCircuit.def("TOF", &BitCircuit::TOF, "Toffoli gate", py::arg("i"),
                      py::arg("j"), py::arg("k"));
    pyBitCircuit.def("X", &BitCircuit::X, "NOT gate (bit flip)",
                      py::arg("i"));

    pyBitCircuit.def(py::self == py::self);
    pyBitCircuit.def(py::self != py::self);
    pyBitCircuit.def(
        "__copy__", [](const BitCircuit& self) { return BitCircuit(self); });
    pyBitCircuit.def("__deepcopy__", [](const BitCircuit& self, py::dict) {
        return BitCircuit(self);
    });
    pyBitCircuit.def("__repr__", [](const BitCircuit& self) {
        std::ostringstream oss;
        oss << self;
        return oss.str();
    });
}

#endif /* PYQPP_CLASSES_REVERSIBLE_BIND_HPP_ */
