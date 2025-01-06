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

#ifndef PYQPP_CLASSES_REVERSIBLE_BIND_HPP_
#define PYQPP_CLASSES_REVERSIBLE_BIND_HPP_

#include "pyqpp/pyqpp_common.hpp"

/* qpp::Dynamic_bitset and qpp::Bit_circuit */
inline void init_classes_reversible(py::module_& m) {
    using namespace qpp;

    /* qpp::Dynamic_bitset */
    auto pyDynamic_bitset = py::class_<Dynamic_bitset>(m, "Dynamic_bitset");

    pyDynamic_bitset.def(py::init<idx>(), py::arg("n"));
    pyDynamic_bitset.def(py::init<std::string, char, char>(), py::arg("str"),
                         py::arg("zero") = '0', py::arg("one") = '1');

    pyDynamic_bitset.def("all", &Dynamic_bitset::all,
                         "True if all of the bits are set");
    pyDynamic_bitset.def("any", &Dynamic_bitset::any,
                         "True if any of the bits is set");
    pyDynamic_bitset.def(
        "count", &Dynamic_bitset::count,
        "Number of bits set to one in the bitset (Hamming weight)");
    pyDynamic_bitset.def("flip", py::overload_cast<>(&Dynamic_bitset::flip),
                         "Flips all bits");
    pyDynamic_bitset.def("flip", py::overload_cast<idx>(&Dynamic_bitset::flip),
                         "Flips the bit at position pos", py::arg("pos"));
    pyDynamic_bitset.def("get", &Dynamic_bitset::get,
                         "The value of the bit at position pos",
                         py::arg("pos"));
    pyDynamic_bitset.def("none", &Dynamic_bitset::none,
                         "True if none of the bits are set");
    pyDynamic_bitset.def(
        "rand", py::overload_cast<realT>(&Dynamic_bitset::rand),
        "Sets all bits according to a Bernoulli(p) distribution",
        py::arg("p") = 0.5);
    pyDynamic_bitset.def(
        "rand", py::overload_cast<idx, realT>(&Dynamic_bitset::rand),
        "Sets the bit at position pos according to a Bernoulli(p) "
        "distribution",
        py::arg("pos"), py::arg("p") = 0.5);
    pyDynamic_bitset.def("reset", py::overload_cast<>(&Dynamic_bitset::reset),
                         "Sets all bits to false");
    pyDynamic_bitset.def(
        "reset", py::overload_cast<idx>(&Dynamic_bitset::reset),
        "Sets the bit at position pos to false", py::arg("pos"));
    pyDynamic_bitset.def("set", py::overload_cast<>(&Dynamic_bitset::set),
                         "Sets all bits to true");
    pyDynamic_bitset.def("set",
                         py::overload_cast<idx, bool>(&Dynamic_bitset::set),
                         "Sets the bit at position pos", py::arg("pos"),
                         py::arg("value") = true);
    pyDynamic_bitset.def("size", &Dynamic_bitset::size,
                         "Number of bits stored in the bitset");
    pyDynamic_bitset.def(
        "storage_size", &Dynamic_bitset::storage_size,
        "Size of the underlying storage space (in units of "
        "qpp::Dynamic_bitset::value_type, unsigned int by default)");
    pyDynamic_bitset.def("to_string", &Dynamic_bitset::to_string,
                         "String representation", py::arg("zero") = '0',
                         py::arg("one") = '1');

    pyDynamic_bitset.def(py::self == py::self);
    pyDynamic_bitset.def(py::self != py::self);
    pyDynamic_bitset.def("__copy__", [](const Dynamic_bitset& self) {
        return Dynamic_bitset(self);
    });
    pyDynamic_bitset.def("__deepcopy__",
                         [](const Dynamic_bitset& self, py::dict) {
                             return Dynamic_bitset(self);
                         });
    pyDynamic_bitset.def("__repr__", [](const Dynamic_bitset& self) {
        std::ostringstream oss;
        oss << self;
        return oss.str();
    });
    pyDynamic_bitset.def(
        "__sub__", &Dynamic_bitset::operator-,
        "Number of places the two bitsets differ (Hamming distance)");

    /* qpp::Bit_circuit */
    auto pyBit_circuit =
        py::class_<Bit_circuit, Dynamic_bitset>(m, "Bit_circuit");

    pyBit_circuit.def(py::init<idx>(), py::arg("n"));
    pyBit_circuit.def(py::init<std::string, char, char>(), py::arg("str"),
                      py::arg("zero") = '0', py::arg("one") = '1');
    pyBit_circuit.def(py::init<const Dynamic_bitset&>(),
                      py::keep_alive<1, 2>());

    pyBit_circuit.def("CNOT", &Bit_circuit::CNOT, "Controlled-NOT gate",
                      py::arg("ctrl"), py::arg("target"));
    pyBit_circuit.def("FRED", &Bit_circuit::FRED,
                      " Fredkin gate (Controlled-SWAP)", py::arg("i"),
                      py::arg("j"), py::arg("k"));
    pyBit_circuit.def(
        "get_gate_count", &Bit_circuit::get_gate_count,
        "(Total) Bit circuit gate count. Possible names are NOT (X), "
        "CNOT, SWAP, TOF, FRED",
        py::arg("name") = std::nullopt);
    pyBit_circuit.def(
        "get_gate_depth", &Bit_circuit::get_gate_depth,
        "(Total) Bit circuit gate depth. Possible names are NOT (X), "
        "CNOT, SWAP, TOF, FRED",
        py::arg("name") = std::nullopt);
    pyBit_circuit.def("NOT", &Bit_circuit::NOT, "NOT gate (bit flip)",
                      py::arg("i"));
    pyBit_circuit.def("reset", &Bit_circuit::reset,
                      "Resets the circuit to all-zero, clears all gates");
    pyBit_circuit.def("SWAP", &Bit_circuit::SWAP, "Swap gate", py::arg("i"),
                      py::arg("j"));
    pyBit_circuit.def("to_JSON", &Bit_circuit::to_JSON,
                      "Displays the bit circuit in JSON format",
                      py::arg("enclosed_in_curly_brackets") = true);
    pyBit_circuit.def("to_string", &Bit_circuit::to_string,
                      "String representation", py::arg("zero") = '0',
                      py::arg("one") = '1');
    pyBit_circuit.def("TOF", &Bit_circuit::TOF, "Toffoli gate", py::arg("i"),
                      py::arg("j"), py::arg("k"));
    pyBit_circuit.def("X", &Bit_circuit::X, "NOT gate (bit flip)",
                      py::arg("i"));

    pyBit_circuit.def(py::self == py::self);
    pyBit_circuit.def(py::self != py::self);
    pyBit_circuit.def(
        "__copy__", [](const Bit_circuit& self) { return Bit_circuit(self); });
    pyBit_circuit.def("__deepcopy__", [](const Bit_circuit& self, py::dict) {
        return Bit_circuit(self);
    });
    pyBit_circuit.def("__repr__", [](const Bit_circuit& self) {
        std::ostringstream oss;
        oss << self;
        return oss.str();
    });
}

#endif /* PYQPP_CLASSES_REVERSIBLE_BIND_HPP_ */
