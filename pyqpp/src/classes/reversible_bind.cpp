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

#include "pyqpp_common.h"

/* qpp::Dynamic_bitset and qpp::Bit_circuit */
void init_classes_reversible(py::module_& m) {
    using namespace qpp;

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
}
