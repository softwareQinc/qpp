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

#ifndef PYQPP_CLASSES_QDUMMY_ENGINE_BIND_HPP_
#define PYQPP_CLASSES_QDUMMY_ENGINE_BIND_HPP_

#include "pyqpp/pyqpp_common.h"

/* qpp::QKetDummyEngine and qpp::QDensityDummyEngine instantiator */
template <typename T>
void declare_QDummyEngine(py::module& m) {
    using namespace qpp;
    static_assert(std::is_same_v<T, ket> || std::is_same_v<T, cmat>,
                  "The underlying type must be qpp::ket or qpp::cmat");

    std::string pyname;
    if constexpr (std::is_same_v<T, ket>) {
        pyname = "_QKetDummyEngine";
    } else {
        pyname = "_QDensityDummyEngine";
    }

    using DummyEngineType = QDummyEngine<T, QCircuit>;
    py::class_<DummyEngineType>(m, pyname.c_str())
        .def(py::init<const QCircuit&>(), py::keep_alive<1, 2>())
        .def("execute", py::overload_cast<idx>(&DummyEngineType::execute),
             "Executes the entire quantum circuit description",
             py::arg("reps") = 1)
        .def(
            "get_circuit",
            [](const DummyEngineType& qe) { return qe.get_circuit(); },
            "Underlying quantum circuit description")
        .def("get_state", &DummyEngineType::get_state,
             "Underlying quantum state")
        .def("set_state", &DummyEngineType::set_state,
             "Sets the underlying quantum state", py::arg("state"))
        .def("to_JSON", &DummyEngineType::to_JSON,
             "State of the engine in JSON format",
             py::arg("enclosed_in_curly_brackets") = true)
        .def("traits_get_name", &DummyEngineType::traits_get_name,
             "Engine name")
        .def("traits_is_noisy", &DummyEngineType::traits_is_noisy,
             "Noisy engine?")
        .def("traits_is_pure", &DummyEngineType::traits_is_pure,
             "Pure state engine?")

        .def("__repr__",
             [](const DummyEngineType& self) {
                 std::ostringstream oss;
                 oss << self;
                 return oss.str();
             })
        .def("__copy__",
             [](const DummyEngineType& self) { return DummyEngineType(self); })
        .def("__deepcopy__", [](const DummyEngineType& self, py::dict) {
            return DummyEngineType(self);
        });

    if constexpr (std::is_same_v<T, ket>) {
        m.def(
            "QKetDummyEngine",
            [](const QCircuit& qc) { return DummyEngineType(qc); },
            py::keep_alive<0, 1>());
    } else {
        m.def(
            "QDensityDummyEngine",
            [](const QCircuit& qc) { return DummyEngineType(qc); },
            py::keep_alive<0, 1>());
    }
}

inline void init_classes_qdummy_engine(py::module_& m) {
    using namespace qpp;

    /* qpp::QKetDummyEngine */
    declare_QDummyEngine<ket>(m);
    /* qpp::QDensityDummyEngine */
    declare_QDummyEngine<cmat>(m);
}

#endif /* PYQPP_CLASSES_QDUMMY_ENGINE_BIND_HPP_ */
