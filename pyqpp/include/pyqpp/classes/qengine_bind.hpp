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

#ifndef PYQPP_CLASSES_QENGINE_BIND_HPP_
#define PYQPP_CLASSES_QENGINE_BIND_HPP_

#include "pyqpp/pyqpp_common.h"

/* qpp::QEngineT instantiator */
template <typename T>
void declare_QEngineT(py::module& m) {
    using namespace qpp;
    static_assert(std::is_same_v<T, ket> || std::is_same_v<T, cmat>,
                  "The underlying type must be qpp::ket or qpp::cmat");

    std::string pyname;
    if constexpr (std::is_same_v<T, ket>) {
        pyname = "_QKetEngine";
    } else {
        pyname = "_QDensityEngine";
    }

    py::class_<QEngineT<T>>(m, pyname.c_str())
        .def(py::init<const QCircuit&>(), py::keep_alive<1, 2>())

        .def("execute", py::overload_cast<idx>(&QEngineT<T>::execute),
             "Executes the entire quantum circuit description",
             py::arg("reps") = 1)
        .def(
            "get_circuit",
            [](const QEngineT<T>& qe) { return qe.get_circuit(); },
            "Underlying quantum circuit description")
        .def("get_dit", &QEngineT<T>::get_dit,
             "Underlying classical dit at position i", py::arg("i"))
        .def("get_dits", &QEngineT<T>::get_dits, "Underlying classical dits")
        .def("get_measured_destructively",
             &QEngineT<T>::get_measured_destructively,
             "Vector of already destructively measured qudit indexes")
        .def("get_non_measured_destructively",
             &QEngineT<T>::get_non_measured_destructively,
             "Vector of qudit indexes that were not measured destructively")
        .def("get_probs", &QEngineT<T>::get_probs,
             "Underlying measurement outcome probabilities")
        .def("get_state", &QEngineT<T>::get_state, "Underlying quantum state")
        .def(
            "get_stats",
            [](const QEngineT<T>& qe) {
                std::map<std::string, idx> result;
                const auto& stats = qe.get_stats();
                for (auto&& elem : stats.data()) {
                    std::stringstream ss;
                    ss << disp(elem.first, IOManipContainerOpts{}
                                               .set_sep("")
                                               .set_left("")
                                               .set_right(""));
                    result[ss.str()] = elem.second;
                }
                return result;
            },
            "Measurement statistics for multiple runs")
        .def("post_select_ok", &QEngineT<T>::post_select_ok,
             "Successful post-selection")
        .def("reset", &QEngineT<T>::reset, "Resets the engine",
             py::arg("reset_stats") = true,
             py::arg("ensure_post_selection") = true)
        .def("reset_stats", &QEngineT<T>::reset_stats,
             "Resets the collected measurement statistics hash table")
        .def("set_dit", &QEngineT<T>::set_dit,
             "Sets the classical dit at position i", py::arg("i"),
             py::arg("value"))
        .def("set_dits", &QEngineT<T>::set_dits, "Set the classical dits",
             py::arg("dits"))
        .def("set_state", &QEngineT<T>::set_state,
             "Sets the underlying quantum state", py::arg("state"))
        .def("to_JSON", &QEngineT<T>::to_JSON,
             "State of the engine in JSON format",
             py::arg("enclosed_in_curly_brackets") = true)
        .def("was_measured_destructively",
             &QEngineT<T>::was_measured_destructively,
             "Whether qudit i was already measured destructively", py::arg("i"))

        .def("traits_get_name", &QEngineT<T>::traits_get_name, "Engine name")
        .def("traits_is_noisy", &QEngineT<T>::traits_is_noisy, "Noisy engine?")
        .def("traits_is_pure", &QEngineT<T>::traits_is_pure,
             "Pure state engine?")

        .def("__repr__",
             [](const QEngineT<T>& self) {
                 std::ostringstream oss;
                 oss << self;
                 return oss.str();
             })
        .def("__copy__",
             [](const QEngineT<T>& self) { return QEngineT<T>(self); })
        .def("__deepcopy__", [](const QEngineT<T>& self, py::dict) {
            return QEngineT<T>(self);
        });

    if constexpr (std::is_same_v<T, ket>) {
        m.def(
            "QKetEngine", [](const QCircuit& qc) { return QEngineT<T>(qc); },
            py::keep_alive<0, 1>());
        // backwards compatibility
        m.attr("QEngine") = m.attr("QKetEngine");
    } else {
        m.def(
            "QDensityEngine",
            [](const QCircuit& qc) { return QEngineT<T>(qc); },
            py::keep_alive<0, 1>());
    }
}
/* qpp::QEngine */
inline void init_classes_qengine(py::module_& m) {
    using namespace qpp;

    /* qpp::QEngineT instantiation, pure */
    declare_QEngineT<ket>(m);
    /* qpp::QEngineT instantiation, mixed */
    declare_QEngineT<cmat>(m);
}

#endif /* PYQPP_CLASSES_QENGINE_BIND_HPP_ */
