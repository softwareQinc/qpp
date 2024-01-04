/*
 * This file is part of pyqpp.
 *
 * Copyright (c) 2019 - 2024 softwareQ Inc. All rights reserved.
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

#ifndef PYQPP_CLASSES_CIRCUITS_ENGINES_BIND_HPP_
#define PYQPP_CLASSES_CIRCUITS_ENGINES_BIND_HPP_

/* qpp::QNoisyEngine instantiator */
template <typename NoiseModel, typename... CtorTypeList>
void declare_noisy_engine(py::module& m, const std::string& type) {
    using namespace qpp;

    py::class_<NoiseModel>(m, type.c_str())
        .def(py::init<CtorTypeList...>())

        .def("get_d", &NoiseModel::get_d, "Qudit dimension")
        .def("get_Ks", &NoiseModel::get_Ks, "Vector of noise operators")
        .def("get_last_idx", &NoiseModel::get_last_idx,
             "Index of the last occurring noise element")
        .def("get_last_K", &NoiseModel::get_last_K,
             "Last occurring noise element")
        .def("get_last_p", &NoiseModel::get_last_p,
             "Probability of the last occurring noise element")
        .def("get_probs", &NoiseModel::get_probs,
             "Vector of probabilities corresponding to each noise operator");

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

/* qpp::QEngine */
inline void init_classes_circuits_engines(py::module_& m) {
    using namespace qpp;

    py::class_<QEngine>(m, "QEngine")
        .def(py::init<const QCircuit&>(), py::keep_alive<1, 2>())

        .def("execute", py::overload_cast<idx, bool>(&QEngine::execute),
             "Executes the entire quantum circuit description",
             py::arg("reps") = 1, py::arg("try_sampling") = true)
        .def(
            "get_circuit", [](const QEngine& qe) { return qe.get_circuit(); },
            "Underlying quantum circuit description")
        .def("get_dit", &QEngine::get_dit,
             "Underlying classical dit at position i", py::arg("i"))
        .def("get_dits", &QEngine::get_dits, "Underlying classical dits")
        .def("get_measured", &QEngine::get_measured,
             "Vector of already measured qudit indexes")
        .def("get_non_measured", &QEngine::get_non_measured,
             "Non-measured qudit indexes")
        .def("get_probs", &QEngine::get_probs,
             "Underlying measurement outcome probabilities")
        .def("get_psi", &QEngine::get_psi, "Underlying quantum state")
        .def(
            "get_stats",
            [](const QEngine& qe) {
                std::map<std::string, idx> result;
                const auto& stats = qe.get_stats();
                for (auto&& elem : stats.data()) {
                    std::stringstream ss;
                    ss << qpp::disp(elem.first, IOManipContainerOpts{}
                                                    .set_sep("")
                                                    .set_left("")
                                                    .set_right(""));
                    result[ss.str()] = elem.second;
                }
                return result;
            },
            "Measurement statistics for multiple runs")
        .def("is_noisy", &QEngine::is_noisy, "Whether the engine is noisy")
        .def("reset", &QEngine::reset, "Resets the engine",
             py::arg("reset_stats") = true)
        .def("reset_stats", &QEngine::reset_stats,
             "Resets the collected measurement statistics hash table")
        .def("set_dit", &QEngine::set_dit,
             "Sets the classical dit at position i", py::arg("i"),
             py::arg("value"))
        .def("set_dits", &QEngine::set_dits, "Set the classical dits",
             py::arg("dits"))
        .def("set_psi", &QEngine::set_psi, "Sets the underlying quantum state",
             py::arg("psi"))
        .def("to_JSON", &QEngine::to_JSON, "State of the engine in JSON format",
             py::arg("enclosed_in_curly_brackets") = true)
        .def("was_measured", &QEngine::was_measured,
             "Whether qudit i was already measured destructively", py::arg("i"))

        .def("__repr__",
             [](const QEngine& self) {
                 std::ostringstream oss;
                 oss << self;
                 return oss.str();
             })
        .def("__copy__", [](const QEngine& self) { return QEngine(self); })
        .def("__deepcopy__",
             [](const QEngine& self, py::dict) { return QEngine(self); });

    /* qpp::QNoisyEngine instantiations with different noise models */
    declare_noisy_engine<qpp::QubitBitFlipNoise, realT>(m, "QubitBitFlipNoise");
    declare_noisy_engine<qpp::QubitBitPhaseFlipNoise, realT>(
        m, "QubitBitPhaseFlipNoise");
    declare_noisy_engine<qpp::QubitDepolarizingNoise, realT>(
        m, "QubitDepolarizingNoise");
    declare_noisy_engine<qpp::QubitPhaseFlipNoise, realT>(
        m, "QubitPhaseFlipNoise");
    declare_noisy_engine<qpp::QuditDepolarizingNoise, realT, idx>(
        m, "QuditDepolarizingNoise");
}

#endif /* PYQPP_CLASSES_CIRCUITS_ENGINES_BIND_HPP_ */
