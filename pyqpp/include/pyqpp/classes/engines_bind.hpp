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

#ifndef PYQPP_CLASSES_ENGINES_BIND_HPP_
#define PYQPP_CLASSES_ENGINES_BIND_HPP_

#include "pyqpp/pyqpp_common.h"

/* qpp::QEngineT instantiator */
template <typename T>
void declare_ideal_engine(py::module& m) {
    using namespace qpp;
    static_assert(std::is_same_v<T, ket> || std::is_same_v<T, cmat>,
                  "The underlying type must be qpp::ket or qpp::cmat");

    std::string pyname;
    if constexpr (std::is_same_v<T, ket>) {
        pyname = "_QKetEngine";
    } else {
        pyname = "_QDensityEngine";
    }

    py::class_<qpp::QEngineT<T>>(m, pyname.c_str())
        .def(py::init<const QCircuit&>(), py::keep_alive<1, 2>())

        .def("execute", py::overload_cast<idx, bool>(&QEngineT<T>::execute),
             "Executes the entire quantum circuit description",
             py::arg("reps") = 1, py::arg("try_sampling") = true)
        .def(
            "get_circuit",
            [](const QEngineT<T>& qe) { return qe.get_circuit(); },
            "Underlying quantum circuit description")
        .def("get_dit", &QEngineT<T>::get_dit,
             "Underlying classical dit at position i", py::arg("i"))
        .def("get_dits", &QEngineT<T>::get_dits, "Underlying classical dits")
        .def("get_measured", &QEngineT<T>::get_measured,
             "Vector of already measured qudit indexes")
        .def("get_non_measured", &QEngineT<T>::get_non_measured,
             "Non-measured qudit indexes")
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
                    ss << qpp::disp(elem.first, IOManipContainerOpts{}
                                                    .set_sep("")
                                                    .set_left("")
                                                    .set_right(""));
                    result[ss.str()] = elem.second;
                }
                return result;
            },
            "Measurement statistics for multiple runs")
        .def("reset", &QEngineT<T>::reset, "Resets the engine",
             py::arg("reset_stats") = true)
        .def("reset_stats", &QEngineT<T>::reset_stats,
             "Resets the collected measurement statistics hash table")
        .def("set_dit", &QEngineT<T>::set_dit,
             "Sets the classical dit at position i", py::arg("i"),
             py::arg("value"))
        .def("set_dits", &QEngineT<T>::set_dits, "Set the classical dits",
             py::arg("dits"))
        .def("set_state", &QEngineT<T>::set_state,
             "Sets the underlying quantum state", py::arg("psi"))
        .def("to_JSON", &QEngineT<T>::to_JSON,
             "State of the engine in JSON format",
             py::arg("enclosed_in_curly_brackets") = true)
        .def("was_measured", &QEngineT<T>::was_measured,
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
            "QKetEngine",
            [](const QCircuit& qc) { return qpp::QEngineT<T>(qc); },
            py::keep_alive<0, 1>());
        // backwards compatibility
        m.attr("QEngine") = m.attr("QKetEngine");
    } else {
        m.def(
            "QDensityEngine",
            [](const QCircuit& qc) { return qpp::QEngineT<T>(qc); },
            py::keep_alive<0, 1>());
    }
}

/** Noise models instantiators */
template <typename NoiseModel, typename... CtorTypeList>
void declare_noise_model(py::module& m, const std::string& type) {
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
}

/* qpp::QNoisyEngineT instantiator */
template <typename T, typename NoiseModel>
void declare_noisy_engine(py::module& m, const std::string& type) {
    using namespace qpp;
    static_assert(std::is_same_v<T, ket> || std::is_same_v<T, cmat>,
                  "The underlying type must be qpp::ket or qpp::cmat");

    std::string pyname;
    if constexpr (std::is_same_v<T, ket>) {
        pyname = "_QKetNoisyEngine_" + type;
    } else {
        pyname = "_QDensityNoisyEngine_" + type;
    }
    py::class_<qpp::QNoisyEngineT<T, NoiseModel>, QEngineT<T>>(m,
                                                               pyname.c_str())
        .def(py::init<const QCircuit&, const NoiseModel&>(),
             py::keep_alive<1, 2>())

        .def(
            "get_noise_results",
            &qpp::QNoisyEngineT<T, NoiseModel>::get_noise_results,
            "Vector of noise results obtained before every step in the circuit")

        .def("__copy__",
             [](const qpp::QNoisyEngineT<T, NoiseModel>& self) {
                 return qpp::QNoisyEngineT<T, NoiseModel>(self);
             })
        .def("__deepcopy__",
             [](const qpp::QNoisyEngineT<T, NoiseModel>& self, py::dict) {
                 return qpp::QNoisyEngineT<T, NoiseModel>(self);
             });

    if constexpr (std::is_same_v<T, ket>) {
        m.def(
            "QKetNoisyEngine",
            [](const QCircuit& qc, const NoiseModel& nm) {
                return qpp::QNoisyEngineT<T, NoiseModel>(qc, nm);
            },
            py::keep_alive<0, 1>());
        // backwards compatibility
        m.attr("QNoisyEngine") = m.attr("QKetNoisyEngine");
    } else {
        m.def(
            "QDensityNoisyEngine",
            [](const QCircuit& qc, const NoiseModel& nm) {
                return qpp::QNoisyEngineT<T, NoiseModel>(qc, nm);
            },
            py::keep_alive<0, 1>());
    }
}

/* qpp::QEngine */
inline void init_classes_engines(py::module_& m) {
    using namespace qpp;

    /* qpp::QEngineT instantiation, pure */
    declare_ideal_engine<qpp::ket>(m);
    /* qpp::QEngineT instantiation, mixed */
    declare_ideal_engine<qpp::cmat>(m);

    /* Noise models instantiators */
    declare_noise_model<qpp::QubitBitFlipNoise, realT>(m, "QubitBitFlipNoise");
    declare_noise_model<qpp::QubitPhaseFlipNoise, realT>(m,
                                                         "QubitPhaseFlipNoise");
    declare_noise_model<qpp::QubitBitPhaseFlipNoise, realT>(
        m, "QubitBitPhaseFlipNoise");
    declare_noise_model<qpp::QubitDepolarizingNoise, realT>(
        m, "QubitDepolarizingNoise");
    declare_noise_model<qpp::QubitAmplitudeDampingNoise, realT>(
        m, "QubitAmplitudeDampingNoise");
    declare_noise_model<qpp::QubitPhaseDampingNoise, realT>(
        m, "QubitPhaseDampingNoise");
    declare_noise_model<qpp::QuditDepolarizingNoise, realT, idx>(
        m, "QuditDepolarizingNoise");

    /* qpp::QNoisyEngineT instantiations with different noise models, pure */
    declare_noisy_engine<qpp::ket, qpp::QubitBitFlipNoise>(m,
                                                           "QubitBitFlipNoise");
    declare_noisy_engine<qpp::ket, qpp::QubitBitPhaseFlipNoise>(
        m, "QubitBitPhaseFlipNoise");
    declare_noisy_engine<qpp::ket, qpp::QubitDepolarizingNoise>(
        m, "QubitDepolarizingNoise");
    declare_noisy_engine<qpp::ket, qpp::QubitPhaseFlipNoise>(
        m, "QubitPhaseFlipNoise");
    declare_noisy_engine<qpp::ket, qpp::QubitAmplitudeDampingNoise>(
        m, "QubitAmplitudeDampingNoise");
    declare_noisy_engine<qpp::ket, qpp::QubitPhaseDampingNoise>(
        m, "QubitPhaseDampingNoise");
    declare_noisy_engine<qpp::ket, qpp::QuditDepolarizingNoise>(
        m, "QuditDepolarizingNoise");

    /* qpp::QNoisyEngineT instantiations with different noise models, mixed */
    declare_noisy_engine<qpp::cmat, qpp::QubitBitFlipNoise>(
        m, "QubitBitFlipNoise");
    declare_noisy_engine<qpp::cmat, qpp::QubitBitPhaseFlipNoise>(
        m, "QubitBitPhaseFlipNoise");
    declare_noisy_engine<qpp::cmat, qpp::QubitDepolarizingNoise>(
        m, "QubitDepolarizingNoise");
    declare_noisy_engine<qpp::cmat, qpp::QubitPhaseFlipNoise>(
        m, "QubitPhaseFlipNoise");
    declare_noisy_engine<qpp::cmat, qpp::QubitAmplitudeDampingNoise>(
        m, "QubitAmplitudeDampingNoise");
    declare_noisy_engine<qpp::cmat, qpp::QubitPhaseDampingNoise>(
        m, "QubitPhaseDampingNoise");
    declare_noisy_engine<qpp::cmat, qpp::QuditDepolarizingNoise>(
        m, "QuditDepolarizingNoise");
}

#endif /* PYQPP_CLASSES_ENGINES_BIND_HPP_ */
