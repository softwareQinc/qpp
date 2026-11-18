
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
 * @file <pyqpp/classes/qnoisy_engine_bind.hpp>
 * @brief Bindings for <qpp/classes/qnoisy_engine.hpp>
 */

#ifndef PYQPP_CLASSES_QNOISY_ENGINE_BIND_HPP_
#define PYQPP_CLASSES_QNOISY_ENGINE_BIND_HPP_

#include "pyqpp/pyqpp_common.hpp"

/* Noise models instantiators */
template <typename NoiseModel, typename... CtorTypeList>
void declare_noise_model(py::module& m, const std::string& type) {
    using namespace qpp;
    auto pyNoiseModel = py::class_<NoiseModel>(m, type.c_str());

    pyNoiseModel.def(py::init<CtorTypeList...>());

    pyNoiseModel.def("get_d", &NoiseModel::get_d, "Qudit dimension");
    pyNoiseModel.def("get_Ks", &NoiseModel::get_Ks,
                     "Vector of noise operators");
    pyNoiseModel.def("get_last_idx", &NoiseModel::get_last_idx,
                     "Index of the last occurring noise element");
    pyNoiseModel.def("get_last_K", &NoiseModel::get_last_K,
                     "Last occurring noise element");
    pyNoiseModel.def("get_last_p", &NoiseModel::get_last_p,
                     "Probability of the last occurring noise element");
    pyNoiseModel.def(
        "get_probs", &NoiseModel::get_probs,
        "Vector of probabilities corresponding to each noise operator");
}

/* qpp::QNoisyEngineT instantiator */
template <typename T, typename NoiseModel>
void declare_QNoisyEngineT(py::module& m, const std::string& type) {
    using namespace qpp;
    static_assert(std::is_same_v<T, ket> || std::is_same_v<T, cmat>,
                  "The underlying type must be qpp::ket or qpp::cmat");

    std::string pyname;
    if constexpr (std::is_same_v<T, ket>) {
        pyname = "_QKetNoisyEngine_" + type;
    } else {
        pyname = "_QDensityNoisyEngine_" + type;
    }

    auto pyQNoisyEngineT =
        py::class_<QNoisyEngineT<T, NoiseModel>, QEngineT<T>>(m,
                                                              pyname.c_str());

    // keep the circuit alive for the lifetime of the engine
    pyQNoisyEngineT.def(py::init<const QCircuit&, const NoiseModel&>(),
                        py::keep_alive<1, 2>());

    pyQNoisyEngineT.def(
        "get_noise_results", &QNoisyEngineT<T, NoiseModel>::get_noise_results,
        "Vector of noise results obtained before every step in the circuit");

    // copies retain the original engine, which keeps the circuit alive
    pyQNoisyEngineT.def(
        "__copy__",
        [](const QNoisyEngineT<T, NoiseModel>& self) {
            return std::make_unique<QNoisyEngineT<T, NoiseModel>>(self);
        },
        py::keep_alive<0, 1>());

    pyQNoisyEngineT.def(
        "__deepcopy__",
        [](const QNoisyEngineT<T, NoiseModel>& self, py::dict) {
            return std::make_unique<QNoisyEngineT<T, NoiseModel>>(self);
        },
        py::keep_alive<0, 1>());
}

/* Noisy engine factory */
inline void declare_noisy_engine_factory(py::module& m, const std::string& name,
                                         const std::string& prefix) {
    // map each noise model type to its registered engine class
    py::dict classes;

    for (const char* type :
         {"QubitBitFlipNoise", "QubitPhaseFlipNoise", "QubitBitPhaseFlipNoise",
          "QubitDepolarizingNoise", "QubitAmplitudeDampingNoise",
          "QubitPhaseDampingNoise", "QuditDepolarizingNoise"}) {
        classes[m.attr(type)] = m.attr((prefix + type).c_str());
    }

    // construct through the bound class to preserve keep_alive
    m.def(
        name.c_str(),
        [classes](py::object qc, py::object noise) {
            py::object type = py::type::of(noise);

            if (!classes.contains(type)) {
                throw py::type_error("Unsupported noise model");
            }

            py::object cls = py::reinterpret_borrow<py::object>(classes[type]);

            return cls(qc, noise);
        },
        py::arg("qc"), py::arg("noise"));
}

/* qpp::QNoisyEngine */
inline void init_classes_qnoisy_engine(py::module_& m) {
    using namespace qpp;

    /* Noise models instantiators */
    declare_noise_model<QubitBitFlipNoise, realT>(m, "QubitBitFlipNoise");
    declare_noise_model<QubitPhaseFlipNoise, realT>(m, "QubitPhaseFlipNoise");
    declare_noise_model<QubitBitPhaseFlipNoise, realT>(
        m, "QubitBitPhaseFlipNoise");
    declare_noise_model<QubitDepolarizingNoise, realT>(
        m, "QubitDepolarizingNoise");
    declare_noise_model<QubitAmplitudeDampingNoise, realT>(
        m, "QubitAmplitudeDampingNoise");
    declare_noise_model<QubitPhaseDampingNoise, realT>(
        m, "QubitPhaseDampingNoise");
    declare_noise_model<QuditDepolarizingNoise, realT, idx>(
        m, "QuditDepolarizingNoise");

    /* qpp::QNoisyEngineT instantiations with different noise models, pure */
    declare_QNoisyEngineT<ket, QubitBitFlipNoise>(m, "QubitBitFlipNoise");
    declare_QNoisyEngineT<ket, QubitBitPhaseFlipNoise>(
        m, "QubitBitPhaseFlipNoise");
    declare_QNoisyEngineT<ket, QubitDepolarizingNoise>(
        m, "QubitDepolarizingNoise");
    declare_QNoisyEngineT<ket, QubitPhaseFlipNoise>(m, "QubitPhaseFlipNoise");
    declare_QNoisyEngineT<ket, QubitAmplitudeDampingNoise>(
        m, "QubitAmplitudeDampingNoise");
    declare_QNoisyEngineT<ket, QubitPhaseDampingNoise>(
        m, "QubitPhaseDampingNoise");
    declare_QNoisyEngineT<ket, QuditDepolarizingNoise>(
        m, "QuditDepolarizingNoise");

    /* qpp::QNoisyEngineT instantiations with different noise models, mixed */
    declare_QNoisyEngineT<cmat, QubitBitFlipNoise>(m, "QubitBitFlipNoise");
    declare_QNoisyEngineT<cmat, QubitBitPhaseFlipNoise>(
        m, "QubitBitPhaseFlipNoise");
    declare_QNoisyEngineT<cmat, QubitDepolarizingNoise>(
        m, "QubitDepolarizingNoise");
    declare_QNoisyEngineT<cmat, QubitPhaseFlipNoise>(m, "QubitPhaseFlipNoise");
    declare_QNoisyEngineT<cmat, QubitAmplitudeDampingNoise>(
        m, "QubitAmplitudeDampingNoise");
    declare_QNoisyEngineT<cmat, QubitPhaseDampingNoise>(
        m, "QubitPhaseDampingNoise");
    declare_QNoisyEngineT<cmat, QuditDepolarizingNoise>(
        m, "QuditDepolarizingNoise");

    /* Public factory functions */
    declare_noisy_engine_factory(m, "QKetNoisyEngine", "_QKetNoisyEngine_");
    declare_noisy_engine_factory(m, "QDensityNoisyEngine",
                                 "_QDensityNoisyEngine_");

    // backwards compatibility
    m.attr("QNoisyEngine") = m.attr("QKetNoisyEngine");
}

#endif /* PYQPP_CLASSES_QNOISY_ENGINE_BIND_HPP_ */
