/*
 * This file is part of Quantum++.
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

/**
 * \file qpp/classes/qdummy_engine.hpp
 * \brief No-op (dummy) quantum engines
 */

#ifndef QPP_CLASSES_QDUMMY_ENGINE_HPP_
#define QPP_CLASSES_QDUMMY_ENGINE_HPP_

#include <string>
#include <type_traits>

#include "qpp/types.hpp"

#include "qpp/classes/qbase_engine.hpp"
#include "qpp/classes/qcircuit.hpp"

namespace qpp {

/**
 * \class qpp::QDummyEngine
 * \brief No-op (dummy) quantum engine
 * \see qpp::QBaseEngine
 *
 * \tparam T Engine's state underlying type
 * \tparam QCT Circuit underlying type
 */
template <typename T, typename QCT>
struct QDummyEngine : public QBaseEngine<T, QCT> {
    using QBaseEngine<T, QCT>::QBaseEngine;
    // traits
    /**
     * \brief qpp::IQEngineTraits::traits_get_name() override
     */
    std::string traits_get_name() const override { return "QDummyEngine"; }

    /**
     * \brief qpp::IQEngineTraits::traits_is_noisy() override
     */
    bool traits_is_noisy() const override { return false; }

    /**
     * \brief qpp::IQEngineTraits::traits_is_pure() override
     */
    bool traits_is_pure() const override {
        if (std::is_same_v<T, cmat>) {
            return false;
        } else if (std::is_same_v<T, ket>) {
            return true;
        }
        // default, we assume pure states
        return true;
    }
    // end traits
}; /* struct QDummyEngine */

/**
 * \class qpp::QKetDummyEngine
 * \brief Pure state no-op (dummy) quantum engine for qpp::QCircuit
 * \see qpp::QDummyEngine
 */
struct QKetDummyEngine : public QDummyEngine<ket, QCircuit> {
    using QDummyEngine<ket, QCircuit>::QDummyEngine;
    // traits
    std::string traits_get_name() const override { return "QKetDummyEngine"; }
    // end traits
};

/**
 * \class qpp::QDensityDummyEngine
 * \brief Mixed state no-op (dummy) quantum engine for qpp::QCircuit
 * \see qpp::QDummyEngine
 */
struct QDensityDummyEngine : public QDummyEngine<cmat, QCircuit> {
    using QDummyEngine<cmat, QCircuit>::QDummyEngine;
    // traits
    std::string traits_get_name() const override {
        return "QDensityDummyEngine";
    }
    // end traits
};
} /* namespace qpp */

#endif /* QPP_CLASSES_QDUMMY_ENGINE_HPP_ */
