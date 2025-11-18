/*
 * This file is part of Quantum++.
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

/**
 * \file qpp/classes/qnoisy_engine.hpp
 * \brief Noisy quantum engines
 */

#ifndef QPP_CLASSES_ENGINES_HPP_
#define QPP_CLASSES_ENGINES_HPP_

#include <optional>
#include <string>
#include <vector>

#include "qpp/types.hpp"

#include "qpp/classes/exception.hpp"
#include "qpp/classes/qcircuit.hpp"
#include "qpp/classes/qengine.hpp"

namespace qpp {
/**
 * \class qpp::QNoisyEngineT
 * \brief Noisy quantum circuit engine, executes qpp::QCircuit
 * \see qpp::QEngineT, qpp::NoiseBase, qpp::QCircuit
 *
 * Assumes an uncorrelated noise model that is applied to each non-measured
 * qubit before every non-measurement step in the logical circuit. To add
 * noise before a measurement, insert a no-op via qpp::QCircuit::nop().
 *
 * \tparam T Engine's state underlying type, qpp::ket or qpp::cmat
 * \tparam NoiseModel Quantum noise model, should be derived from
 * qpp::NoiseBase
 */
template <typename T, typename NoiseModel>
class QNoisyEngineT : public QEngineT<T> {
    NoiseModel noise_;                            ///< quantum noise model
    std::vector<std::vector<idx>> noise_results_; ///< noise results
  public:
    using QBaseEngine<T, QCircuit>::execute;
    /**
     * \brief Constructs a noisy quantum engine out of a quantum circuit
     * description
     *
     * \param qc Quantum circuit description
     * \param noise Quantum noise model
     */
    explicit QNoisyEngineT(const QCircuit& qc, const NoiseModel& noise)
        : QEngineT<T>{qc}, noise_{noise}, noise_results_(qc.get_step_count()) {
        // EXCEPTION CHECKS
        // check noise has the correct dimensionality
        if (qc.get_d() != noise.get_d()) {
            throw exception::DimsNotEqual("qpp::QNoisyEngineT::QNoisyEngineT()",
                                          "noise");
        }
        // END EXCEPTION CHECKS
    }

    // traits
    /**
     * \brief qpp::IQEngineTraits::traits_get_name() override
     */
    std::string traits_get_name() const override { return "QNoisyEngineT"; }

    /**
     * \brief qpp::IQEngineTraits::traits_is_noisy() override
     */
    bool traits_is_noisy() const override { return true; }
    // end traits

    /**
     * \brief Executes one step in the quantum circuit description
     *
     * \note Override only this member function in every derived class to
     * achieve the desired behaviour
     *
     * \param it Iterator pointing to the step to be executed
     * \return Reference to the current instance
     */
    QNoisyEngineT&
    execute(typename QCircuitTraits<QCircuit>::iterator_type& it) override {
        // get the relative position of the target
        std::vector<idx> target_rel_pos =
            this->get_relative_pos_(this->get_non_measured_d());

        // apply the noise
        for (idx i : target_rel_pos) {
            this->qeng_st_.qstate_ = noise_(this->qeng_st_.qstate_, i);
            // record the Kraus operator that occurred
            noise_results_[it->get_ip()].emplace_back(noise_.get_last_idx());
        }

        // execute the circuit step
        (void)QEngineT<T>::execute(it);

        return *this;
    }

    /**
     * \brief Executes the entire quantum circuit description
     *
     * \param reps Number of repetitions
     * \return Reference to the current instance
     */
    QNoisyEngineT& execute(idx reps = 1) override {
        // EXCEPTION CHECKS
        if (reps == 0) {
            throw exception::OutOfRange("qpp::QNoisyEngineT::execute()",
                                        "reps");
        }
        // END EXCEPTION CHECKS
        auto steps = internal::circuit_as_iterators(*this->qc_ptr_);
        this->execute_prj_steps_no_sample_(steps, 0, reps);

        return *this;
    }

    /**
     * \brief Resets the engine
     *
     * Re-initializes everything to zero and sets the initial state to
     * \f$|0\rangle^{\otimes n}\f$
     *
     * \param reset_stats Optional (true by default), resets the collected
     * measurement statistics hash table
     * \return Reference to the current instance
     */
    QNoisyEngineT& reset(std::optional<T> qstate = std::nullopt,
                         bool reset_stats = true) override {
        QEngineT<T>::reset(qstate, reset_stats);
        noise_results_ = {};

        return *this;
    }

    // getters
    /**
     * \brief Vector of noise results obtained before every step in the
     * circuit
     *
     * \note The first vector contains the noise measurement results
     * obtained before applying the first step in the circuit, and so on,
     * ordered by non-measured qudits. That is, the first element in the
     * vector corresponding to noise obtained before a given step in the
     * circuit represents the noise result obtained on the first
     * non-measured qudit etc.
     *
     * \return Vector of noise results
     */
    std::vector<std::vector<idx>> get_noise_results() const {
        return noise_results_;
    }
    // end getters
}; /* class QNoisyEngineT */

/**
 * \class qpp::QNoisyEngine
 * \brief Pure state noisy quantum engine
 * \note Kept for backwards compatibility, use qpp::QKetNoisyEngine
 * \see qpp::QNoisyEngineT
 */
template <typename NoiseModel>
struct QNoisyEngine : public QNoisyEngineT<ket, NoiseModel> {
    using QNoisyEngineT<ket, NoiseModel>::QNoisyEngineT;
    // traits
    /**
     * \brief qpp::IQEngineTraits::traits_get_name() override
     */
    std::string traits_get_name() const override { return "QNoisyEngine"; }
    // end traits
};
// template deduction rule
template <class NoiseModel>
QNoisyEngine(const qpp::QCircuit& qc, const NoiseModel& noise)
    -> QNoisyEngine<NoiseModel>;

/**
 * \class qpp::QKetNoisyEngine
 * \brief Pure state noisy quantum engine
 * \see qpp::QNoisyEngineT
 */
template <typename NoiseModel>
struct QKetNoisyEngine : public QNoisyEngineT<ket, NoiseModel> {
    using QNoisyEngineT<ket, NoiseModel>::QNoisyEngineT;
    // traits
    /**
     * \brief qpp::IQEngineTraits::traits_get_name() override
     */
    std::string traits_get_name() const override { return "QKetNoisyEngine"; }
    // end traits
};
// template deduction rule
template <class NoiseModel>
QKetNoisyEngine(const qpp::QCircuit& qc, const NoiseModel& noise)
    -> QKetNoisyEngine<NoiseModel>;

/**
 * \class qpp::QDensityNoisyEngine
 * \brief Mixed state noisy quantum engine
 * \see qpp::QNoisyEngineT
 */
template <typename NoiseModel>
struct QDensityNoisyEngine : public QNoisyEngineT<cmat, NoiseModel> {
    using QNoisyEngineT<cmat, NoiseModel>::QNoisyEngineT;
    // traits
    /**
     * \brief qpp::IQEngineTraits::traits_get_name() override
     */
    std::string traits_get_name() const override {
        return "QDensityNoisyEngine";
    }
    // end traits
};
// template deduction rule
template <class NoiseModel>
QDensityNoisyEngine(const qpp::QCircuit& qc, const NoiseModel& noise)
    -> QDensityNoisyEngine<NoiseModel>;

} /* namespace qpp */

#endif /* QPP_CLASSES_ENGINES_HPP_ */
