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
 * \file qpp/internal/classes/qengine_state.hpp
 * \brief Quantum engine state
 */

#ifndef QPP_INTERNAL_CLASSES_QENGINE_STATE_HPP_
#define QPP_INTERNAL_CLASSES_QENGINE_STATE_HPP_

#include <numeric>
#include <optional>
#include <type_traits>
#include <vector>

#include "qpp/classes/qcircuit.hpp"
#include "qpp/classes/states.hpp"

namespace qpp {
namespace internal {
/**
 * \class qpp::internal::QEngineState
 * \brief Current state of qpp::QEngineT and derived
 *
 * \tparam T State underlying type, qpp::ket or qpp::cmat
 */
template <typename T>
struct QEngineState {
    static_assert(std::is_same_v<T, ket> || std::is_same_v<T, cmat>,
                  "The underlying type must be qpp::ket or qpp::cmat");

    const QCircuit* qc_ptr_;     ///< non-owning pointer to the parent
                                 ///< const quantum circuit description
    T qstate_{};                 ///< state vector or density matrix
    std::vector<realT> probs_{}; ///< measurement probabilities
    std::vector<idx> dits_{};    ///< classical dits (where measurement
                                 ///< results are usually stored)
    std::vector<idx> subsys_{};  ///< keeps track of the destructively measured
                                 ///< subsystems, re-label them after
                                 ///< measurements
    bool post_select_ok_ =
        true; ///< flag that becomes false if/when post-selection fails

    /**
     * \brief Constructor
     *
     * \param qc_ptr Non-owning pointer to the parent const quantum circuit
     * description
     */
    explicit QEngineState(const QCircuit* qc_ptr) : qc_ptr_{qc_ptr} {
        // EXCEPTION CHECKS
        if (qc_ptr->get_nq() == 0) {
            throw exception::ZeroSize(
                "qpp::internal::QEngineState::QEngineState()", "nq");
        }
        // END EXCEPTION CHECKS
        reset();
    }

    // silence -Weffc++ class has pointer data members
    /**
     * \brief Default copy constructor
     */
    QEngineState(const QEngineState&) = default;

    // silence -Weffc++ class has pointer data members
    /**
     * \brief Default copy assignment operator
     *
     * \return Reference to the current instance
     */
    QEngineState& operator=(const QEngineState&) = default;

    /**
     * \brief Resets the engine state
     *
     * \param qstate Optional engine's initial quantum state
     */
    void reset(std::optional<T> qstate = std::nullopt) {
        if (qstate.has_value()) {
            idx D = internal::safe_pow(qc_ptr_->get_d(), qc_ptr_->get_nq());
            // EXCEPTION CHECKS
            if constexpr (std::is_same_v<T, ket>) {
                if (static_cast<idx>(qstate.value().rows()) != D) {
                    throw exception::DimsNotEqual(
                        "qpp::internal::QEngineState::reset()", "state");
                }
            } else {
                if (!internal::check_square_mat(qstate.value())) {
                    throw exception::MatrixNotSquare("qpp::QEngineState::reset",
                                                     "state");
                }
                if (static_cast<idx>(qstate.value().rows()) != D) {
                    throw exception::DimsNotEqual("qpp::QEngineState::reset()",
                                                  "state");
                }
            }
            // END EXCEPTION CHECKS
            qstate_ = qstate.value();
        } else {
            if constexpr (std::is_same_v<T, ket>) {
                qstate_ = States::get_no_thread_local_instance().zero(
                    qc_ptr_->get_nq(), qc_ptr_->get_d());
            } else {
                qstate_ = prj(States::get_no_thread_local_instance().zero(
                    qc_ptr_->get_nq(), qc_ptr_->get_d()));
            }
        }
        probs_ = std::vector<realT>(qc_ptr_->get_nc(), 0);
        dits_ = std::vector<idx>(qc_ptr_->get_nc(), 0);
        subsys_ = std::vector<idx>(qc_ptr_->get_nq(), 0);
        std::iota(subsys_.begin(), subsys_.end(), 0);
        post_select_ok_ = true;
    }
}; /* struct QEngineState */

} /* namespace internal */
} /* namespace qpp */

#endif /* QPP_INTERNAL_CLASSES_QENGINE_STATE_HPP_ */
