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
 * \file qpp/classes/qcircuit.hpp
 * \brief Quantum circuits
 */

#ifndef QPP_CLASSES_QCIRCUIT_HPP_
#define QPP_CLASSES_QCIRCUIT_HPP_

#include <algorithm>
#include <cstddef>
#include <initializer_list>
#include <iomanip>
#include <iterator>
#include <optional>
#include <ostream>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <variant>
#include <vector>

#include "qpp/functions.hpp"
#include "qpp/input_output.hpp"
#include "qpp/types.hpp"

#include "qpp/classes/exception.hpp"
#include "qpp/classes/gates.hpp"
#include "qpp/classes/idisplay.hpp"
#include "qpp/classes/ijson.hpp"

#include "qpp/internal/classes/qcircuit_conditional_step.hpp"
#include "qpp/internal/classes/qcircuit_gate_step.hpp"
#include "qpp/internal/classes/qcircuit_measurement_step.hpp"
#include "qpp/internal/classes/qcircuit_nop_step.hpp"
#include "qpp/internal/classes/qcircuit_resources.hpp"
#include "qpp/internal/util.hpp"

namespace qpp {
// forward declaration
namespace internal {
class QCircuitIterator;
}
/**
 * \class qpp::QCircuit
 * \brief Quantum circuit description
 * \see qpp::QEngineT
 */
class QCircuit : public IDisplay, public IJSON {
    friend class internal::QCircuitIterator;
    idx nq_;                          ///< number of qudits
    idx nc_;                          ///< number of classical "dits"
    idx d_;                           ///< qudit dimension
    std::optional<std::string> name_; ///< optional circuit name

    std::vector<bool>
        measured_d_; ///< keeps track of the destructively measured qudits
    std::vector<bool>
        measured_nd_; ///< keeps track of the non-destructively measured qudits
    std::vector<bool> clean_qudits_; ///< keeps track of clean (unused) qudits

    std::vector<bool>
        clean_dits_; ///< keeps track of clean (unused) classical dits
    std::vector<bool> measurement_dits_; ///< keeps track of classical dits used
    ///< in measurements

    std::unordered_map<std::size_t, cmat>
        cmat_hash_tbl_{}; ///< hash table with the matrices used in the circuit,
                          ///< with [Key = std::size_t, Value = cmat]
    std::unordered_map<std::size_t, idx> gate_count_{}; ///< gate counts
    std::unordered_map<std::size_t, idx>
        measurement_count_{}; ///< measurement counts

    /**
     * \brief Adds matrix to the hash table
     *
     * \note Throws if a hash collision is detected, i.e., if two different
     * matrices have the same hash value
     *
     * \param hashU Hash value of U
     * \param U Complex matrix
     */
    void add_hash_(std::size_t hashU, const cmat& U) {
        // EXCEPTION CHECKS
        auto search = cmat_hash_tbl_.find(hashU);
        static internal::EqualEigen equal_eigen;
        if (search != cmat_hash_tbl_.end()) { // found the hash in the table
            // have a hash collision
            if (!equal_eigen(search->second, U)) {
                throw exception::CustomException("qpp::QCircuit::add_hash_()",
                                                 "Matrix hash collision");
            }
        }
        // END EXCEPTION CHECKS
        cmat_hash_tbl_.insert({hashU, U});
    }

    using VarStep =
        std::variant<internal::QCircuitConditionalStep,
                     internal::QCircuitGateStep,
                     internal::QCircuitMeasurementStep,
                     internal::QCircuitNOPStep>; ///< circuit step variant
    std::vector<VarStep> circuit_; ///<< quantum circuit representation

  protected:
    /**
     * \brief Visits a qpp::QCircuit::VarStep variant
     *
     * \tparam Fns Callable(s) matching the variant
     *
     * \param circuit_step qpp::QCircuit::VarStep variant
     * \param fns Callable(s) matching the variant
     */
    template <class... Fns>
    void visit_circuit_step_(VarStep& circuit_step, Fns&&... fns);

    /**
     * \brief Visits a qpp::QCircuit::VarStep variant
     *
     * \tparam Fns Callable(s) matching the variant
     *
     * \param circuit_step qpp::QCircuit::VarStep variant
     * \param fns Callable(s) matching the variant
     */
    template <class... Fns>
    void visit_circuit_step_(const VarStep& circuit_step, Fns&&... fns) const;

    /**
     * \brief Visits a vector of qpp::QCircuit::VarStep variants
     *
     * \tparam Fns Callable(s) matching the variant
     *
     * \param circuit Vector of qpp::QCircuit::VarStep variants
     * \param fns Callable(s) matching the variant
     */
    template <class... Fns>
    void visit_circuit_(std::vector<VarStep>& circuit, Fns&&... fns);

    /**
     * \brief Visits a vector of qpp::QCircuit::VarStep variants
     *
     * \tparam Fns Callable(s) matching the variant
     *
     * \param circuit Vector of qpp::QCircuit::VarStep variants
     * \param fns Callable(s) matching the variant
     */
    template <class... Fns>
    void visit_circuit_(const std::vector<VarStep>& circuit,
                        Fns&&... fns) const;

  public:
    ///< both iterators are const_iterators
    using iterator = internal::QCircuitIterator;
    using const_iterator = iterator;

    /**
     * \brief Iterator to the first element
     *
     * \return Iterator to the first element
     */
    iterator begin() noexcept;

    /**
     * \brief Iterator to the next to the last element
     *
     * \return Iterator to the next to the last element
     */
    iterator end() noexcept;

    /**
     * \brief Constant iterator to the first element
     *
     * \return Constant iterator to the first element
     */
    const_iterator begin() const noexcept;

    /**
     * \brief Constant iterator to the next to the last element
     *
     * \return Constant iterator to the next to the last element
     */
    const_iterator end() const noexcept;

    /**
     * \brief Constant iterator to the first element
     *
     * \return Constant iterator to the first element
     */
    const_iterator cbegin() const noexcept;

    /**
     * \brief Constant iterator to the next to the last element
     *
     * \return Constant iterator to the next to the last element
     */
    const_iterator cend() const noexcept;

    /**
     * \brief Constructs a quantum circuit description
     *
     * \note The measurement results can only be stored in the classical
     * dits of which number is specified by \a nc
     *
     * \param nq Number of qudits (defaults to 1 so qpp::QCircuit is
     * default-constructible)
     * \param nc Number of classical dits
     * \param d Subsystem dimensions (default is qubit, i.e., \a d =2)
     * \param name Optional circuit name
     */
    explicit QCircuit(idx nq = 1, idx nc = 0, idx d = 2,
                      std::optional<std::string> name = std::nullopt)
        : nq_{nq}, nc_{nc}, d_{d}, name_{std::move(name)},
          measured_d_(nq, false), measured_nd_(nq, false),
          clean_qudits_(nq_, true), clean_dits_(nc_, true),
          measurement_dits_(nc_, false), circuit_{} {
        // EXCEPTION CHECKS
        // if (nq == 0)
        //    throw exception::ZeroSize("qpp::QCircuit::QCircuit()", "nq");
        if (d < 2) {
            throw exception::OutOfRange("qpp::QCircuit::QCircuit()", "d");
        }
        // END EXCEPTION CHECKS
    }

    /**
     * \brief Default virtual destructor
     */
    ~QCircuit() override = default;

    // static member functions
    /**
     * \brief Checks whether a gate step is a controlled gate
     *
     * \return True if the gate step is a controlled gate, false otherwise
     */
    inline static bool is_CTRL(const internal::QCircuitGateStep& gate_step) {
        switch (gate_step.gate_type_) {
            case internal::QCircuitGateStep::Type::CTRL:
            case internal::QCircuitGateStep::Type::CTRL_FAN:
                return true;
            default:
                return false;
        }
    }

    /**
     * \brief Checks whether a gate step is a classically-controlled gate
     *
     * \return True if the gate step is a classically-controlled gate, false
     * otherwise
     */
    inline static bool is_cCTRL(const internal::QCircuitGateStep& gate_step) {
        switch (gate_step.gate_type_) {
            case internal::QCircuitGateStep::Type::cCTRL:
            case internal::QCircuitGateStep::Type::cCTRL_FAN:
                return true;
            default:
                return false;
        }
    }

    // getters
    /**
     * \brief Hash table with the matrices used in the circuit
     *
     * \return Hash table with the matrices used in the circuit
     */
    const std::unordered_map<std::size_t, cmat>&
    get_cmat_hash_tbl() const noexcept {
        return cmat_hash_tbl_;
    }

    /**
     * \brief Total number of qudits in the circuit
     *
     * \return Total number of qudits
     */
    idx get_nq() const noexcept { return nq_; }

    /**
     * \brief Total number of classical dits in the circuit
     *
     * \return Total number of classical dits
     */
    idx get_nc() const noexcept { return nc_; }

    /**
     * \brief Dimension of the comprising qudits
     *
     * \return Qudit dimension
     */
    idx get_d() const noexcept { return d_; }

    /**
     * \brief Quantum circuit description name
     *
     * \return Quantum circuit description name
     */
    std::optional<std::string> get_name() const { return name_; }

    /**
     * \brief Check whether qudit \a i was already measured (destructively)
     *
     * \param i Qudit index
     * \return True if qudit \a i was already measured (destructively),
     * false otherwise
     */
    bool was_measured_d(idx i) const {
        // EXCEPTION CHECKS
        if (i >= nq_) {
            throw exception::OutOfRange("qpp::QCircuit::was_measured_d()", "i");
        }
        // END EXCEPTION CHECKS

        return measured_d_[i];
    }

    /**
     * \brief Vector of already measured (destructively) qudit indexes
     *
     * \return Vector of already measured (destructively) qudit indexes
     */
    std::vector<idx> get_measured_d() const {
        std::vector<idx> result;
        for (idx i = 0; i < nq_; ++i) {
            if (was_measured_d(i)) {
                result.emplace_back(i);
            }
        }

        return result;
    }

    /**
     * \brief Check whether qudit \a i was already measured (non-destructively)
     *
     * \param i Qudit index
     * \return True if qudit \a i was already measured (non-destructively),
     * false otherwise
     */
    bool was_measured_nd(idx i) const {
        // EXCEPTION CHECKS
        if (i >= nq_) {
            throw exception::OutOfRange("qpp::QCircuit::was_measured_nd()",
                                        "i");
        }
        // END EXCEPTION CHECKS

        return measured_nd_[i];
    }

    /**
     * \brief Vector of already measured (non-destructively) qudit indexes
     *
     * \return Vector of already measured (non-destructively) qudit indexes
     */
    std::vector<idx> get_measured_nd() const {
        std::vector<idx> result;
        for (idx i = 0; i < nq_; ++i) {
            if (was_measured_nd(i)) {
                result.emplace_back(i);
            }
        }

        return result;
    }

    /**
     * \brief Vector of non-measured (destructively) qudit indexes
     *
     * \return Vector of non-measured (destructively) qudit indexes
     */
    std::vector<idx> get_non_measured_d() const {
        std::vector<idx> result;
        for (idx i = 0; i < nq_; ++i) {
            if (!was_measured_d(i)) {
                result.emplace_back(i);
            }
        }

        return result;
    }

    /**
     * \brief Vector of classical dits that were used to store results of
     * measurements (either destructive or non-destructive)
     *
     * \return Vector of classical dits that participate in measurements
     */
    std::vector<idx> get_measurement_dits() const {
        std::vector<idx> result;
        for (idx i = 0; i < nc_; ++i) {
            if (is_measurement_dit(i)) {
                result.emplace_back(i);
            }
        }

        return result;
    }

    /**
     * \brief Quantum circuit description gate count
     *
     * \param U Optional, gate. If absent (default), this function computes the
     * total circuit gate count.
     * \return Gate count
     */
    idx get_gate_count(std::optional<cmat> U = std::nullopt) const {
        // EXCEPTION CHECKS
        // square matrix
        if (U.has_value()) {
            if (!internal::check_square_mat(U.value())) {
                throw exception::MatrixNotSquare("qpp::get_gate_count()", "U");
            }
        }
        // END EXCEPTION CHECKS

        if (U.has_value()) {
            idx result = 0;
            std::size_t hashU = hash_eigen(U.value());
            try {
                result = gate_count_.at(hashU);
            } // gate hash not found in the hash table
            catch (...) {
                return 0;
            }
            return result;
        }

        idx result = 0;
        for (auto&& elem : gate_count_) {
            result += elem.second;
        }

        return result;
    }

    /**
     * \brief Quantum circuit description gate depth
     *
     * \param U Optional, gate. If absent (default), this function computes the
     * total circuit gate depth.
     * \return Gate depth
     */
    idx get_gate_depth(std::optional<cmat> U = std::nullopt) const {
        // EXCEPTION CHECKS
        // square matrix
        if (U.has_value() && !internal::check_square_mat(U.value())) {
            throw exception::MatrixNotSquare("qpp::get_gate_depth()", "U");
        }
        // END EXCEPTION CHECKS

        bool found = false;
        std::vector<idx> heights(nc_ + nq_, 0);

        std::size_t hashU = U.has_value() ? hash_eigen(U.value()) : 0;
        // iterate over all steps in the circuit

        for (auto&& step : circuit_) {
            if (std::holds_alternative<internal::QCircuitGateStep>(step)) {
                auto gate_step = std::get<internal::QCircuitGateStep>(step);
                if (U.has_value() && gate_step.gate_hash_ != hashU) {
                    continue; // we skip this gate_step elem
                }

                found = true; // gate_step was found in the circuit

                std::vector<idx> ctrl =
                    gate_step.ctrl_.value_or(std::vector<idx>{});
                std::vector<idx> target = gate_step.target_;
                std::vector<idx> ctrl_target;
                ctrl_target.reserve(ctrl.size() + target.size());
                ctrl_target.insert(ctrl_target.end(), ctrl.begin(), ctrl.end());
                ctrl_target.insert(ctrl_target.end(), target.begin(),
                                   target.end());

                idx max_height = 0;

                if (is_cCTRL(gate_step)) {
                    // compute the "height" of the to-be-placed
                    // gate_step
                    for (idx i : ctrl) {
                        if (heights[i] > max_height) {
                            max_height = heights[i];
                        }
                    }
                    for (idx i : target) {
                        if (heights[nc_ + i] > max_height) {
                            max_height = heights[nc_ + i];
                        }
                    }
                    // apply classical ctrl
                    for (idx i : ctrl) {
                        heights[i] = max_height + 1;
                    }
                    // apply gate_step
                    for (idx i : target) {
                        heights[nc_ + i] = max_height + 1;
                    }
                } else {
                    // compute the "height" of the to-be-placed
                    // gate_step
                    for (idx i : ctrl_target) {
                        if (heights[nc_ + i] > max_height) {
                            max_height = heights[nc_ + i];
                        }
                    }
                    // apply (ctrl) gate_step
                    for (idx i : ctrl_target) {
                        heights[nc_ + i] = max_height + 1;
                    }
                }
            } // end if
        } // end for

        return found ? *std::max_element(heights.begin(), heights.end()) : 0;
    }

    /**
     * \brief Quantum circuit description measurement depth
     *
     * \param V Optional, orthonormal basis or rank-1 projectors specified by
     * the columns of matrix \a V. If absent (default), this function computes
     * the total circuit measurement depth.
     * \return Measurement depth
     */
    idx get_measurement_depth(std::optional<cmat> V = std::nullopt) const {
        bool found = false;
        std::vector<idx> heights(nc_ + nq_, 0);

        std::size_t hashV = V.has_value() ? hash_eigen(V.value()) : 0;

        // TODO check this

        // iterate over all steps in the circuit
        for (auto&& step : circuit_) {
            if (std::holds_alternative<internal::QCircuitMeasurementStep>(
                    step)) {
                auto measurement_step =
                    std::get<internal::QCircuitMeasurementStep>(step);
                if (V.has_value() && measurement_step.mats_hash_[0] != hashV) {
                    continue; // we skip this measurement_step elem
                }

                found = true; // measurement_step was found in the circuit

                std::vector<idx> target = measurement_step.target_;
                idx c_reg = measurement_step.c_reg_;

                idx max_height = 0;
                switch (measurement_step.measurement_type_) {
                    case internal::QCircuitMeasurementStep::Type::NONE:

                    case internal::QCircuitMeasurementStep::Type::MEASURE_V:
                    case internal::QCircuitMeasurementStep::Type::
                        MEASURE_V_JOINT:
                    case internal::QCircuitMeasurementStep::Type::MEASURE_V_ND:
                    case internal::QCircuitMeasurementStep::Type::
                        MEASURE_V_JOINT_ND:
                    case internal::QCircuitMeasurementStep::Type::POST_SELECT_V:
                    case internal::QCircuitMeasurementStep::Type::
                        POST_SELECT_V_JOINT:
                    case internal::QCircuitMeasurementStep::Type::
                        POST_SELECT_V_ND:
                    case internal::QCircuitMeasurementStep::Type::
                        POST_SELECT_V_JOINT_ND:
                        // compute the "height" of the to-be-placed
                        // measurement_step
                        if (heights[c_reg] > max_height) {
                            max_height = heights[c_reg];
                        }
                        for (idx i : target) {
                            if (heights[nc_ + i] > max_height) {
                                max_height = heights[nc_ + i];
                            }
                        }
                        // apply measurement_step
                        heights[c_reg] = max_height + 1;
                        for (idx i : target) {
                            heights[nc_ + i] = max_height + 1;
                        }
                        break;

                    case internal::QCircuitMeasurementStep::Type::MEASURE:
                    case internal::QCircuitMeasurementStep::Type::MEASURE_MANY:
                    case internal::QCircuitMeasurementStep::Type::MEASURE_ND:
                    case internal::QCircuitMeasurementStep::Type::
                        MEASURE_MANY_ND:
                    case internal::QCircuitMeasurementStep::Type::POST_SELECT:
                    case internal::QCircuitMeasurementStep::Type::
                        POST_SELECT_MANY:
                    case internal::QCircuitMeasurementStep::Type::
                        POST_SELECT_ND:
                    case internal::QCircuitMeasurementStep::Type::
                        POST_SELECT_MANY_ND:
                        // compute the "height" of the to-be-placed
                        // measurement_step
                        for (idx i = 0; i < static_cast<idx>(target.size());
                             ++i) {
                            if (heights[c_reg + i] > max_height) {
                                max_height = heights[c_reg + i];
                            }
                        }
                        for (idx i : target) {
                            if (heights[nc_ + i] > max_height) {
                                max_height = heights[nc_ + i];
                            }
                        }
                        // apply measurement_step
                        for (idx i = 0; i < static_cast<idx>(target.size());
                             ++i) {
                            heights[c_reg + i] = max_height + 1;
                        }

                        for (idx i : target) {
                            heights[nc_ + i] = max_height + 1;
                        }
                        break;

                    case internal::QCircuitMeasurementStep::Type::RESET:
                    case internal::QCircuitMeasurementStep::Type::RESET_MANY:
                    case internal::QCircuitMeasurementStep::Type::DISCARD:
                    case internal::QCircuitMeasurementStep::Type::DISCARD_MANY:
                        for (idx i : target) {
                            if (heights[nc_ + i] > max_height) {
                                max_height = heights[nc_ + i];
                            }
                        }
                        // apply reset/discard
                        for (idx i : target) {
                            heights[nc_ + i] = max_height + 1;
                        }
                        break;
                } // end switch
            } // end if
        } // end for

        return found ? *std::max_element(heights.begin(), heights.end()) : 0;
    }

    // computes the depth greedily, measuring the "height" (depth) of the
    // "pieces" (gates) placed in a Tetris-like style
    /**
     * \brief Quantum circuit description total depth
     *
     * \return Gate/measurement total depth
     */
    idx get_depth() const { return get_gate_depth() + get_measurement_depth(); }

    /**
     * \brief Quantum circuit description measurement count
     *
     * \param V Optional, orthonormal basis or rank-1 projectors specified by
     * the columns of matrix \a V. If absent (default), this function computes
     * the total circuit measurement count.
     * \return Measurement count
     */
    idx get_measurement_count(std::optional<cmat> V = std::nullopt) const {
        // EXCEPTION CHECKS
        if (V.has_value()) {
            idx result = 0;
            std::size_t hashV = hash_eigen(V.value());
            // basis matrix hash not found in the hash table
            try {
                result = measurement_count_.at(hashV);
            } catch (...) {
                return 0;
            }
            return result;
        }

        idx result = 0;
        for (auto&& elem : measurement_count_) {
            result += elem.second;
        }

        return result;
    }

    /**
     * \brief Quantum circuit total steps count, i.e., the sum of gate count
     * and measurement count
     *
     * \return Total (gates + measurements) count
     */
    idx get_step_count() const noexcept { return circuit_.size(); }

    /**
     * \brief No-op count
     *
     * \return No-op count
     */
    idx get_nop_count() const {
        return std::count_if(circuit_.begin(), circuit_.end(), [](auto&& arg) {
            return std::holds_alternative<internal::QCircuitNOPStep>(arg);
        });
    }

    /**
     * \brief Quantum circuit resources
     *
     * \return Instance of \a internal::QCircuitResources
     */
    internal::QCircuitResources get_resources() const {
        internal::QCircuitResources result;
        result.nq = get_nq();
        result.nc = get_nc();
        result.d = get_d();
        result.name = get_name();
        result.step_count = get_step_count();
        result.gate_count = get_gate_count();
        result.gate_depth = get_gate_depth();
        result.measurement_count = get_measurement_count();
        result.measurement_depth = get_measurement_depth();
        result.total_depth = get_depth();

        return result;
    }

    // end getters

    // setters
    /**
     * \brief Sets name for the quantum circuit description
     *
     * \param name Quantum circuit description name
     * \return Reference to the current instance
     */
    QCircuit& set_name(const std::string& name) {
        name_ = name;
        return *this;
    }
    // end setters

    /**
     * \brief Adds \a n additional qudits before qudit \a pos
     *
     * \note Qudits with indexes greater or equal than the newly inserted
     * ones have their indexes automatically incremented
     *
     * \param n Number of qudits
     * \param pos Qudit index
     * \return Reference to the current instance
     */
    QCircuit& add_qudit(idx n, idx pos) {
        // EXCEPTION CHECKS
        if (pos > nq_) {
            throw exception::OutOfRange("qpp::QCircuit::add_qudit()", "pos");
        }
        // END EXCEPTION CHECKS

        nq_ += n;

        // update gate indexes
        auto gate_step_visitor =
            [&](internal::QCircuitGateStep& gate_step) { // update ctrl indexes
                if (is_CTRL(gate_step)) {
                    for (idx& elem : gate_step.ctrl_.value()) {
                        if (elem >= pos) {
                            elem += n;
                        }
                    }
                }
                // update target indexes
                for (idx& elem : gate_step.target_) {
                    if (elem >= pos) {
                        elem += n;
                    }
                }
            };
        // update measurement indexes
        auto measurement_step_visitor =
            [&](internal::QCircuitMeasurementStep& measurement_step) {
                for (idx& elem : measurement_step.target_) {
                    if (elem >= pos) {
                        elem += n;
                    }
                }
            };
        auto nop_step_visitor = [&](internal::QCircuitNOPStep&) {};
        auto conditional_step_visitor =
            [&](const internal::QCircuitConditionalStep&) {};

        visit_circuit_(circuit_, conditional_step_visitor, gate_step_visitor,
                       measurement_step_visitor, nop_step_visitor);

        // update the destructively measured qudits
        measured_d_.insert(
            std::next(measured_d_.begin(), static_cast<std::ptrdiff_t>(pos)), n,
            false);

        // update the non-destructively measured qudits
        measured_nd_.insert(
            std::next(measured_nd_.begin(), static_cast<std::ptrdiff_t>(pos)),
            n, false);

        // update (enlarge) the clean qudits vector
        clean_qudits_.insert(
            std::next(clean_qudits_.begin(), static_cast<std::ptrdiff_t>(pos)),
            n, true);

        return *this;
    }

    /**
     * \brief Adds \a n additional qudits after the last qudit
     *
     * \param n Number of qudits
     * \return Reference to the current instance
     */
    QCircuit& add_qudit(idx n = 1) { return add_qudit(n, nq_); }

    /**
     * \brief Adds \a n additional classical dits before dit \a pos
     *
     * \note Classical dits with indexes greater or equal than the newly
     * inserted ones have their indexes automatically incremented
     *
     * \param n Number of classical dits
     * \param pos Classical dit index
     * \return Reference to the current instance
     */
    QCircuit& add_dit(idx n, idx pos) {
        // EXCEPTION CHECKS
        if (pos > nc_) {
            throw exception::OutOfRange("qpp::QCircuit::add_dit()", "pos");
        }
        // END EXCEPTION CHECKS

        nc_ += n;

        // update cCTRL gate indexes
        auto gate_step_visitor = [&](internal::QCircuitGateStep& gate_step) {
            if (is_cCTRL(gate_step)) {
                for (idx& elem : gate_step.ctrl_.value()) {
                    if (elem >= pos) {
                        elem += n;
                    }
                }
            }
        };
        // update measurement indexes
        auto measurement_step_visitor =
            [&](internal::QCircuitMeasurementStep& measurement_step) {
                if (measurement_step.c_reg_ >= pos) {
                    measurement_step.c_reg_ += n;
                }
            };
        auto nop_step_visitor = [&](internal::QCircuitNOPStep&) {};
        auto conditional_step_visitor =
            [&](const internal::QCircuitConditionalStep&) {};

        visit_circuit_(circuit_, conditional_step_visitor, gate_step_visitor,
                       measurement_step_visitor, nop_step_visitor);

        // update (enlarge) the clean dits vector
        clean_dits_.insert(
            std::next(clean_dits_.begin(), static_cast<std::ptrdiff_t>(pos)), n,
            true);

        // update (enlarge) the measurement dits vector
        measurement_dits_.insert(std::next(measurement_dits_.begin(),
                                           static_cast<std::ptrdiff_t>(pos)),
                                 n, false);

        return *this;
    }

    /**
     * \brief Adds \a n additional classical dits after the last dit
     *
     * \param n Number of classical dits
     * \return Reference to the current instance
     */
    QCircuit& add_dit(idx n = 1) { return add_dit(n, nc_); }

    QCircuit& cond_if(std::function<bool(std::vector<idx>)> cond_func) {
        circuit_.emplace_back(internal::QCircuitConditionalStep{
            internal::QCircuitConditionalStep::Type::IF, cond_func});

        return *this;
    }

    QCircuit& cond_else() {
        circuit_.emplace_back(internal::QCircuitConditionalStep{
            internal::QCircuitConditionalStep::Type::ELSE});

        return *this;
    }

    QCircuit& cond_endif() {
        circuit_.emplace_back(internal::QCircuitConditionalStep{
            internal::QCircuitConditionalStep::Type::ENDIF});

        return *this;
    }

    /**
     * \brief Applies the single qudit gate \a U on single qudit \a i
     *
     * \param U Single qudit quantum gate
     * \param i Qudit index
     * \param name Optional gate name
     * \return Reference to the current instance
     */
    QCircuit& gate(const cmat& U, idx i,
                   std::optional<std::string> name = std::nullopt) {
        // EXCEPTION CHECKS
        std::string context{"Step " + std::to_string(get_step_count())};

        // check valid target
        if (i >= nq_) {
            throw exception::OutOfRange("qpp::QCircuit::gate()",
                                        context + ": i");
        }
        // check not measured before
        if (was_measured_d(i)) {
            throw exception::QuditAlreadyMeasured("qpp::QCircuit::gate()",
                                                  context + ": i");
        }

        // check square matrix for the gate
        if (!internal::check_square_mat(U)) {
            throw exception::MatrixNotSquare("qpp::QCircuit::gate()",
                                             context + ": U");
        }
        // check correct dimension
        if (static_cast<idx>(U.rows()) != d_) {
            throw exception::MatrixMismatchSubsys("qpp::QCircuit::gate()",
                                                  context + ": U");
        }
        // END EXCEPTION CHECKS

        if (!name.has_value()) {
            name = qpp::Gates::get_no_thread_local_instance().get_name(U);
        };
        std::size_t hashU = hash_eigen(U);
        add_hash_(hashU, U);

        circuit_.emplace_back(internal::QCircuitGateStep{
            internal::QCircuitGateStep::Type::SINGLE, hashU, std::nullopt,
            std::vector<idx>{i}, std::nullopt, name});

        ++gate_count_[hashU];

        clean_qudits_[i] = false;

        return *this;
    }

    /**
     * \brief Applies the two qudit gate \a U on qudits \a i and \a j
     *
     * \param U Two qudit quantum gate
     * \param i Qudit index
     * \param j Qudit index
     * \param name Optional gate name
     * \return Reference to the current instance
     */
    QCircuit& gate(const cmat& U, idx i, idx j,
                   std::optional<std::string> name = std::nullopt) {
        // EXCEPTION CHECKS
        std::string context{"Step " + std::to_string(get_step_count())};

        // check valid target
        if (i >= nq_ || j >= nq_ || i == j) {
            throw exception::OutOfRange("qpp::QCircuit::gate()",
                                        context + ": i/j");
        }
        if (was_measured_d(i) || was_measured_d(j)) {
            throw exception::QuditAlreadyMeasured("qpp::QCircuit::gate()",
                                                  context + ": i/j");
        }

        // check square matrix for the gate
        if (!internal::check_square_mat(U)) {
            throw exception::MatrixNotSquare("qpp::QCircuit::gate()",
                                             context + ": U");
        }
        // check correct dimension
        if (static_cast<idx>(U.rows()) != d_ * d_) {
            throw exception::MatrixMismatchSubsys("qpp::QCircuit::gate()",
                                                  context + ": U");
        }
        // END EXCEPTION CHECKS

        if (!name.has_value()) {
            name = qpp::Gates::get_no_thread_local_instance().get_name(U);
        }
        std::size_t hashU = hash_eigen(U);
        add_hash_(hashU, U);

        circuit_.emplace_back(internal::QCircuitGateStep{
            internal::QCircuitGateStep::Type::TWO, hashU, std::nullopt,
            std::vector<idx>{i, j}, std::nullopt, name});

        ++gate_count_[hashU];

        clean_qudits_[i] = clean_qudits_[j] = false;

        return *this;
    }

    /**
     * \brief Applies the three qudit gate \a U on qudits \a i, \a j and \a
     * k
     *
     * \param U Three qudit quantum gate
     * \param i Qudit index
     * \param j Qudit index
     * \param k Qudit index
     * \param name Optional gate name
     * \return Reference to the current instance
     */
    QCircuit& gate(const cmat& U, idx i, idx j, idx k,
                   std::optional<std::string> name = std::nullopt) {
        // EXCEPTION CHECKS
        std::string context{"Step " + std::to_string(get_step_count())};

        // check valid target
        if (i >= nq_ || j >= nq_ || k >= nq_ || (i == j) || (i == k) ||
            (j == k)) {
            throw exception::OutOfRange("qpp::QCircuit::gate()",
                                        context + ": i/j/k");
        }
        if (was_measured_d(i) || was_measured_d(j) || was_measured_d(k)) {
            throw exception::QuditAlreadyMeasured("qpp::QCircuit::gate()",
                                                  context + ": i/j/k");
        }

        // check square matrix for the gate
        if (!internal::check_square_mat(U)) {
            throw exception::MatrixNotSquare("qpp::QCircuit::gate()",
                                             context + ": U");
        }
        // check correct dimension
        if (static_cast<idx>(U.rows()) != d_ * d_ * d_) {
            throw exception::MatrixMismatchSubsys("qpp::QCircuit::gate()",
                                                  context + ": U");
        }
        // END EXCEPTION CHECKS

        if (!name.has_value()) {
            name = qpp::Gates::get_no_thread_local_instance().get_name(U);
        }
        std::size_t hashU = hash_eigen(U);
        add_hash_(hashU, U);

        circuit_.emplace_back(internal::QCircuitGateStep{
            internal::QCircuitGateStep::Type::THREE, hashU, std::nullopt,
            std::vector<idx>{i, j, k}, std::nullopt, name});

        ++gate_count_[hashU];

        clean_qudits_[i] = clean_qudits_[j] = clean_qudits_[k] = false;

        return *this;
    }

    /**
     * \brief Applies the single qudit gate \a U on every qudit listed in
     * \a target
     *
     * \param U Single qudit quantum gate
     * \param target Target qudit indexes; the gate \a U is applied on every one
     * of them
     * \param name Optional gate name
     * \return Reference to the current instance
     */
    QCircuit& gate_fan(const cmat& U, const std::vector<idx>& target,
                       std::optional<std::string> name = std::nullopt) {
        // EXCEPTION CHECKS
        std::string context{"Step " + std::to_string(get_step_count())};

        // check valid target
        if (target.empty()) {
            throw exception::ZeroSize("qpp::QCircuit::gate_fan()",
                                      context + ": target");
        }
        for (idx elem : target) {
            if (elem >= nq_) {
                throw exception::OutOfRange("qpp::QCircuit::gate_fan()",
                                            context + ": target");
            }
            // check target was not measured before
            if (was_measured_d(elem)) {
                throw exception::QuditAlreadyMeasured(
                    "qpp::QCircuit::gate_fan()", context + ": target");
            }
        }
        // check no duplicates target
        if (!internal::check_no_duplicates(target)) {
            throw exception::Duplicates("qpp::QCircuit::gate_fan()",
                                        context + ": target");
        }

        // check square matrix for the gate
        if (!internal::check_square_mat(U)) {
            throw exception::MatrixNotSquare("qpp::QCircuit::gate_fan()",
                                             context + ": U");
        }
        // check correct dimension
        if (static_cast<idx>(U.rows()) != d_) {
            throw exception::MatrixMismatchSubsys("qpp::QCircuit::gate_fan()",
                                                  context + ": U");
        }
        // END EXCEPTION CHECKS

        if (!name.has_value()) {
            name = qpp::Gates::get_no_thread_local_instance().get_name(U);
        }
        std::size_t hashU = hash_eigen(U);
        add_hash_(hashU, U);

        circuit_.emplace_back(internal::QCircuitGateStep{
            internal::QCircuitGateStep::Type::FAN, hashU, std::nullopt, target,
            std::nullopt, name});

        gate_count_[hashU] += target.size();

        for (idx elem : target) {
            clean_qudits_[elem] = false;
        }

        return *this;
    }

    // std::initializer_list overload, avoids ambiguity for 2-element lists,
    // see
    // https://stackoverflow.com/questions/26750039/ambiguity-when-using-initializer-list-as-parameter
    /**
     * \brief Applies the single qudit gate \a U on every qudit listed in
     * \a target
     *
     * \param U Single qudit quantum gate
     * \param target Target qudit indexes; the gate \a U is applied on every
     * one of them
     * \param name Optional gate name
     * \return Reference to the current instance
     */
    QCircuit& gate_fan(const cmat& U, const std::initializer_list<idx>& target,
                       std::optional<std::string> name = std::nullopt) {
        return gate_fan(U, std::vector<idx>(target), std::move(name));
    }

    /**
     * \brief Applies the single qudit gate \a U on all of the remaining
     * non-measured qudits
     *
     * \param U Single qudit quantum gate
     * \param name Optional gate name
     * \return Reference to the current instance
     */
    QCircuit& gate_fan(const cmat& U,
                       std::optional<std::string> name = std::nullopt) {
        // EXCEPTION CHECKS
        std::string context{"Step " + std::to_string(get_step_count())};

        // check non-empty target
        std::vector<idx> target = get_non_measured_d();
        if (target.empty()) {
            throw exception::QuditAlreadyMeasured(
                "qpp::QCircuit::gate_fan()",
                context + ": all qudits have been already measured");
        }
        // END EXCEPTION CHECKS

        if (!name.has_value()) {
            name = qpp::Gates::get_no_thread_local_instance().get_name(U);
        }

        return gate_fan(U, target, name);
    }

    /**
     * \brief Jointly applies the multiple-qudit gate \a U on the qudit
     * indexes specified by \a target
     *
     * \param U Multiple qudit quantum gate
     * \param target Subsystem indexes where the gate \a U is applied
     * \param name Optional gate name
     * \return Reference to the current instance
     */
    QCircuit& gate(const cmat& U, const std::vector<idx>& target,
                   std::optional<std::string> name = std::nullopt) {
        // EXCEPTION CHECKS
        idx n = static_cast<idx>(target.size());
        idx D = internal::safe_pow(d_, n);

        std::string context{"Step " + std::to_string(get_step_count())};

        // check valid target
        if (target.empty()) {
            throw exception::ZeroSize("qpp::QCircuit::gate()",
                                      context + ": target");
        }
        for (idx elem : target) {
            if (elem >= nq_) {
                throw exception::OutOfRange("qpp::QCircuit::gate()",
                                            context + ": target");
            }
            // check target was not measured before
            if (was_measured_d(elem)) {
                throw exception::QuditAlreadyMeasured("qpp::QCircuit::gate()",
                                                      context + ": target");
            }
        }
        // check no duplicates target
        if (!internal::check_no_duplicates(target)) {
            throw exception::Duplicates("qpp::QCircuit::gate()",
                                        context + ": target");
        }

        // check square matrix for the gate
        if (!internal::check_square_mat(U)) {
            throw exception::MatrixNotSquare("qpp::QCircuit::gate()",
                                             context + ": U");
        }
        // check correct dimension
        if (static_cast<idx>(U.rows()) != D) {
            throw exception::MatrixMismatchSubsys("qpp::QCircuit::gate()",
                                                  context + ": U");
        }
        // END EXCEPTION CHECKS

        if (!name.has_value()) {
            name = qpp::Gates::get_no_thread_local_instance().get_name(U);
        }
        std::size_t hashU = hash_eigen(U);
        add_hash_(hashU, U);

        circuit_.emplace_back(internal::QCircuitGateStep{
            internal::QCircuitGateStep::Type::JOINT, hashU, std::nullopt,
            target, std::nullopt, name});

        ++gate_count_[hashU];

        for (idx elem : target) {
            clean_qudits_[elem] = false;
        }

        return *this;
    }

    /**
     * \brief Applies the quantum Fourier transform (as a series of gates)
     * on the qudit indexes specified by \a target
     *
     * \param target Subsystem indexes where the quantum Fourier transform
     * is applied
     * \param swap Swaps the qubits at the end (true by default)
     * \return Reference to the current instance
     */
    QCircuit& QFT(const std::vector<idx>& target, bool swap = true) {
        // EXCEPTION CHECKS
        std::string context{"Step " + std::to_string(get_step_count())};

        // check valid target
        if (target.empty()) {
            throw exception::ZeroSize("qpp::QCircuit::QFT()",
                                      context + ": target");
        }
        for (idx elem : target) {
            if (elem >= nq_) {
                throw exception::OutOfRange("qpp::QCircuit::QFT()",
                                            context + ": target");
            }
            // check target was not measured before
            if (was_measured_d(elem)) {
                throw exception::QuditAlreadyMeasured("qpp::QCircuit::QFT()",
                                                      context + ": target");
            }
        }
        // check no duplicates target
        if (!internal::check_no_duplicates(target)) {
            throw exception::Duplicates("qpp::QCircuit::QFT()",
                                        context + ": target");
        }
        // END EXCEPTION CHECKS

        idx n_subsys = target.size();
        if (d_ == 2) // qubits
        {
            for (idx i = 0; i < n_subsys; ++i) {
                // apply Hadamard on qubit i
                gate(Gates::get_no_thread_local_instance().H, target[i]);
                // apply controlled rotations
                for (idx j = 2; j <= n_subsys - i; ++j) {
                    // construct Rj
                    cmat Rj(2, 2);
                    Rj << 1, 0, 0,
                        std::exp(static_cast<cplx::value_type>(2.0 * pi) * 1_i /
                                 static_cast<cplx::value_type>(std::pow(2, j)));
                    CTRL(Rj, target[i + j - 1], target[i], {},
                         "CTRL-R" + std::to_string(j));
                }
            }
            if (swap) {
                // we have the qubits in reversed order, we must swap them
                for (idx i = 0; i < n_subsys / 2; ++i) {
                    gate(Gates::get_no_thread_local_instance().SWAP, target[i],
                         target[n_subsys - i - 1]);
                }
            }

        } else { // qudits
            for (idx i = 0; i < n_subsys; ++i) {
                // apply qudit Fourier on qudit i
                gate(Gates::get_no_thread_local_instance().Fd(d_), target[i],
                     "Fd");
                // apply controlled rotations
                for (idx j = 2; j <= n_subsys - i; ++j) {
                    // construct Rj
                    cmat Rj = cmat::Zero(d_, d_);
                    for (idx m = 0; m < d_; ++m) {
                        Rj(m, m) = std::exp(
                            static_cast<cplx::value_type>(2.0 * pi * m) * 1_i /
                            static_cast<cplx::value_type>(std::pow(d_, j)));
                    }
                    CTRL(Rj, target[i + j - 1], target[i], {},
                         "CTRL-R" + std::to_string(j) + "d");
                }
            }
            if (swap) {
                // we have the qudits in reversed order, we must swap them
                for (idx i = 0; i < n_subsys / 2; ++i) {
                    gate(Gates::get_no_thread_local_instance().SWAPd(d_),
                         target[i], target[n_subsys - i - 1], "SWAPd");
                }
            }
        }

        return *this;
    }

    // std::initializer_list overload, avoids ambiguity for {idx} -> bool
    /**
     * \brief Applies the quantum Fourier transform (as a series of gates)
     * on the qudit indexes specified by \a target
     *
     * \param target Subsystem indexes where the quantum Fourier transform
     * is applied
     * \param swap Swaps the qubits at the end (true by default)
     * \return Reference to the current instance
     */
    QCircuit& QFT(const std::initializer_list<idx>& target, bool swap = true) {
        return QFT(std::vector<idx>(target), swap);
    }

    /**
     * \brief Applies the quantum Fourier transform (as a series of gates)
     * on all of remaining non-measured qudits
     *
     * \param swap Swaps the qubits at the end (true by default)
     * \return Reference to the current instance
     */
    QCircuit& QFT(bool swap = true) { return QFT(get_non_measured_d(), swap); }

    /**
     * \brief Applies the inverse quantum Fourier transform (as a series of
     * gates) on the qudit indexes specified by \a target
     *
     * \param target Subsystem indexes where the inverse quantum Fourier
     * transform is applied
     * \param swap Swaps the qubits at the end (true by default)
     * \return Reference to the current instance
     */
    QCircuit& TFQ(const std::vector<idx>& target,
                  [[maybe_unused]] bool swap = true) {
        // EXCEPTION CHECKS
        std::string context{"Step " + std::to_string(get_step_count())};

        // check valid target
        if (target.empty()) {
            throw exception::ZeroSize("qpp::QCircuit::TFQ()",
                                      context + ": target");
        }
        for (idx elem : target) {
            if (elem >= nq_) {
                throw exception::OutOfRange("qpp::QCircuit::TFQ()",
                                            context + ": target");
            }
            // check target was not measured before
            if (was_measured_d(elem)) {
                throw exception::QuditAlreadyMeasured("qpp::QCircuit::TFQ()",
                                                      context + ": target");
            }
        }
        // check no duplicates target
        if (!internal::check_no_duplicates(target)) {
            throw exception::Duplicates("qpp::QCircuit::TFQ()",
                                        context + ": target");
        }
        // END EXCEPTION CHECKS

        idx n_subsys = target.size();
        if (d_ == 2) // qubits
        {
            if (swap) {
                // we have the qubits in reversed order, we must swap them
                for (idx i = n_subsys / 2; i-- > 0;) {
                    gate(Gates::get_no_thread_local_instance().SWAP, target[i],
                         target[n_subsys - i - 1]);
                }
            }
            for (idx i = n_subsys; i-- > 0;) {
                // apply controlled rotations
                for (idx j = n_subsys - i + 1; j-- > 2;) {
                    // construct Rj
                    cmat Rj(2, 2);
                    Rj << 1, 0, 0,
                        std::exp(static_cast<cplx::value_type>(-2.0 * pi) *
                                 1_i /
                                 static_cast<cplx::value_type>(std::pow(2, j)));
                    CTRL(Rj, target[i + j - 1], target[i], {},
                         "CTRL-R" + std::to_string(j) + "+");
                }
                // apply Hadamard on qubit i
                gate(Gates::get_no_thread_local_instance().H, target[i]);
            }
        } else { // qudits
            if (swap) {
                // we have the qudits in reversed order, we must swap them
                for (idx i = n_subsys / 2; i-- > 0;) {
                    gate(Gates::get_no_thread_local_instance().SWAPd(d_),
                         target[i], target[n_subsys - i - 1], "SWAPd");
                }
            }
            for (idx i = n_subsys; i-- > 0;) {
                // apply controlled rotations
                for (idx j = n_subsys - i + 1; j-- > 2;) {
                    // construct Rj
                    cmat Rj = cmat::Zero(d_, d_);
                    for (idx m = 0; m < d_; ++m) {
                        Rj(m, m) = std::exp(
                            static_cast<cplx::value_type>(-2.0 * pi * m) * 1_i /
                            static_cast<cplx::value_type>(std::pow(d_, j)));
                    }
                    CTRL(Rj, target[i + j - 1], target[i], {},
                         "CTRL-R" + std::to_string(j) + "d+");
                }
                // apply adjoint qudit Fourier on qudit i
                gate(qpp::adjoint(Gates::get_no_thread_local_instance().Fd(d_)),
                     target[i], "Fd+");
            }
        }

        return *this;
    }

    // std::initializer_list overload, avoids ambiguity for {idx} -> bool
    /**
     * \brief Applies the inverse quantum Fourier transform (as a series of
     * gates) on the qudit indexes specified by \a target
     *
     * \param target Subsystem indexes where the inverse quantum Fourier
     * transform is applied
     * \param swap Swaps the qubits at the end (true by default)
     * \return Reference to the current instance
     */
    QCircuit& TFQ(const std::initializer_list<idx>& target, bool swap = true) {
        return TFQ(std::vector<idx>(target), swap);
    }

    /**
     * \brief Applies the inverse quantum Fourier transform (as a series of
     * gates) on all of remaining non-measured qudits
     *
     * \param swap Swaps the qubits at the end (true by default)
     * \return Reference to the current instance
     */
    QCircuit& TFQ(bool swap = true) { return TFQ(get_non_measured_d(), swap); }

    // single ctrl multiple targets
    /**
     * \brief Applies the single qudit controlled gate \a U with control
     * qudit \a ctrl on every qudit listed in \a target, i.e.,
     * CTRL-U-U-...-U
     *
     * \param U Single qudit quantum gate
     * \param ctrl Control qudit index
     * \param target Target qudit indexes; the gate \a U is applied on every
     * one of them depending on the values of the control qudits
     * \param shift Optional, performs the control as if the \a ctrl qudit state
     * was \f$X\f$-incremented by \a shift
     * \param name Optional gate name
     * \return Reference to the current instance
     */
    QCircuit& CTRL_fan(const cmat& U, idx ctrl, const std::vector<idx>& target,
                       std::optional<idx> shift = std::nullopt,
                       std::optional<std::string> name = std::nullopt) {
        // EXCEPTION CHECKS
        std::string context{"Step " + std::to_string(get_step_count())};

        // check valid ctrl
        if (ctrl >= nq_) {
            throw exception::OutOfRange("qpp::QCircuit::CTRL_fan()",
                                        context + ": ctrl");
        }
        if (was_measured_d(ctrl)) {
            throw exception::QuditAlreadyMeasured("qpp::QCircuit::CTRL_fan()",
                                                  context + ": ctrl");
        }

        // check valid target
        if (target.empty()) {
            throw exception::ZeroSize("qpp::QCircuit::CTRL_fan()",
                                      context + ": target");
        }
        for (idx elem : target) {
            if (elem >= nq_) {
                throw exception::OutOfRange("qpp::QCircuit::CTRL_fan()",
                                            context + ": target");
            }
            // check target was not measured before
            if (was_measured_d(elem)) {
                throw exception::QuditAlreadyMeasured(
                    "qpp::QCircuit::CTRL_fan()", context + ": target");
            }
        }
        // check no duplicates target
        if (!internal::check_no_duplicates(target)) {
            throw exception::Duplicates("qpp::QCircuit::CTRL_fan()",
                                        context + ": target");
        }

        // check that ctrl and target don't overlap
        for (idx elem : target) {
            if (elem == ctrl) {
                throw exception::OutOfRange("qpp::QCircuit::CTRL_fan()",
                                            context + ": ctrl/target");
            }
        }

        // check square matrix for the gate
        if (!internal::check_square_mat(U)) {
            throw exception::MatrixNotSquare("qpp::QCircuit::CTRL_fan()",
                                             context + ": U");
        }
        // check correct dimension
        if (static_cast<idx>(U.rows()) != d_) {
            throw exception::MatrixMismatchSubsys("qpp::QCircuit::CTRL_fan()",
                                                  context + ": U");
        }

        // check shift
        if (shift.has_value() && shift.value() >= d_) {
            throw exception::OutOfRange("qpp::QCircuit::CTRL_fan()",
                                        context + ": shift");
        }
        // END EXCEPTION CHECKS

        if (!name.has_value()) {
            auto gate_name =
                qpp::Gates::get_no_thread_local_instance().get_name(U);
            name = gate_name.has_value() ? "CTRL_fan-" + gate_name.value()
                                         : "CTRL_fan";
        }
        std::size_t hashU = hash_eigen(U);
        add_hash_(hashU, U);

        std::optional<std::vector<idx>> shift_vec;
        if (shift.has_value()) {
            shift_vec = std::vector<idx>{shift.value()};
        }
        circuit_.emplace_back(internal::QCircuitGateStep{
            internal::QCircuitGateStep::Type::CTRL_FAN, hashU,
            std::vector<idx>{ctrl}, target, shift_vec, name});

        gate_count_[hashU] += target.size();

        clean_qudits_[ctrl] = false;
        for (idx elem : target) {
            clean_qudits_[elem] = false;
        }

        return *this;
    }

    // multiple ctrl multiple targets
    /**
     * \brief Applies the single qudit controlled gate \a U with multiple
     * control qudits listed in \a ctrl on every qudit listed in \a target,
     * i.e., CTRL-CTRL-...-CTRL-U-U-...-U
     *
     * \param U Single qudit quantum gate
     * \param ctrl Control qudit indexes
     * \param target Target qudit indexes; the gate \a U is applied on every
     * one of them depending on the values of the control qudits
     * \param shift Performs the control as if the \a ctrl qudit states were
     * \f$X\f$-incremented component-wise by \a shift. If non-empty
     * (default), the size of \a shift must be the same as the size of \a
     * ctrl.
     * \param name Optional gate name
     * \return Reference to the current instance
     */
    QCircuit& CTRL_fan(const cmat& U, const std::vector<idx>& ctrl,
                       const std::vector<idx>& target,
                       std::optional<std::vector<idx>> shift = std::nullopt,
                       std::optional<std::string> name = std::nullopt) {
        // EXCEPTION CHECKS
        std::string context{"Step " + std::to_string(get_step_count())};

        // check valid ctrl
        if (ctrl.empty()) {
            throw exception::ZeroSize("qpp::QCircuit::CTRL_fan()",
                                      context + ": ctrl");
        }
        for (idx elem : ctrl) {
            if (elem >= nq_) {
                throw exception::OutOfRange("qpp::QCircuit::CTRL_fan()",
                                            context + ": ctrl");
            }
            // check ctrl was not measured before
            if (was_measured_d(elem)) {
                throw exception::QuditAlreadyMeasured(
                    "qpp::QCircuit::CTRL_fan()", context + ": ctrl");
            }
        }
        // check no duplicates ctrl
        if (!internal::check_no_duplicates(ctrl)) {
            throw exception::Duplicates("qpp::QCircuit::CTRL_fan()",
                                        context + ": ctrl");
        }

        // check valid target
        if (target.empty()) {
            throw exception::ZeroSize("qpp::QCircuit::CTRL_fan()",
                                      context + ": target");
        }
        for (idx elem : target) {
            if (elem >= nq_) {
                throw exception::OutOfRange("qpp::QCircuit::CTRL_fan()",
                                            context + ": target");
            }
            // check target was not measured before
            if (was_measured_d(elem)) {
                throw exception::QuditAlreadyMeasured(
                    "qpp::QCircuit::CTRL_fan()", context + ": target");
            }
        }
        // check no duplicates target
        if (!internal::check_no_duplicates(target)) {
            throw exception::Duplicates("qpp::QCircuit::CTRL_fan()",
                                        context + ": target");
        }

        // check that ctrl and target don't overlap
        for (idx elem_ctrl : ctrl) {
            for (idx elem_target : target) {
                if (elem_ctrl == elem_target) {
                    throw exception::OutOfRange("qpp::QCircuit::CTRL_fan()",
                                                context + ": ctrl/target");
                }
            }
        }

        // check square matrix for the gate
        if (!internal::check_square_mat(U)) {
            throw exception::MatrixNotSquare("qpp::QCircuit::CTRL_fan()",
                                             context + ": U");
        }
        // check correct dimension
        if (static_cast<idx>(U.rows()) != d_) {
            throw exception::MatrixMismatchSubsys("qpp::QCircuit::CTRL_fan()",
                                                  context + ": U");
        }

        // check shift
        if (shift.has_value() && (shift.value().size() != ctrl.size())) {
            throw exception::SizeMismatch("qpp::QCircuit::CTRL_fan()",
                                          context + ": ctrl/shift");
        }
        if (shift.has_value()) {
            for (idx elem : shift.value()) {
                if (elem >= d_) {
                    throw exception::OutOfRange("qpp::QCircuit::CTRL_fan()",
                                                context + ": shift");
                }
            }
        }
        // END EXCEPTION CHECKS

        if (!name.has_value()) {
            auto gate_name =
                qpp::Gates::get_no_thread_local_instance().get_name(U);
            name = gate_name.has_value() ? "CTRL_fan-" + gate_name.value()
                                         : "CTRL_fan";
        }
        std::size_t hashU = hash_eigen(U);
        add_hash_(hashU, U);

        circuit_.emplace_back(internal::QCircuitGateStep{
            internal::QCircuitGateStep::Type::CTRL_FAN, hashU, ctrl, target,
            shift, name});

        gate_count_[hashU] += target.size();

        for (idx elem : ctrl) {
            clean_qudits_[elem] = false;
        }
        for (idx elem : target) {
            clean_qudits_[elem] = false;
        }

        return *this;
    }

    // multiple ctrl joint target
    /**
     * \brief Jointly applies the multiple-qudit controlled gate \a U with
     * multiple control qudits listed in \a ctrl on the qudit indexes
     * specified by \a target, i.e., CTRL-CTRL-...-CTRL-U_{joint}
     *
     * \param U Multiple-qudit quantum gate
     * \param ctrl Control qudit indexes
     * \param target Target qudit indexes where the gate \a U is applied
     * depending on the values of the control qudits
     * \param shift Optional, performs the control as if the \a ctrl qudit
     * states were \f$X\f$-incremented component-wise by \a shift. If present,
     * the size of \a shift must be the same as the size of \a ctrl.
     * \param name Optional gate name
     * \return Reference to the current instance
     */
    QCircuit& CTRL(const cmat& U, const std::vector<idx>& ctrl,
                   const std::vector<idx>& target,
                   std::optional<std::vector<idx>> shift = std::nullopt,
                   std::optional<std::string> name = std::nullopt) {
        // EXCEPTION CHECKS
        std::string context{"Step " + std::to_string(get_step_count())};

        idx D_target = internal::safe_pow(d_, target.size());

        // check valid ctrl
        if (ctrl.empty()) {
            throw exception::ZeroSize("qpp::QCircuit::CTRL()",
                                      context + ": ctrl");
        }
        for (idx elem : ctrl) {
            if (elem >= nq_) {
                throw exception::OutOfRange("qpp::QCircuit::CTRL()",
                                            context + ": ctrl");
            }
            // check ctrl was not measured before
            if (was_measured_d(elem)) {
                throw exception::QuditAlreadyMeasured("qpp::QCircuit::CTRL()",
                                                      context + ": ctrl");
            }
        }
        // check no duplicates ctrl
        if (!internal::check_no_duplicates(ctrl)) {
            throw exception::Duplicates("qpp::QCircuit::CTRL()",
                                        context + ": ctrl");
        }

        // check valid target
        if (target.empty()) {
            throw exception::ZeroSize("qpp::QCircuit::CTRL()",
                                      context + ": target");
        }
        for (idx elem : target) {
            if (elem >= nq_) {
                throw exception::OutOfRange("qpp::QCircuit::CTRL()",
                                            context + ": target");
            }
            // check target was not measured before
            if (was_measured_d(elem)) {
                throw exception::QuditAlreadyMeasured("qpp::QCircuit::CTRL()",
                                                      context + ": target");
            }
        }

        // check that ctrl and target don't overlap
        for (idx elem_ctrl : ctrl) {
            for (idx elem_target : target) {
                if (elem_ctrl == elem_target) {
                    throw exception::OutOfRange("qpp::QCircuit::CTRL()",
                                                context + ": ctrl/target");
                }
            }
        }

        // check square matrix for the gate
        if (!internal::check_square_mat(U)) {
            throw exception::MatrixNotSquare("qpp::QCircuit::CTRL()",
                                             context + ": U");
        }
        // check correct dimension
        if (static_cast<idx>(U.rows()) != D_target) {
            throw exception::MatrixMismatchSubsys("qpp::QCircuit::CTRL()",
                                                  context + ": U");
        }

        // check shift
        if (shift.has_value() && (shift.value().size() != ctrl.size())) {
            throw exception::SizeMismatch("qpp::QCircuit::CTRL()",
                                          context + ": ctrl/shift");
        }
        if (shift.has_value()) {
            for (idx elem : shift.value()) {
                if (elem >= d_) {
                    throw exception::OutOfRange("qpp::QCircuit::CTRL()",
                                                context + ": shift");
                }
            }
        }
        // END EXCEPTION CHECKS

        if (!name.has_value()) {
            auto gate_name =
                qpp::Gates::get_no_thread_local_instance().get_name(U);
            name = gate_name.has_value() ? "CTRL-" + gate_name.value() : "CTRL";
        }
        std::size_t hashU = hash_eigen(U);
        add_hash_(hashU, U);

        circuit_.emplace_back(
            internal::QCircuitGateStep{internal::QCircuitGateStep::Type::CTRL,
                                       hashU, ctrl, target, shift, name});

        ++gate_count_[hashU];

        for (idx elem : ctrl) {
            clean_qudits_[elem] = false;
        }
        for (idx elem : target) {
            clean_qudits_[elem] = false;
        }

        return *this;
    }

    // multiple ctrl single target
    /**
     * \brief Applies the multiple-qudit controlled gate \a U with
     * multiple control qudits listed in \a ctrl on the target qudit
     * specified by \a target, i.e., CTRL-CTRL-...-CTRL-U
     *
     * \param U Single qudit quantum gate
     * \param ctrl Control qudit indexes
     * \param target Target qudit index
     * \param shift Optional, performs the control as if the \a ctrl qudit
     * states were \f$X\f$-incremented component-wise by \a shift. If present,
     * the size of \a shift must be the same as the size of \a ctrl.
     * \param name Optional gate name
     * \return Reference to the current instance
     */
    QCircuit& CTRL(const cmat& U, const std::vector<idx>& ctrl, idx target,
                   std::optional<std::vector<idx>> shift = std::nullopt,
                   std::optional<std::string> name = std::nullopt) {
        // EXCEPTION CHECKS
        std::string context{"Step " + std::to_string(get_step_count())};

        // check valid ctrl
        if (ctrl.empty()) {
            throw exception::ZeroSize("qpp::QCircuit::CTRL()",
                                      context + ": ctrl");
        }
        for (idx elem : ctrl) {
            if (elem >= nq_) {
                throw exception::OutOfRange("qpp::QCircuit::CTRL()",
                                            context + ": ctrl");
            }
            // check ctrl was not measured before
            if (was_measured_d(elem)) {
                throw exception::QuditAlreadyMeasured("qpp::QCircuit::CTRL()",
                                                      context + ": ctrl");
            }
        }
        // check no duplicates ctrl
        if (!internal::check_no_duplicates(ctrl)) {
            throw exception::Duplicates("qpp::QCircuit::CTRL()",
                                        context + ": ctrl");
        }

        // check valid target
        if (target >= nq_) {
            throw exception::OutOfRange("qpp::QCircuit::CTRL()",
                                        context + ": target");
        }
        if (was_measured_d(target)) {
            throw exception::QuditAlreadyMeasured("qpp::QCircuit::CTRL()",
                                                  context + ": target");
        }

        // check that ctrl and target don't overlap
        for (idx elem : ctrl) {
            if (elem == target) {
                throw exception::OutOfRange("qpp::QCircuit::CTRL()",
                                            context + ": ctrl/target");
            }
        }

        // check square matrix for the gate
        if (!internal::check_square_mat(U)) {
            throw exception::MatrixNotSquare("qpp::QCircuit::CTRL()",
                                             context + ": U");
        }
        // check correct dimension
        if (static_cast<idx>(U.rows()) != d_) {
            throw exception::MatrixMismatchSubsys("qpp::QCircuit::CTRL()",
                                                  context + ": U");
        }

        // check shift
        if (shift.has_value() && (shift.value().size() != ctrl.size())) {
            throw exception::SizeMismatch("qpp::QCircuit::CTRL()",
                                          context + ": ctrl/shift");
        }
        if (shift.has_value()) {
            for (idx elem : shift.value()) {
                if (elem >= d_) {
                    throw exception::OutOfRange("qpp::QCircuit::CTRL()",
                                                context + ": shift");
                }
            }
        }
        // END EXCEPTION CHECKS

        if (!name.has_value()) {
            auto gate_name =
                qpp::Gates::get_no_thread_local_instance().get_name(U);
            name = gate_name.has_value() ? "CTRL-" + gate_name.value() : "CTRL";
        }
        std::size_t hashU = hash_eigen(U);
        add_hash_(hashU, U);

        circuit_.emplace_back(internal::QCircuitGateStep{
            internal::QCircuitGateStep::Type::CTRL, hashU, ctrl,
            std::vector<idx>{target}, shift, name});

        ++gate_count_[hashU];

        for (idx elem : ctrl) {
            clean_qudits_[elem] = false;
        }
        clean_qudits_[target] = false;

        return *this;
    }

    // single ctrl joint target
    /**
     * \brief Jointly applies the single qudit controlled gate \a U with
     * control qudit \a ctrl on the qudit indexes specified by \a target,
     * i.e., CTRL-U_{joint}
     *
     * \param U Multiple-qudit quantum gate
     * \param ctrl Control qudit index
     * \param target Target qudit indexes; the gate \a U is applied on every
     * one of them depending on the values of the control qudits
     * \param shift Optional, performs the control as if the \a ctrl qudit state
     * was \f$X\f$-incremented by \a shift
     * \param name Optional gate name
     * \return Reference to the current instance
     */
    QCircuit& CTRL(const cmat& U, idx ctrl, const std::vector<idx>& target,
                   std::optional<idx> shift = std::nullopt,
                   std::optional<std::string> name = std::nullopt) {
        // EXCEPTION CHECKS
        std::string context{"Step " + std::to_string(get_step_count())};

        idx D_target = internal::safe_pow(d_, target.size());

        // check valid ctrl
        if (ctrl >= nq_) {
            throw exception::OutOfRange("qpp::QCircuit::CTRL()",
                                        context + ": ctrl");
        }
        if (was_measured_d(ctrl)) {
            throw exception::QuditAlreadyMeasured("qpp::QCircuit::CTRL()",
                                                  context + ": ctrl");
        }

        // check valid target
        if (target.empty()) {
            throw exception::ZeroSize("qpp::QCircuit::CTRL()",
                                      context + ": target");
        }
        for (idx elem : target) {
            if (elem >= nq_) {
                throw exception::OutOfRange("qpp::QCircuit::CTRL()",
                                            context + ": target");
            }
            // check target was not measured before
            if (was_measured_d(elem)) {
                throw exception::QuditAlreadyMeasured("qpp::QCircuit::CTRL()",
                                                      context + ": target");
            }
        }
        // check no duplicates target
        if (!internal::check_no_duplicates(target)) {
            throw exception::Duplicates("qpp::QCircuit::CTRL()",
                                        context + ": target");
        }

        // check that ctrl and target don't overlap
        for (idx elem : target) {
            if (elem == ctrl) {
                throw exception::OutOfRange("qpp::QCircuit::CTRL()",
                                            context + ": ctrl/target");
            }
        }

        // check square matrix for the gate
        if (!internal::check_square_mat(U)) {
            throw exception::MatrixNotSquare("qpp::QCircuit::CTRL()",
                                             context + ": U");
        }
        // check correct dimension
        if (static_cast<idx>(U.rows()) != D_target) {
            throw exception::MatrixMismatchSubsys("qpp::QCircuit::CTRL()",
                                                  context + ": U");
        }

        // check shift
        if (shift.has_value() && shift.value() >= d_) {
            throw exception::OutOfRange("qpp::QCircuit::CTRL()",
                                        context + ": shift");
        }
        // END EXCEPTION CHECKS

        if (!name.has_value()) {
            auto gate_name =
                qpp::Gates::get_no_thread_local_instance().get_name(U);
            name = gate_name.has_value() ? "CTRL-" + gate_name.value() : "CTRL";
        }
        std::size_t hashU = hash_eigen(U);
        add_hash_(hashU, U);

        std::optional<std::vector<idx>> shift_vec;
        if (shift.has_value()) {
            shift_vec = std::vector<idx>{shift.value()};
        }
        circuit_.emplace_back(internal::QCircuitGateStep{
            internal::QCircuitGateStep::Type::CTRL, hashU,
            std::vector<idx>{ctrl}, target, shift_vec, name});

        ++gate_count_[hashU];

        clean_qudits_[ctrl] = false;
        for (idx elem : target) {
            clean_qudits_[elem] = false;
        }

        return *this;
    }

    // single ctrl single target
    /**
     * \brief Applies the single qudit controlled gate \a U with control
     * qudit \a ctrl and target qudit \a target, i.e., CTRL-U
     *
     * \param U Single qudit quantum gate
     * \param ctrl Control qudit index
     * \param target Target qudit index
     * \param shift Optional, performs the control as if the \a ctrl qudit state
     * was \f$X\f$-incremented by \a shift
     * \param name Optional gate name
     * \return Reference to the current instance
     */
    QCircuit& CTRL(const cmat& U, idx ctrl, idx target,
                   std::optional<idx> shift = std::nullopt,
                   std::optional<std::string> name = std::nullopt) {
        // EXCEPTION CHECKS
        std::string context{"Step " + std::to_string(get_step_count())};

        // check valid ctrl and target
        if (ctrl >= nq_ || target >= nq_ || ctrl == target) {
            throw exception::OutOfRange("qpp::QCircuit::CTRL()",
                                        context + ": ctrl/target");
        }
        if (was_measured_d(ctrl) || was_measured_d(target)) {
            throw exception::QuditAlreadyMeasured("qpp::QCircuit::CTRL()",
                                                  context + ": ctrl/target");
        }

        // check square matrix for the gate
        if (!internal::check_square_mat(U)) {
            throw exception::MatrixNotSquare("qpp::QCircuit::CTRL()",
                                             context + ": U");
        }
        // check correct dimension
        if (static_cast<idx>(U.rows()) != d_) {
            throw exception::MatrixMismatchSubsys("qpp::QCircuit::CTRL()",
                                                  context + ": U");
        }

        // check shift
        if (shift.has_value() && shift.value() >= d_) {
            throw exception::OutOfRange("qpp::QCircuit::CTRL()",
                                        context + ": shift");
        }
        // END EXCEPTION CHECKS

        if (!name.has_value()) {
            auto gate_name =
                qpp::Gates::get_no_thread_local_instance().get_name(U);
            name = gate_name.has_value() ? "CTRL-" + gate_name.value() : "CTRL";
        }
        std::size_t hashU = hash_eigen(U);
        add_hash_(hashU, U);

        std::optional<std::vector<idx>> shift_vec;
        if (shift.has_value()) {
            shift_vec = std::vector<idx>{shift.value()};
        }
        circuit_.emplace_back(internal::QCircuitGateStep{
            internal::QCircuitGateStep::Type::CTRL, hashU,
            std::vector<idx>{ctrl}, std::vector<idx>{target}, shift_vec, name});

        ++gate_count_[hashU];

        clean_qudits_[ctrl] = false;
        clean_qudits_[target] = false;

        return *this;
    }

    // single cctrl multiple targets
    /**
     * \brief Applies the single qudit controlled gate \a U with classical
     * control dit \a ctrl on every qudit listed in \a target, i.e.,
     * cCTRL-U-U-...-U
     *
     * \param U Single qudit quantum gate
     * \param ctrl_dit Classical control dit index
     * \param target Target qudit indexes; the gate \a U is applied on every
     * one of them depending on the values of the classical control dits
     * \param shift Optional, performs the control as if the \a ctrl_dit
     * classical dit was incremented by \a shift
     * \param name Optional gate name
     * \return Reference to the current instance
     */
    QCircuit& cCTRL_fan(const cmat& U, idx ctrl_dit,
                        const std::vector<idx>& target,
                        std::optional<idx> shift = std::nullopt,
                        std::optional<std::string> name = std::nullopt) {
        // EXCEPTION CHECKS
        std::string context{"Step " + std::to_string(get_step_count())};

        // check valid ctrl_dit
        if (ctrl_dit >= nc_) {
            throw exception::OutOfRange("qpp::QCircuit::cCTRL_fan()",
                                        context + ": ctrl_dit");
        }

        // check valid target
        if (target.empty()) {
            throw exception::ZeroSize("qpp::QCircuit::cCTRL_fan()",
                                      context + ": target");
        }
        for (idx elem : target) {
            if (elem >= nq_) {
                throw exception::OutOfRange("qpp::QCircuit::cCTRL_fan()",
                                            context + ": target");
            }
            // check target was not measured before
            if (was_measured_d(elem)) {
                throw exception::QuditAlreadyMeasured(
                    "qpp::QCircuit::cCTRL_fan()", context + ": target");
            }
        }
        // check no duplicates target
        if (!internal::check_no_duplicates(target)) {
            throw exception::Duplicates("qpp::QCircuit::cCTRL_fan()",
                                        context + ": target");
        }

        // check square matrix for the gate
        if (!internal::check_square_mat(U)) {
            throw exception::MatrixNotSquare("qpp::QCircuit::cCTRL_fan()",
                                             context + ": U");
        }
        // check correct dimension
        if (static_cast<idx>(U.rows()) != d_) {
            throw exception::MatrixMismatchSubsys("qpp::QCircuit::cCTRL_fan()",
                                                  context + ": U");
        }

        // check shift
        if (shift.has_value() && shift >= d_) {
            throw exception::OutOfRange("qpp::QCircuit::cCTRL_fan()",
                                        context + ": shift");
        }
        // END EXCEPTION CHECKS

        if (!name.has_value()) {
            auto gate_name =
                qpp::Gates::get_no_thread_local_instance().get_name(U);
            name = gate_name.has_value() ? "cCTRL_fan-" + gate_name.value()
                                         : "cCTRL_fan";
        }
        std::size_t hashU = hash_eigen(U);
        add_hash_(hashU, U);

        std::optional<std::vector<idx>> shift_vec;
        if (shift.has_value()) {
            shift_vec = std::vector<idx>{shift.value()};
        }
        circuit_.emplace_back(internal::QCircuitGateStep{
            internal::QCircuitGateStep::Type::cCTRL_FAN, hashU,
            std::vector<idx>{ctrl_dit}, target, shift_vec, name});

        gate_count_[hashU] += target.size();

        clean_dits_[ctrl_dit] = false;
        for (idx elem : target) {
            clean_qudits_[elem] = false;
        }

        return *this;
    }

    // multiple cctrl multiple targets
    /**
     * \brief Applies the single qudit controlled gate \a U with multiple
     * classical control dits listed in \a ctrl on every qudit listed in
     * \a target, i.e., cCTRL-cCTRL-...-cCTRL-U-U-...-U
     *
     * \param U Single qudit quantum gate
     * \param ctrl_dits Classical control dits indexes
     * \param target Target qudit indexes; the gate \a U is applied on every
     * one of them depending on the values of the classical control dits
     * \param shift Optional, performs the control as if the \a ctrl_dits
     * classical dits were incremented component-wise by \a shift. If present,
     * the size of \a shift must be the same as the size of \a ctrl_dits.
     * \param name Optional gate name
     * \return Reference to the current instance
     */
    QCircuit& cCTRL_fan(const cmat& U, const std::vector<idx>& ctrl_dits,
                        const std::vector<idx>& target,
                        std::optional<std::vector<idx>> shift = std::nullopt,
                        std::optional<std::string> name = std::nullopt) {
        // EXCEPTION CHECKS
        std::string context{"Step " + std::to_string(get_step_count())};

        // check valid ctrl_dits
        if (ctrl_dits.empty()) {
            throw exception::ZeroSize("qpp::QCircuit::cCTRL_fan()",
                                      context + ": ctrl_dits");
        }
        for (idx elem : ctrl_dits) {
            if (elem >= nc_) {
                throw exception::OutOfRange("qpp::QCircuit::cCTRL_fan()",
                                            context + ": ctrl_dits");
            }
        }
        // check no duplicates ctrl_dits
        if (!internal::check_no_duplicates(ctrl_dits)) {
            throw exception::Duplicates("qpp::QCircuit::cCTRL_fan()",
                                        context + ": ctrl_dits");
        }

        // check valid target
        if (target.empty()) {
            throw exception::ZeroSize("qpp::QCircuit::cCTRL_fan()",
                                      context + ": target");
        }
        for (idx elem : target) {
            if (elem >= nq_) {
                throw exception::OutOfRange("qpp::QCircuit::cCTRL_fan()",
                                            context + ": target");
            }
            // check target was not measured before
            if (was_measured_d(elem)) {
                throw exception::QuditAlreadyMeasured(
                    "qpp::QCircuit::cCTRL_fan()", context + ": target");
            }
        }
        // check no duplicates target
        if (!internal::check_no_duplicates(target)) {
            throw exception::Duplicates("qpp::QCircuit::cCTRL_fan()",
                                        context + ": target");
        }

        // check square matrix for the gate
        if (!internal::check_square_mat(U)) {
            throw exception::MatrixNotSquare("qpp::QCircuit::cCTRL_fan()",
                                             context + ": U");
        }
        // check correct dimension
        if (static_cast<idx>(U.rows()) != d_) {
            throw exception::MatrixMismatchSubsys("qpp::QCircuit::cCTRL_fan()",
                                                  context + ": U");
        }

        // check shift
        if (shift.has_value() && (shift.value().size() != ctrl_dits.size())) {
            throw exception::SizeMismatch("qpp::QCircuit::cCTRL_fan()",
                                          context + ": ctrl_dits/shift");
        }
        if (shift.has_value()) {
            for (idx elem : shift.value()) {
                if (elem >= d_) {
                    throw exception::OutOfRange("qpp::QCircuit::cCTRL_fan()",
                                                context + ": shift");
                }
            }
        }
        // END EXCEPTION CHECKS

        if (!name.has_value()) {
            auto gate_name =
                qpp::Gates::get_no_thread_local_instance().get_name(U);
            name = gate_name.has_value() ? "cCTRL_fan-" + gate_name.value()
                                         : "cCTRL_fan";
        }
        std::size_t hashU = hash_eigen(U);
        add_hash_(hashU, U);

        circuit_.emplace_back(internal::QCircuitGateStep{
            internal::QCircuitGateStep::Type::cCTRL_FAN, hashU, ctrl_dits,
            std::vector<idx>{target}, shift, name});

        gate_count_[hashU] += target.size();

        for (idx elem : ctrl_dits) {
            clean_dits_[elem] = false;
        }
        for (idx elem : target) {
            clean_qudits_[elem] = false;
        }

        return *this;
    }

    // multiple cctrl joint target
    /**
     * \brief Jointly applies the multiple-qudit controlled gate \a U with
     * multiple classical control dits listed in \a ctrl on the qudit
     * indexes specified by \a target, i.e., cCTRL-cCTRL-...-cCTRL-U_{joint}
     *
     * \param U Multiple-qudit quantum gate
     * \param ctrl_dits Classical control dits indexes
     * \param target Target qudit indexes where the gate \a U is applied
     * depending on the values of the classical control dits
     * \param shift Optional, performs the control as if the \a ctrl_dits
     * classical dits were incremented component-wise by \a shift. If present,
     * the size of \a shift must be the same as the size of \a ctrl_dits.
     * \param name Optional gate name
     * \return Reference to the current instance
     */
    QCircuit& cCTRL(const cmat& U, const std::vector<idx>& ctrl_dits,
                    const std::vector<idx>& target,
                    std::optional<std::vector<idx>> shift = std::nullopt,
                    std::optional<std::string> name = std::nullopt) {
        // EXCEPTION CHECKS
        std::string context{"Step " + std::to_string(get_step_count())};

        idx D_target = internal::safe_pow(d_, target.size());

        // check valid ctrl_dits
        if (ctrl_dits.empty()) {
            throw exception::ZeroSize("qpp::QCircuit::cCTRL()",
                                      context + ": ctrl_dits");
        }
        for (idx elem : ctrl_dits) {
            if (elem >= nc_) {
                throw exception::OutOfRange("qpp::QCircuit::cCTRL()",
                                            context + ": ctrl_dits");
            }
        }
        // check no duplicates ctrl_dits
        if (!internal::check_no_duplicates(ctrl_dits)) {
            throw exception::Duplicates("qpp::QCircuit::cCTRL()",
                                        context + ": ctrl_dits");
        }

        // check valid target
        if (target.empty()) {
            throw exception::ZeroSize("qpp::QCircuit::cCTRL()",
                                      context + ": target");
        }
        for (idx elem : target) {
            if (elem >= nq_) {
                throw exception::OutOfRange("qpp::QCircuit::cCTRL()",
                                            context + ": target");
            }
            // check target was not measured before
            if (was_measured_d(elem)) {
                throw exception::QuditAlreadyMeasured("qpp::QCircuit::cCTRL()",
                                                      context + ": target");
            }
        }
        // check no duplicates target
        if (!internal::check_no_duplicates(target)) {
            throw exception::Duplicates("qpp::QCircuit::cCTRL()",
                                        context + ": target");
        }

        // check square matrix for the gate
        if (!internal::check_square_mat(U)) {
            throw exception::MatrixNotSquare("qpp::QCircuit::cCTRL()",
                                             context + ": U");
        }
        // check correct dimension
        if (static_cast<idx>(U.rows()) != D_target) {
            throw exception::MatrixMismatchSubsys("qpp::QCircuit::cCTRL()",
                                                  context + ": U");
        }

        // check shift
        if (shift.has_value() && (shift.value().size() != ctrl_dits.size())) {
            throw exception::SizeMismatch("qpp::QCircuit::cCTRL()",
                                          context + ": ctrl_dits/shift");
        }
        if (shift.has_value()) {
            for (idx elem : shift.value()) {
                if (elem >= d_) {
                    throw exception::OutOfRange("qpp::QCircuit::cCTRL()",
                                                context + ": shift");
                }
            }
        }
        // END EXCEPTION CHECKS

        if (!name.has_value()) {
            auto gate_name =
                qpp::Gates::get_no_thread_local_instance().get_name(U);
            name =
                gate_name.has_value() ? "cCTRL-" + gate_name.value() : "cCTRL";
        }
        std::size_t hashU = hash_eigen(U);
        add_hash_(hashU, U);

        circuit_.emplace_back(
            internal::QCircuitGateStep{internal::QCircuitGateStep::Type::cCTRL,
                                       hashU, ctrl_dits, target, shift, name});

        ++gate_count_[hashU];

        for (idx elem : ctrl_dits) {
            clean_dits_[elem] = false;
        }
        for (idx elem : target) {
            clean_qudits_[elem] = false;
        }

        return *this;
    }

    // multiple cctrl single target
    /**
     * \brief Applies the single qudit controlled gate \a U with multiple
     * classical control dits listed in \a ctrl on the target qudit \a
     * target, i.e., cCTRL-cCTRL-...-cCTRL-U
     *
     * \param U Single qudit quantum gate
     * \param ctrl_dits Classical control dits indexes
     * \param target Target qudit index
     * \param shift Optional, performs the control as if the \a ctrl_dits
     * classical dits were incremented component-wise by \a shift. If present,
     * the size of \a shift must be the same as the size of \a ctrl_dits.
     * \param name Optional gate name
     * \return Reference to the current instance
     */
    QCircuit& cCTRL(const cmat& U, const std::vector<idx>& ctrl_dits,
                    idx target,
                    std::optional<std::vector<idx>> shift = std::nullopt,
                    std::optional<std::string> name = std::nullopt) {
        // EXCEPTION CHECKS
        std::string context{"Step " + std::to_string(get_step_count())};

        // check valid ctrl_dits
        if (ctrl_dits.empty()) {
            throw exception::ZeroSize("qpp::QCircuit::cCTRL()",
                                      context + ": ctrl_dits");
        }
        for (idx elem : ctrl_dits) {
            if (elem >= nc_) {
                throw exception::OutOfRange("qpp::QCircuit::cCTRL()",
                                            context + ": ctrl_dits");
            }
        }
        // check no duplicates ctrl_dits
        if (!internal::check_no_duplicates(ctrl_dits)) {
            throw exception::Duplicates("qpp::QCircuit::cCTRL()",
                                        context + ": ctrl_dits");
        }

        // check valid target
        if (target >= nq_) {
            throw exception::OutOfRange("qpp::QCircuit::cCTRL()",
                                        context + ": target");
        }
        // check target was not measured before
        if (was_measured_d(target)) {
            throw exception::QuditAlreadyMeasured("qpp::QCircuit::cCTRL()",
                                                  context + ": target");
        }

        // check square matrix for the gate
        if (!internal::check_square_mat(U)) {
            throw exception::MatrixNotSquare("qpp::QCircuit::cCTRL()",
                                             context + ": U");
        }
        // check correct dimension
        if (static_cast<idx>(U.rows()) != d_) {
            throw exception::MatrixMismatchSubsys("qpp::QCircuit::cCTRL()",
                                                  context + ": U");
        }

        // check shift
        if (shift.has_value() && (shift.value().size() != ctrl_dits.size())) {
            throw exception::SizeMismatch("qpp::QCircuit::cCTRL()",
                                          context + ": ctrl_dits/shift");
        }
        if (shift.has_value()) {
            for (idx elem : shift.value()) {
                if (elem >= d_) {
                    throw exception::OutOfRange("qpp::QCircuit::cCTRL()",
                                                context + ": shift");
                }
            }
        }
        // END EXCEPTION CHECKS

        if (!name.has_value()) {
            auto gate_name =
                qpp::Gates::get_no_thread_local_instance().get_name(U);
            name =
                gate_name.has_value() ? "cCTRL-" + gate_name.value() : "cCTRL";
        }
        std::size_t hashU = hash_eigen(U);
        add_hash_(hashU, U);

        circuit_.emplace_back(internal::QCircuitGateStep{
            internal::QCircuitGateStep::Type::cCTRL, hashU, ctrl_dits,
            std::vector<idx>{target}, shift, name});

        ++gate_count_[hashU];

        for (idx elem : ctrl_dits) {
            clean_dits_[elem] = false;
        }
        clean_qudits_[target] = false;

        return *this;
    }

    // single cctrl joint target
    /**
     * \brief Jointly applies the single qudit controlled gate \a U with
     * classical control dit \a ctrl on the qudit indexes specified by \a
     * target, i.e., cCTRL-U_{joint}
     *
     * \param U Multiple-qudit quantum gate
     * \param ctrl_dit Classical control dit index
     * \param target Target qudit indexes; the gate \a U is applied on every
     * one of them depending on the values of the classical control dits
     * \param shift Optional, performs the control as if the \a ctrl_dit
     * classical dit was incremented by \a shift
     * \param name Optional gate name
     * \return Reference to the current instance
     */
    QCircuit& cCTRL(const cmat& U, idx ctrl_dit, const std::vector<idx>& target,
                    std::optional<idx> shift = std::nullopt,
                    std::optional<std::string> name = std::nullopt) {
        // EXCEPTION CHECKS
        std::string context{"Step " + std::to_string(get_step_count())};

        idx D_target = internal::safe_pow(d_, target.size());

        // check valid ctrl_dit
        if (ctrl_dit >= nc_) {
            throw exception::OutOfRange("qpp::QCircuit::cCTRL()",
                                        context + ": ctrl_dit");
        }

        // check valid target
        if (target.empty()) {
            throw exception::ZeroSize("qpp::QCircuit::cCTRL()",
                                      context + ": target");
        }
        for (idx elem : target) {
            if (elem >= nq_) {
                throw exception::OutOfRange("qpp::QCircuit::cCTRL()",
                                            context + ": target");
            }
            // check target was not measured before
            if (was_measured_d(elem)) {
                throw exception::QuditAlreadyMeasured("qpp::QCircuit::cCTRL()",
                                                      context + ": target");
            }
        }
        // check no duplicates target
        if (!internal::check_no_duplicates(target)) {
            throw exception::Duplicates("qpp::QCircuit::cCTRL()",
                                        context + ": target");
        }

        // check square matrix for the gate
        if (!internal::check_square_mat(U)) {
            throw exception::MatrixNotSquare("qpp::QCircuit::cCTRL()",
                                             context + ": U");
        }
        // check correct dimension
        if (static_cast<idx>(U.rows()) != D_target) {
            throw exception::MatrixMismatchSubsys("qpp::QCircuit::cCTRL()",
                                                  context + ": U");
        }

        // check shift
        if (shift.has_value() && shift.value() >= d_) {
            throw exception::OutOfRange("qpp::QCircuit::cCTRL()",
                                        context + ": shift");
        }
        // END EXCEPTION CHECKS

        if (!name.has_value()) {
            auto gate_name =
                qpp::Gates::get_no_thread_local_instance().get_name(U);
            name =
                gate_name.has_value() ? "cCTRL-" + gate_name.value() : "cCTRL";
        }
        std::size_t hashU = hash_eigen(U);
        add_hash_(hashU, U);

        std::optional<std::vector<idx>> shift_vec;
        if (shift.has_value()) {
            shift_vec = std::vector<idx>{shift.value()};
        }
        circuit_.emplace_back(internal::QCircuitGateStep{
            internal::QCircuitGateStep::Type::cCTRL, hashU,
            std::vector<idx>{ctrl_dit}, target, shift_vec, name});

        ++gate_count_[hashU];

        clean_dits_[ctrl_dit] = false;
        for (idx elem : target) {
            clean_qudits_[elem] = false;
        }

        return *this;
    }

    // single cctrl single target
    /**
     * \brief Applies the single qubit controlled gate \a U with classical
     * control dit \a ctrl and target qudit \a target, i.e., cCTRL-U
     *
     * \param U Single qudit quantum gate
     * \param ctrl_dit Classical control dit index
     * \param target Target qudit index
     * \param shift Optional, performs the control as if the \a ctrl_dit
     * classical dit was incremented by \a shift
     * \param name Optional gate name
     * \return Reference to the current instance
     */
    QCircuit& cCTRL(const cmat& U, idx ctrl_dit, idx target,
                    std::optional<idx> shift = std::nullopt,
                    std::optional<std::string> name = std::nullopt) {
        // EXCEPTION CHECKS
        std::string context{"Step " + std::to_string(get_step_count())};

        // check valid ctrl_dit and target
        if (ctrl_dit >= nc_ || target >= nq_) {
            throw exception::OutOfRange("qpp::QCircuit::cCTRL()",
                                        context + ": ctrl/target");
        }
        if (was_measured_d(target)) {
            throw exception::QuditAlreadyMeasured("qpp::QCircuit::cCTRL()",
                                                  context + ": target");
        }

        // check square matrix for the gate
        if (!internal::check_square_mat(U)) {
            throw exception::MatrixNotSquare("qpp::QCircuit::cCTRL()",
                                             context + ": U");
        }
        // check correct dimension
        if (static_cast<idx>(U.rows()) != d_) {
            throw exception::MatrixMismatchSubsys("qpp::QCircuit::cCTRL()",
                                                  context + ": U");
        }

        // check shift
        if (shift.has_value() && shift.value() >= d_) {
            throw exception::OutOfRange("qpp::QCircuit::cCTRL()",
                                        context + ": shift");
        }
        // END EXCEPTION CHECKS

        if (!name.has_value()) {
            auto gate_name =
                qpp::Gates::get_no_thread_local_instance().get_name(U);
            name =
                gate_name.has_value() ? "cCTRL-" + gate_name.value() : "cCTRL";
        }

        std::size_t hashU = hash_eigen(U);
        add_hash_(hashU, U);

        std::optional<std::vector<idx>> shift_vec;
        if (shift.has_value()) {
            shift_vec = std::vector<idx>{shift.value()};
        }
        circuit_.emplace_back(internal::QCircuitGateStep{
            internal::QCircuitGateStep::Type::cCTRL, hashU,
            std::vector<idx>{ctrl_dit}, std::vector<idx>{target}, shift_vec,
            name});

        ++gate_count_[hashU];

        clean_dits_[ctrl_dit] = false;
        clean_qudits_[target] = false;

        return *this;
    }

    // Z measurement of single qudit
    /**
     * \brief Measures single qudit in the computational basis (Z-basis)
     *
     * \param target Target qudit index that is measured
     * \param c_reg Classical register where the value of the measurement is
     * being stored
     * \param destructive Destructive measurement, true by default
     * \param name Optional measurement name, default is "mZ"
     * \return Reference to the current instance
     */
    QCircuit& measure(idx target, idx c_reg, bool destructive = true,
                      std::optional<std::string> name = "mZ") {
        // EXCEPTION CHECKS
        std::string context{"Step " + std::to_string(get_step_count())};

        // measuring non-existing qudit
        if (target >= nq_) {
            throw exception::OutOfRange("qpp::QCircuit::measure()",
                                        context + ": target");
        }
        // trying to put the result into a non-existing classical slot
        if (c_reg >= nc_) {
            throw exception::OutOfRange("qpp::QCircuit::measure()",
                                        context + ": c_reg");
        }
        // qudit was measured before
        if (was_measured_d(target)) {
            throw exception::QuditAlreadyMeasured("qpp::QCircuit::measure()",
                                                  context + ": target");
        }
        // END EXCEPTION CHECKS

        std::size_t m_hash = hash_eigen(Gates::get_instance().Zd(d_));

        if (destructive) {
            measured_d_[target] = true;

            circuit_.emplace_back(internal::QCircuitMeasurementStep{
                internal::QCircuitMeasurementStep::Type::MEASURE,
                std::vector<std::size_t>{m_hash}, std::vector<idx>{target},
                c_reg, name});

        } else {
            measured_nd_[target] = true;

            circuit_.emplace_back(internal::QCircuitMeasurementStep{
                internal::QCircuitMeasurementStep::Type::MEASURE_ND,
                std::vector<std::size_t>{m_hash}, std::vector<idx>{target},
                c_reg, name});
        }
        ++measurement_count_[m_hash];

        clean_qudits_[target] = false;

        clean_dits_[c_reg] = false;
        measurement_dits_[c_reg] = true;

        return *this;
    }

    // Z measurement of multiple qudits
    /**
     * \brief Measures multiple qudits in the computational basis (Z-basis)
     *
     * \param target Target qudit indexes that are measured
     * \param c_reg Classical register where the values of the measurements
     * are being stored in order, that is, the measurement result
     * corresponding to target[0] is stored in \a c_reg, the one
     * corresponding to target[1] is stored in \a c_reg + 1 and so on
     * \param destructive Destructive measurement, true by default
     * \param name Optional measurement name, default is "mZ"
     * \return Reference to the current instance
     */
    QCircuit& measure(const std::vector<idx>& target, idx c_reg = 0,
                      bool destructive = true,
                      std::optional<std::string> name = "mZ") {
        // EXCEPTION CHECKS
        std::string context{"Step " + std::to_string(get_step_count())};

        // check valid target
        if (target.empty()) {
            throw exception::ZeroSize("qpp::QCircuit::measure()",
                                      context + ": target");
        }
        for (idx elem : target) {
            // measuring non-existing qudit
            if (elem >= nq_) {
                throw exception::OutOfRange("qpp::QCircuit::measure()",
                                            context + ": target");
            }
            // qudit was measured before
            if (was_measured_d(elem)) {
                throw exception::QuditAlreadyMeasured(
                    "qpp::QCircuit::measure()", context + ": target");
            }
        }
        // check no duplicates target
        if (!internal::check_no_duplicates(target)) {
            throw exception::Duplicates("qpp::QCircuit::measure()",
                                        context + ": target");
        }
        // not enough dits to store the result
        if (static_cast<idx>(target.size()) > nc_ ||
            c_reg > static_cast<idx>(nc_ - target.size())) {
            throw exception::OutOfRange("qpp::QCircuit::measure()",
                                        context + ": c_reg, target");
        }
        // END EXCEPTION CHECKS

        std::size_t m_hash = hash_eigen(Gates::get_instance().Zd(d_));

        if (destructive) {
            for (idx elem : target) {
                measured_d_[elem] = true;
            }
            circuit_.emplace_back(internal::QCircuitMeasurementStep{
                internal::QCircuitMeasurementStep::Type::MEASURE_MANY,
                std::vector<std::size_t>{m_hash}, target, c_reg, name});
        } else {
            for (idx elem : target) {
                measured_nd_[elem] = true;
            }
            circuit_.emplace_back(internal::QCircuitMeasurementStep{
                internal::QCircuitMeasurementStep::Type::MEASURE_MANY_ND,
                std::vector<std::size_t>{m_hash}, target, c_reg, name});
        }
        ++measurement_count_[m_hash];

        for (idx elem : target) {
            clean_qudits_[elem] = false;
        }

        for (idx i = 0; i < static_cast<idx>(target.size()); ++i) {
            clean_dits_[c_reg + i] = false;
            measurement_dits_[c_reg + i] = true;
        }

        return *this;
    }

    /**
     * \brief Measures all remaining measurable qudits in the computational
     * basis (Z-basis)
     *
     * \param c_reg Classical register where the values of the measurements
     * are being stored in order, that is, the measurement result
     * corresponding to the first measurable qudit is stored in \a c_reg,
     * the one corresponding to the second measurable qudit is stored in \a
     * c_reg + 1 and so on
     * \param destructive Destructive measurement, true by default
     * \param name Optional measurement name, default is "mZ"
     * \return Reference to the current instance
     */
    QCircuit& measure_all(idx c_reg = 0, bool destructive = true,
                          std::optional<std::string> name = "mZ") {
        // EXCEPTION CHECKS
        std::string context{"Step " + std::to_string(get_step_count())};

        std::vector<idx> non_measured = get_non_measured_d();
        if (non_measured.empty()) {
            throw exception::QuditAlreadyMeasured(
                "qpp::QCircuit::measure_all()",
                context + ": all qudits measured already");
        }
        if (static_cast<idx>(non_measured.size()) > nc_ ||
            c_reg > static_cast<idx>(nc_ - non_measured.size())) {
            throw exception::OutOfRange("qpp::QCircuit::measure_all()",
                                        context + ": c_reg");
        }
        // END EXCEPTION CHECKS

        return measure(non_measured, c_reg, destructive, std::move(name));
    }

    // measurement of single qudit in the orthonormal basis or rank-1
    // projectors specified by the columns of matrix V
    /**
     * \brief Measures single qudit in the orthonormal basis or rank-1
     * projectors specified by the columns of matrix \a V
     *
     * \param V Orthonormal basis or rank-1 projectors specified by the
     * columns of matrix \a V
     * \param target Target qudit index that is measured
     * \param c_reg Classical register where the value of the measurement is
     * stored
     * \param destructive Destructive measurement, true by default
     * \param name Optional measurement name
     * \return Reference to the current instance
     */
    QCircuit& measureV(const cmat& V, idx target, idx c_reg,
                       bool destructive = true,
                       std::optional<std::string> name = std::nullopt) {
        // EXCEPTION CHECKS
        std::string context{"Step " + std::to_string(get_step_count())};

        // measuring non-existing qudit
        if (target >= nq_) {
            throw exception::OutOfRange("qpp::QCircuit::measureV()",
                                        context + ": target");
        }
        // trying to put the result into a non-existing classical slot
        if (c_reg >= nc_) {
            throw exception::OutOfRange("qpp::QCircuit::measureV()",
                                        context + ": c_reg");
        }
        // qudit was measured before
        if (was_measured_d(target)) {
            throw exception::QuditAlreadyMeasured("qpp::QCircuit::measureV()",
                                                  context + ": target");
        }
        // check matrix V
        if (V.cols() == 0) {
            throw exception::ZeroSize("qpp::QCircuit::measureV()",
                                      context + ": V");
        }
        if (static_cast<idx>(V.rows()) != d_) {
            throw exception::DimsNotEqual("qpp::QCircuit::measureV()",
                                          context + ": V");
        }
        // END EXCEPTION CHECKS

        if (!name.has_value()) {
            auto gate_name =
                qpp::Gates::get_no_thread_local_instance().get_name(V);
            name = gate_name.has_value() ? "mV-" + gate_name.value() : "mV";
        }

        std::size_t hashV = hash_eigen(V);
        add_hash_(hashV, V);

        if (destructive) {
            measured_d_[target] = true;

            circuit_.emplace_back(internal::QCircuitMeasurementStep{
                internal::QCircuitMeasurementStep::Type::MEASURE_V,
                std::vector<std::size_t>{hashV}, std::vector<idx>{target},
                c_reg, name});
        } else {
            measured_nd_[target] = true;

            circuit_.emplace_back(internal::QCircuitMeasurementStep{
                internal::QCircuitMeasurementStep::Type::MEASURE_V_ND,
                std::vector<std::size_t>{hashV}, std::vector<idx>{target},
                c_reg, name});
        }
        ++measurement_count_[hashV];

        clean_qudits_[target] = false;

        clean_dits_[c_reg] = false;
        measurement_dits_[c_reg] = true;

        return *this;
    }

    // measurement of multiple qudits in the orthonormal basis or rank-1
    // projectors specified by the columns of matrix V
    /**
     * \brief Jointly measures multiple qudits in the orthonormal basis or
     * rank-1 projectors specified by the columns of matrix \a V
     *
     * \param V Orthonormal basis or rank-1 projectors specified by the
     * columns of matrix \a V
     * \param target Target qudit indexes that are jointly measured
     * \param c_reg Classical register where the value of the measurement is
     * stored
     * \param destructive Destructive measurement, true by default
     * \param name Optional measurement name
     * \return Reference to the current instance
     */
    QCircuit& measureV(const cmat& V, const std::vector<idx>& target, idx c_reg,
                       bool destructive = true,
                       std::optional<std::string> name = std::nullopt) {
        // EXCEPTION CHECKS
        std::string context{"Step " + std::to_string(get_step_count())};

        // check valid target
        if (target.empty()) {
            throw exception::ZeroSize("qpp::QCircuit::measureV()",
                                      context + ": target");
        }
        for (idx elem : target) {
            if (elem >= nq_) {
                throw exception::OutOfRange("qpp::QCircuit::measureV()",
                                            context + ": target");
            }
            // check target was not measured before
            if (was_measured_d(elem)) {
                throw exception::QuditAlreadyMeasured(
                    "qpp::QCircuit::measureV()", context + ": target");
            }
        }
        // check no duplicates target
        if (!internal::check_no_duplicates(target)) {
            throw exception::Duplicates("qpp::QCircuit::measureV()",
                                        context + ": target");
        }
        // trying to put the result into a non-existing classical slot
        if (c_reg >= nc_) {
            throw exception::OutOfRange("qpp::QCircuit::measureV()",
                                        context + ": c_reg");
        }
        // qudit was measured before
        for (idx elem : target) {
            if (was_measured_d(elem)) {
                throw exception::QuditAlreadyMeasured(
                    "qpp::QCircuit::measureV()", context + ": target");
            }
        }
        // check matrix V
        if (V.cols() == 0) {
            throw exception::ZeroSize("qpp::QCircuit::measureV()",
                                      context + ": V");
        }
        if (V.rows() !=
            static_cast<Eigen::Index>(std::pow(d_, target.size()))) {
            throw exception::DimsNotEqual("qpp::QCircuit::measureV()",
                                          context + ": V");
        }
        // END EXCEPTION CHECKS

        if (!name.has_value()) {
            auto gate_name =
                qpp::Gates::get_no_thread_local_instance().get_name(V);
            name = gate_name.has_value() ? "mV-" + gate_name.value() : "mV";
        }

        std::size_t hashV = hash_eigen(V);
        add_hash_(hashV, V);

        if (destructive) {
            for (idx elem : target) {
                measured_d_[elem] = true;
            }

            circuit_.emplace_back(internal::QCircuitMeasurementStep{
                internal::QCircuitMeasurementStep::Type::MEASURE_V_JOINT,
                std::vector<std::size_t>{hashV}, target, c_reg, name});
        } else {
            for (idx elem : target) {
                measured_nd_[elem] = true;
            }

            circuit_.emplace_back(internal::QCircuitMeasurementStep{
                internal::QCircuitMeasurementStep::Type::MEASURE_V_JOINT_ND,
                std::vector<std::size_t>{hashV}, target, c_reg, name});
        }
        ++measurement_count_[hashV];

        for (idx elem : target) {
            clean_qudits_[elem] = false;
        }

        clean_dits_[c_reg] = false;
        measurement_dits_[c_reg] = true;

        return *this;
    }

    // Z post-selection of single qudit
    /**
     * \brief Post-selects single qudit in the computational basis (Z-basis)
     *
     * \param target Target qudit index that is post-selected
     * \param ps_val Post-election value
     * \param c_reg Classical register where the value of the measurement is
     * being stored
     * \param destructive Destructive measurement, true by default
     * \param name Optional measurement name, default is "pZ"
     * \return Reference to the current instance
     */
    QCircuit& post_select(idx target, idx ps_val, idx c_reg,
                          bool destructive = true,
                          std::optional<std::string> name = "pZ") {
        // EXCEPTION CHECKS
        std::string context{"Step " + std::to_string(get_step_count())};

        // post-selecting non-existing qudit
        if (target >= nq_) {
            throw exception::OutOfRange("qpp::QCircuit::post_select()",
                                        context + ": target");
        }
        // post-selection value out of range
        if (ps_val >= d_) {
            throw exception::OutOfRange("qpp::QCircuit::post_select()",
                                        context + ": ps_val");
        }
        // trying to put the result into a non-existing classical slot
        if (c_reg >= nc_) {
            throw exception::OutOfRange("qpp::QCircuit::post_select()",
                                        context + ": c_reg");
        }
        // qudit was measured before
        if (was_measured_d(target)) {
            throw exception::QuditAlreadyMeasured(
                "qpp::QCircuit::post_select()", context + ": target");
        }
        // END EXCEPTION CHECKS

        std::size_t m_hash = hash_eigen(Gates::get_instance().Zd(d_));

        if (destructive) {
            measured_d_[target] = true;

            circuit_.emplace_back(internal::QCircuitMeasurementStep{
                internal::QCircuitMeasurementStep::Type::POST_SELECT,
                std::vector<std::size_t>{m_hash}, std::vector<idx>{target},
                c_reg, name, std::vector<idx>{ps_val}});

        } else {
            measured_nd_[target] = true;

            circuit_.emplace_back(internal::QCircuitMeasurementStep{
                internal::QCircuitMeasurementStep::Type::POST_SELECT_ND,
                std::vector<std::size_t>{m_hash}, std::vector<idx>{target},
                c_reg, name, std::vector<idx>{ps_val}});
        }
        ++measurement_count_[m_hash];

        clean_qudits_[target] = false;

        clean_dits_[c_reg] = false;
        measurement_dits_[c_reg] = true;

        return *this;
    }

    /**
     * \brief Post-selects multiple qudits in the computational basis (Z-basis)
     *
     * \param target Target qudit indexes that are post-selected
     * \param ps_vals Post-election values
     * \param c_reg Classical register where the values of the measurements
     * are being stored in order, that is, the measurement result
     * corresponding to target[0] is stored in \a c_reg, the one
     * corresponding to target[1] is stored in \a c_reg + 1 and so on
     * \param destructive Destructive measurement, true by default
     * \param name Optional post-selection name, default is "pZ"
     * \return Reference to the current instance
     */
    QCircuit& post_select(const std::vector<idx>& target,
                          const std::vector<idx>& ps_vals, idx c_reg,
                          bool destructive = true,
                          std::optional<std::string> name = "pZ") {
        // EXCEPTION CHECKS
        std::string context{"Step " + std::to_string(get_step_count())};

        // check valid target
        if (target.empty()) {
            throw exception::ZeroSize("qpp::QCircuit::post_select()",
                                      context + ": target");
        }
        for (idx elem : target) {
            // post-selecting non-existing qudit
            if (elem >= nq_) {
                throw exception::OutOfRange("qpp::QCircuit::post_select()",
                                            context + ": target");
            }
            // qudit was measured before
            if (was_measured_d(elem)) {
                throw exception::QuditAlreadyMeasured(
                    "qpp::QCircuit::post_select()", context + ": target");
            }
        }
        // check no duplicates target
        if (!internal::check_no_duplicates(target)) {
            throw exception::Duplicates("qpp::QCircuit::post_select()",
                                        context + ": target");
        }
        // not enough dits to store the result
        if (static_cast<idx>(target.size()) > nc_ ||
            c_reg > static_cast<idx>(nc_ - target.size())) {
            throw exception::OutOfRange("qpp::QCircuit::post_select()",
                                        context + ": c_reg, target");
        }
        // post-selection value out of range
        for (idx elem : ps_vals) {
            if (elem >= d_) {
                throw exception::OutOfRange("qpp::QCircuit::post_select()",
                                            context + ": ps_vals");
            }
        }
        // END EXCEPTION CHECKS

        std::size_t m_hash = hash_eigen(Gates::get_instance().Zd(d_));

        if (destructive) {
            for (idx elem : target) {
                measured_d_[elem] = true;
            }
            circuit_.emplace_back(internal::QCircuitMeasurementStep{
                internal::QCircuitMeasurementStep::Type::POST_SELECT_MANY,
                std::vector<std::size_t>{m_hash}, std::vector<idx>{target},
                c_reg, name, ps_vals});

        } else {
            for (idx elem : target) {
                measured_nd_[elem] = true;
            }
            circuit_.emplace_back(internal::QCircuitMeasurementStep{
                internal::QCircuitMeasurementStep::Type::POST_SELECT_MANY_ND,
                std::vector<std::size_t>{m_hash}, std::vector<idx>{target},
                c_reg, name, ps_vals});
        }
        ++measurement_count_[m_hash];

        for (idx elem : target) {
            clean_qudits_[elem] = false;
        }

        for (idx i = 0; i < static_cast<idx>(target.size()); ++i) {
            clean_dits_[c_reg + i] = false;
            measurement_dits_[c_reg + i] = true;
        }

        return *this;
    }

    // post-selection of single qudit in the orthonormal basis or rank-1
    // projectors specified by the columns of matrix V
    /**
     * \brief Post-selects single qudit in the orthonormal basis or rank-1
     * projectors specified by the columns of matrix \a V
     *
     * \param V Orthonormal basis or rank-1 projectors specified by the
     * columns of matrix \a V
     * \param target Target qudit index that is post-selected
     * \param ps_val Post-election value
     * \param c_reg Classical register where the value of the measurement is
     * stored
     * \param destructive Destructive measurement, true by default
     * \param name Optional measurement name
     * \return Reference to the current instance
     */
    QCircuit& post_selectV(const cmat& V, idx target, idx ps_val, idx c_reg,
                           bool destructive = true,
                           std::optional<std::string> name = std::nullopt) {
        // EXCEPTION CHECKS
        std::string context{"Step " + std::to_string(get_step_count())};

        // measuring non-existing qudit
        if (target >= nq_) {
            throw exception::OutOfRange("qpp::QCircuit::post_selectV()",
                                        context + ": target");
        }
        // post-selection value out of range
        if (ps_val >= d_) {
            throw exception::OutOfRange("qpp::QCircuit::post_select()",
                                        context + ": ps_val");
        }
        // trying to put the result into a non-existing classical slot
        if (c_reg >= nc_) {
            throw exception::OutOfRange("qpp::QCircuit::post_selectV()",
                                        context + ": c_reg");
        }
        // qudit was measured before
        if (was_measured_d(target)) {
            throw exception::QuditAlreadyMeasured(
                "qpp::QCircuit::post_selectV()", context + ": target");
        }
        // check matrix V
        if (V.cols() == 0) {
            throw exception::ZeroSize("qpp::QCircuit::post_selectV()",
                                      context + ": V");
        }
        if (static_cast<idx>(V.rows()) != d_) {
            throw exception::DimsNotEqual("qpp::QCircuit::post_selectV()",
                                          context + ": V");
        }
        // END EXCEPTION CHECKS

        if (!name.has_value()) {
            auto gate_name =
                qpp::Gates::get_no_thread_local_instance().get_name(V);
            name = gate_name.has_value() ? "pV-" + gate_name.value() : "pV";
        }

        std::size_t hashV = hash_eigen(V);
        add_hash_(hashV, V);

        if (destructive) {
            measured_d_[target] = true;

            circuit_.emplace_back(internal::QCircuitMeasurementStep{
                internal::QCircuitMeasurementStep::Type::POST_SELECT_V,
                std::vector<std::size_t>{hashV}, std::vector<idx>{target},
                c_reg, name, std::vector<idx>{ps_val}});
        } else {
            measured_nd_[target] = true;

            circuit_.emplace_back(internal::QCircuitMeasurementStep{
                internal::QCircuitMeasurementStep::Type::POST_SELECT_V_ND,
                std::vector<std::size_t>{hashV}, std::vector<idx>{target},
                c_reg, name, std::vector<idx>{ps_val}});
        }
        ++measurement_count_[hashV];

        clean_qudits_[target] = false;

        clean_dits_[c_reg] = false;
        measurement_dits_[c_reg] = true;

        return *this;
    }

    // post-selection of multiple qudits in the orthonormal basis or rank-1
    // projectors specified by the columns of matrix V
    /**
     * \brief Jointly post-select multiple qudits in the orthonormal basis or
     * rank-1 projectors specified by the columns of matrix \a V
     *
     * \param V Orthonormal basis or rank-1 projectors specified by the
     * columns of matrix \a V
     * \param target Target qudit indexes that are jointly post-selected
     * \param ps_val Post-election value
     * \param c_reg Classical register where the value of the measurement is
     * stored
     * \param destructive Destructive measurement, true by default
     * \param name Optional post-selection name
     * \return Reference to the current instance
     */
    QCircuit& post_selectV(const cmat& V, const std::vector<idx>& target,
                           idx ps_val, idx c_reg, bool destructive = true,
                           std::optional<std::string> name = std::nullopt) {
        // EXCEPTION CHECKS
        std::string context{"Step " + std::to_string(get_step_count())};

        // check valid target
        if (target.empty()) {
            throw exception::ZeroSize("qpp::QCircuit::post_selectV()",
                                      context + ": target");
        }
        for (idx elem : target) {
            if (elem >= nq_) {
                throw exception::OutOfRange("qpp::QCircuit::post_selectV()",
                                            context + ": target");
            }
            // check target was not measured before
            if (was_measured_d(elem)) {
                throw exception::QuditAlreadyMeasured(
                    "qpp::QCircuit::post_selectV()", context + ": target");
            }
        }
        // check no duplicates target
        if (!internal::check_no_duplicates(target)) {
            throw exception::Duplicates("qpp::QCircuit::post_selectV()",
                                        context + ": target");
        }
        // trying to put the result into a non-existing classical slot
        if (c_reg >= nc_) {
            throw exception::OutOfRange("qpp::QCircuit::post_selectV()",
                                        context + ": c_reg");
        }
        // qudit was measured before
        for (idx elem : target) {
            if (was_measured_d(elem)) {
                throw exception::QuditAlreadyMeasured(
                    "qpp::QCircuit::post_selectV()", context + ": target");
            }
        }
        // post-selection value out of range
        if (ps_val >= static_cast<idx>(V.cols())) {
            throw exception::OutOfRange("qpp::QCircuit::post_selectV()",
                                        context + ": ps_val");
        }
        // trying to put the result into a non-existing classical slot
        if (c_reg >= nc_) {
            throw exception::OutOfRange("qpp::QCircuit::post_selectV()",
                                        context + ": c_reg");
        }
        // check matrix V
        if (V.cols() == 0) {
            throw exception::ZeroSize("qpp::QCircuit::post_selectV()",
                                      context + ": V");
        }
        if (V.rows() !=
            static_cast<Eigen::Index>(std::pow(d_, target.size()))) {
            throw exception::DimsNotEqual("qpp::QCircuit::post_selectV()",
                                          context + ": V");
        }
        // END EXCEPTION CHECKS

        if (!name.has_value()) {
            auto gate_name =
                qpp::Gates::get_no_thread_local_instance().get_name(V);
            name = gate_name.has_value() ? "pV-" + gate_name.value() : "pV";
        }

        std::size_t hashV = hash_eigen(V);
        add_hash_(hashV, V);

        if (destructive) {
            for (idx elem : target) {
                measured_d_[elem] = true;
            }
            circuit_.emplace_back(internal::QCircuitMeasurementStep{
                internal::QCircuitMeasurementStep::Type::POST_SELECT_V_JOINT,
                std::vector<std::size_t>{hashV}, std::vector<idx>{target},
                c_reg, name, std::vector<idx>{ps_val}});
        } else {
            for (idx elem : target) {
                measured_nd_[elem] = true;
            }
            circuit_.emplace_back(internal::QCircuitMeasurementStep{
                internal::QCircuitMeasurementStep::Type::POST_SELECT_V_JOINT_ND,
                std::vector<std::size_t>{hashV}, std::vector<idx>{target},
                c_reg, name, std::vector<idx>{ps_val}});
        }
        ++measurement_count_[hashV];

        for (idx elem : target) {
            clean_qudits_[elem] = false;
        }

        clean_dits_[c_reg] = false;
        measurement_dits_[c_reg] = true;

        return *this;
    }

    /**
     * \brief Discards single qudit by measuring it destructively in the
     * computational basis (Z-basis) and discarding the measurement result
     *
     * \param target Target qudit index that is discarded
     * \param name Optional discard operation name, default is "discard"
     * \return Reference to the current instance
     */
    QCircuit& discard(idx target, std::optional<std::string> name = "discard") {
        // EXCEPTION CHECKS
        std::string context{"Step " + std::to_string(get_step_count())};

        // discarding non-existing qudit
        if (target >= nq_) {
            throw exception::OutOfRange("qpp::QCircuit::discard()",
                                        context + ": target");
        }
        // qudit was measured before
        if (was_measured_d(target)) {
            throw exception::QuditAlreadyMeasured("qpp::QCircuit::discard()",
                                                  context + ": target");
        }
        // END EXCEPTION CHECKS

        std::size_t m_hash = hash_eigen(cmat::Zero(d_, d_));

        measured_d_[target] = true;

        circuit_.emplace_back(internal::QCircuitMeasurementStep{
            internal::QCircuitMeasurementStep::Type::DISCARD,
            std::vector<std::size_t>{m_hash}, std::vector<idx>{target},
            static_cast<idx>(-1), name});
        //++measurement_count_[m_hash];

        clean_qudits_[target] = false;

        return *this;
    }

    /**
     * \brief Discards multiple qudits by measuring them destructively in
     * the computational basis (Z-basis) and discarding the measurement
     * result
     *
     * \param target Target qudit indexes that are discarded
     * \param name Optional discard operation name, default is "discard"
     * \return Reference to the current instance
     */
    QCircuit& discard(const std::vector<idx>& target,
                      std::optional<std::string> name = "discard") {
        // EXCEPTION CHECKS
        std::string context{"Step " + std::to_string(get_step_count())};

        // check valid target
        if (target.empty()) {
            throw exception::ZeroSize("qpp::QCircuit::discard()",
                                      context + ": target");
        }
        for (idx elem : target) {
            // discarding non-existing qudit
            if (elem >= nq_) {
                throw exception::OutOfRange("qpp::QCircuit::discard()",
                                            context + ": target");
            }
            // qudit was measured before
            if (was_measured_d(elem)) {
                throw exception::QuditAlreadyMeasured(
                    "qpp::QCircuit::discard()", context + ": target");
            }
        }
        // END EXCEPTION CHECKS

        std::size_t m_hash = hash_eigen(cmat::Zero(d_, d_));

        for (idx elem : target) {
            measured_d_[elem] = true;
        }

        circuit_.emplace_back(internal::QCircuitMeasurementStep{
            internal::QCircuitMeasurementStep::Type::DISCARD_MANY,
            std::vector<std::size_t>{m_hash}, target, static_cast<idx>(-1),
            name});
        //++measurement_count_[m_hash];

        for (idx elem : target) {
            clean_qudits_[elem] = false;
        }

        return *this;
    }

    /**
     * \brief No operation (no-op)
     *
     * \note If the underlying step is executed on a noisy engine, then the
     * noise acts before it
     *
     * \return Reference to the current instance
     */
    QCircuit& nop() {
        circuit_.emplace_back(internal::QCircuitNOPStep{});

        return *this;
    }

    // reset single qudit
    /**
     * \brief Resets single qudit by first measuring it non-destructively in
     * the computational basis and discarding the measurement result,
     * followed by shifting it back to the \f$|0\rangle\f$ state
     *
     * \param target Target qudit index that is reset
     * \param name Optional reset operation name, default is "reset"
     * \return Reference to the current instance
     */
    QCircuit& reset(idx target, std::optional<std::string> name = "reset") {
        // EXCEPTION CHECKS
        std::string context{"Step " + std::to_string(get_step_count())};

        // resetting non-existing qudit
        if (target >= nq_) {
            throw exception::OutOfRange("qpp::QCircuit::reset()",
                                        context + ": target");
        }
        // qudit was measured before
        if (was_measured_d(target)) {
            throw exception::QuditAlreadyMeasured("qpp::QCircuit::reset()",
                                                  context + ": target");
        }
        // END EXCEPTION CHECKS

        std::size_t m_hash = hash_eigen(mprj({0}, d_));

        circuit_.emplace_back(internal::QCircuitMeasurementStep{
            internal::QCircuitMeasurementStep::Type::RESET,
            std::vector<std::size_t>{m_hash}, std::vector<idx>{target},
            static_cast<idx>(-1), name});
        //++measurement_count_[m_hash];

        return *this;
    }

    // reset multiple qudits
    /**
     * \brief Resets multiple qudits by first measuring them
     * non-destructively in the computational basis and discarding the
     * measurement results, followed by shifting them back to the
     * \f$|0\cdots 0\rangle\f$ state
     *
     * \param target Target qudit indexes that are reset
     * \param name Optional measurement name, default is "reset"
     * \return Reference to the current instance
     */
    QCircuit& reset(const std::vector<idx>& target,
                    std::optional<std::string> name = "reset") {
        // EXCEPTION CHECKS
        std::string context{"Step " + std::to_string(get_step_count())};

        // check valid target
        if (target.empty()) {
            throw exception::ZeroSize("qpp::QCircuit::reset()",
                                      context + ": target");
        }
        for (idx elem : target) {
            // resetting non-existing qudit
            if (elem >= nq_) {
                throw exception::OutOfRange("qpp::QCircuit::reset()",
                                            context + ": target");
            }
            // qudit was measured before
            if (was_measured_d(elem)) {
                throw exception::QuditAlreadyMeasured("qpp::QCircuit::reset()",
                                                      context + ": target");
            }
        }
        // END EXCEPTION CHECKS

        std::size_t m_hash = hash_eigen(mprj({0}, d_));

        circuit_.emplace_back(internal::QCircuitMeasurementStep{
            internal::QCircuitMeasurementStep::Type::RESET_MANY,
            std::vector<std::size_t>{m_hash}, target, static_cast<idx>(-1),
            name});
        //++measurement_count_[m_hash];

        return *this;
    }

    /**
     * \brief Replicates the circuit, in place
     * \note The circuit should not contain any measurement steps
     *
     * \param n Number of repetitions. If \a n == 1, returns the original
     * circuit.
     * \return Reference to the current instance
     */
    QCircuit& replicate(idx n) {
        // EXCEPTION CHECKS
        if (n == 0) {
            throw exception::OutOfRange("qpp::QCircuit::replicate()", "n");
        }
        if (this->removes_qudits()) {
            throw exception::QuditAlreadyMeasured("qpp::QCircuit::replicate()",
                                                  "n");
        }
        // END EXCEPTION CHECKS

        std::vector<VarStep> result;
        result.reserve(n * circuit_.size());
        for (idx i = 0; i < n; ++i) {
            result.insert(result.end(), circuit_.begin(), circuit_.end());
        }
        circuit_ = std::move(result);

        // update gate counts
        for (auto& elem : gate_count_) {
            elem.second *= n;
        }

        // update measurement counts
        for (auto& elem : measurement_count_) {
            elem.second *= n;
        }

        return *this;
    }

    /**
     * \brief Composes (appends) a quantum circuit description to the end of
     * the current one
     * \see qpp::QCircuit::couple_circuit_left() and
     * qpp::QCircuit::couple_circuit_right()
     *
     * \note If the qudit indexes of the added quantum circuit description
     * do not totally overlap with the indexes of the current quantum
     * circuit description, then the required number of additional qudits
     * are automatically added to the current quantum circuit description
     *
     * \param other Quantum circuit description
     * \param pos_qudit The index of the first/top qudit of \a other quantum
     * circuit description relative to the index of the first/top qudit of
     * the current quantum circuit description, with the rest following in
     * order. If negative or greater than the total number of qudits of the
     * current quantum circuit description, then the required number of
     * additional qudits are automatically added to the current quantum
     * circuit description.
     * \param pos_dit Optional, the first classical dit of \a other quantum
     * circuit description is inserted before the \a pos_dit classical dit
     * index of the current quantum circuit description (in the classical
     * dits array), the rest following in order. If absent (default),
     * insertion is performed at the end.
     * \return Reference to the current instance
     */
    QCircuit& compose_circuit(QCircuit other, bigint pos_qudit,
                              std::optional<idx> pos_dit = std::nullopt) {
        // EXCEPTION CHECKS
        // check equal dimensions
        if (other.d_ != d_) {
            throw exception::DimsNotEqual("qpp::QCircuit::compose_circuit()",
                                          "other");
        }
        // check classical dits
        if (!pos_dit.has_value()) {
            pos_dit = nc_;
        } else {
            if (internal::is_negative(pos_dit.value()) ||
                pos_dit.value() > nc_) {
                throw exception::OutOfRange("qpp::QCircuit::compose_circuit()",
                                            "pos_dit");
            }
        }
        // check that overlapping qudits (in the current instance) were not
        // already destructively measured
        if (pos_qudit < 0 &&
            (pos_qudit + static_cast<bigint>(other.nq_)) >= 0) {
            for (idx i = 0;
                 i < std::min(static_cast<idx>(pos_qudit +
                                               static_cast<bigint>(other.nq_)),
                              nq_);
                 ++i) {
                if (was_measured_d(i)) {
                    throw exception::QuditAlreadyMeasured(
                        "qpp::QCircuit::compose_circuit()",
                        "Current qpp::QCircuit instance");
                }
            }
        }
        if (pos_qudit >= 0 && static_cast<idx>(pos_qudit) < nq_) {
            for (idx i = 0;
                 i < std::min(static_cast<idx>(nq_ - pos_qudit), other.nq_);
                 ++i) {
                if (was_measured_d(pos_qudit + i)) {
                    throw exception::QuditAlreadyMeasured(
                        "qpp::QCircuit::compose_circuit()",
                        "Current qpp::QCircuit instance");
                }
            }
        }
        // END EXCEPTION CHECKS

        // STEP 0: add additional qudits (if needed) and classical dits from
        // the to-be-coupled circuit
        if (pos_qudit < 0) {
            // add qudits before beginning
            idx extra_qudits = std::abs(pos_qudit);
            add_qudit(extra_qudits, 0);
            pos_qudit = 0;
        }
        if (pos_qudit >= 0) {
            // add qudits after beginning
            idx tmp = pos_qudit + other.nq_;
            if (tmp > nq_) {
                idx extra_qudits = tmp - nq_;
                add_qudit(extra_qudits);
            }
        }
        add_dit(other.nc_, pos_dit.value());

        // STEP 1: update [c]ctrl and target indexes of other
        auto gate_step_visitor = [&](internal::QCircuitGateStep& gate_step) {
            // update the cctrl indexes
            if (is_cCTRL(gate_step)) {
                for (idx& pos : gate_step.ctrl_.value()) {
                    pos += pos_dit.value();
                }
            }
            // update the ctrl indexes
            if (is_CTRL(gate_step) && pos_qudit >= 0) {
                for (idx& pos : gate_step.ctrl_.value()) {
                    pos += pos_qudit;
                }
            }

            // update the target indexes
            if (pos_qudit >= 0) {
                for (idx& pos : gate_step.target_) {
                    pos += pos_qudit;
                }
            }
        };
        auto measurement_step_visitor =
            [&](internal::QCircuitMeasurementStep& measurement_step) {
                measurement_step.c_reg_ += pos_dit.value();
                if (pos_qudit >= 0) {
                    for (idx& pos : measurement_step.target_) {
                        pos += pos_qudit;
                    }
                }
            };
        auto nop_step_visitor = [&](internal::QCircuitNOPStep&) {};
        auto conditional_step_visitor =
            [&](const internal::QCircuitConditionalStep&) {};

        visit_circuit_(other.circuit_, conditional_step_visitor,
                       gate_step_visitor, measurement_step_visitor,
                       nop_step_visitor);

        // STEP 2
        // replace the corresponding elements of measured_, measured_nd_,
        // and clean_qudits_ with the ones of other
        if (pos_qudit < 0) {
            std::copy_if(other.measured_d_.begin(), other.measured_d_.end(),
                         measured_d_.begin(), [](bool val) { return val; });
            std::copy_if(other.measured_nd_.begin(), other.measured_nd_.end(),
                         measured_nd_.begin(), [](bool val) { return val; });
            std::copy_if(other.clean_qudits_.begin(), other.clean_qudits_.end(),
                         clean_qudits_.begin(), [](bool val) { return !val; });
        } else {
            std::copy_if(other.measured_d_.begin(), other.measured_d_.end(),
                         std::next(measured_d_.begin(), pos_qudit),
                         [](bool val) { return val; });
            std::copy_if(other.measured_nd_.begin(), other.measured_nd_.end(),
                         std::next(measured_nd_.begin(), pos_qudit),
                         [](bool val) { return val; });
            std::copy_if(other.clean_qudits_.begin(), other.clean_qudits_.end(),
                         std::next(clean_qudits_.begin(), pos_qudit),
                         [](bool val) { return !val; });
        }

        // STEP 3
        // replace the corresponding elements of clean_dits_ and
        // measurement_dits_ with the ones of other
        std::copy(other.clean_dits_.begin(), other.clean_dits_.end(),
                  std::next(clean_dits_.begin(),
                            static_cast<std::ptrdiff_t>(pos_dit.value())));
        std::copy(other.measurement_dits_.begin(),
                  other.measurement_dits_.end(),
                  std::next(measurement_dits_.begin(),
                            static_cast<std::ptrdiff_t>(pos_dit.value())));

        // STEP 4: append the copy of other to the current instance
        circuit_.insert(circuit_.end(), other.circuit_.begin(),
                        other.circuit_.end());

        // STEP 5: modify gate counts, hash tables etc. accordingly
        // update matrix hash table
        for (auto& elem : other.cmat_hash_tbl_) {
            cmat_hash_tbl_[elem.first] = elem.second;
        }
        // update gate counts
        for (auto& elem : other.gate_count_) {
            gate_count_[elem.first] += elem.second;
        }
        // update measurement counts
        for (auto& elem : other.measurement_count_) {
            measurement_count_[elem.first] += elem.second;
        }

        return *this;
    }

    /**
     * \brief Couples in-place a quantum circuit description to the current
     * quantum circuit description, with the to-be-added quantum circuit
     * description placed at the left (beginning) of the current quantum
     * circuit description
     * \see qpp::QCircuit::couple_circuit_right() and
     * qpp::QCircuit::compose_circuit()
     *
     * \note The added quantum circuit description should not contain any
     * destructive  measurements and should not be larger than the current
     * quantum circuit description, i.e., all qudit indexes of the added
     * quantum circuit description must match with qudits from the current
     * quantum circuit description (and \a other should not contain any
     * destructive measurements)
     *
     * \note The classical dits are not relabeled
     *
     * \param other Quantum circuit description
     * \param target Qudit indexes of the current circuit description where
     * the qudits of \a other are being coupled, i.e., the first/top qudit
     * of \a other quantum circuit description is coupled with the target[0]
     * qudit of the current circuit description, and so on
     * \param pos_dit Optional, the first classical dit of \a other quantum
     * circuit description is inserted before the \a pos_dit classical dit
     * index of the current quantum circuit description (in the classical
     * dits array), the rest following in order. If absent (default),
     * insertion is performed at the end.
     * \return Reference to the current instance
     */
    QCircuit& couple_circuit_left(QCircuit other,
                                  const std::vector<idx>& target,
                                  std::optional<idx> pos_dit = std::nullopt) {
        // EXCEPTION CHECKS
        // check equal dimensions
        if (other.d_ != d_) {
            throw exception::DimsNotEqual(
                "qpp::QCircuit::couple_circuit_left()", "other");
        }
        // check classical dits
        if (!pos_dit.has_value()) {
            pos_dit = nc_;
        } else if (internal::is_negative(pos_dit.value()) ||
                   pos_dit.value() > nc_) {
            throw exception::OutOfRange("qpp::QCircuit::couple_circuit_left()",
                                        "pos_dit");
        }
        // check no measurement for the coupled circuit
        if (!other.get_measured_d().empty()) {
            throw exception::QuditAlreadyMeasured(
                "qpp::QCircuit::couple_circuit_left()", "other");
        }
        // check valid target
        if (static_cast<idx>(target.size()) != other.nq_) {
            throw exception::OutOfRange("qpp::QCircuit::couple_circuit_left()",
                                        "target");
        }
        if (static_cast<idx>(target.size()) > nq_) {
            throw exception::OutOfRange("qpp::QCircuit::couple_circuit_left()",
                                        "target");
        }
        if (!internal::check_no_duplicates(target)) {
            throw exception::Duplicates("qpp::QCircuit::couple_circuit_left()",
                                        "target");
        }
        for (idx elem : target) {
            if (elem >= nq_) {
                throw exception::OutOfRange(
                    "qpp::QCircuit::couple_circuit_left()", "target");
            }
        }
        // check matching qudits (in the current instance) were not already
        // measured destructively
        for (idx elem : target) {
            if (was_measured_d(elem)) {
                throw exception::QuditAlreadyMeasured(
                    "qpp::QCircuit::couple_circuit_left()", "target");
            }
        }
        // END EXCEPTION CHECKS

        // STEP 0: insert classical dits from the to-be-coupled circuit
        add_dit(other.nc_, pos_dit.value());

        // STEP 1: update [c]ctrl and target indexes of other
        // update gate_step indexes of other
        auto gate_step_visitor = [&](internal::QCircuitGateStep& gate_step) {
            // update the cctrl indexes
            if (is_cCTRL(gate_step)) {
                for (idx& dit : gate_step.ctrl_.value()) {
                    dit += pos_dit.value();
                }
            }
            // update the ctrl indexes
            if (is_CTRL(gate_step)) {
                for (idx& pos : gate_step.ctrl_.value()) {
                    pos = target[pos];
                }
            }
            // update the target indexes
            for (idx& pos : gate_step.target_) {
                pos = target[pos];
            }
        };
        // update measurement_step indexes of other
        auto measurement_step_visitor =
            [&](internal::QCircuitMeasurementStep& measurement_step) {
                measurement_step.c_reg_ += pos_dit.value();
                for (idx& pos : measurement_step.target_) {
                    pos = target[pos];
                }
            };
        auto nop_step_visitor = [&](internal::QCircuitNOPStep&) {};
        auto conditional_step_visitor =
            [&](const internal::QCircuitConditionalStep&) {};

        visit_circuit_(other.circuit_, conditional_step_visitor,
                       gate_step_visitor, measurement_step_visitor,
                       nop_step_visitor);

        // TODO check this

        // replace the corresponding elements of measured_, measured_nd_,
        // clean_qudits_, clean_dits_, and measurement_dits_ with the ones
        // of other
        for (idx i = 0; i < static_cast<idx>(other.measured_d_.size()); ++i) {
            if (other.measured_d_[i]) {
                measured_d_[target[i]] = true;
            }
        }
        for (idx i = 0; i < static_cast<idx>(other.measured_nd_.size()); ++i) {
            if (other.measured_nd_[i]) {
                measured_nd_[target[i]] = true;
            }
        }
        for (idx i = 0; i < static_cast<idx>(other.clean_qudits_.size()); ++i) {
            if (!other.clean_qudits_[i]) {
                clean_qudits_[target[i]] = false;
            }
        }

        for (idx i = 0; i < static_cast<idx>(other.clean_dits_.size()); ++i) {
            if (!other.clean_dits_[i]) {
                clean_dits_[target[i]] = false;
            }
        }
        for (idx i = 0; i < static_cast<idx>(other.measurement_dits_.size());
             ++i) {
            if (other.measurement_dits_[i]) {
                measurement_dits_[target[i]] = true;
            }
        }

        // STEP 2: append the copy of other to the current instance
        circuit_.insert(circuit_.begin(), other.circuit_.begin(),
                        other.circuit_.end());

        // STEP 3: modify gate counts, hash tables etc. accordingly
        // update matrix hash table
        for (auto& elem : other.cmat_hash_tbl_) {
            cmat_hash_tbl_[elem.first] = elem.second;
        }
        // update gate counts
        for (auto& elem : other.gate_count_) {
            gate_count_[elem.first] += elem.second;
        }
        // update measurement counts
        for (auto& elem : other.measurement_count_) {
            measurement_count_[elem.first] += elem.second;
        }

        return *this;
    }

    /**
     * \brief Couples in-place a quantum circuit description to the current
     * quantum circuit description, with the to-be-added quantum circuit
     * description placed at the right (end) of the current quantum circuit
     * description
     * \see qpp::QCircuit::couple_circuit_left() and
     * qpp::QCircuit::compose_circuit()
     *
     * \note The added quantum circuit description cannot be larger than
     * the current quantum circuit description, i.e., all qudit indexes of
     * the added quantum circuit description must match with qudits from the
     * current quantum circuit description (and the latter should not
     * contain any destructive measurements on the coupled qudits)
     *
     * \note The classical dits are not relabeled
     *
     * \param other Quantum circuit description
     * \param target Qudit indexes of the current circuit description where
     * the qudits of \a other are being coupled, i.e., the first/top qudit
     * of \a other quantum circuit description is coupled with the target[0]
     * qudit of the current circuit description, and so on
     * \param pos_dit Optional, the first classical dit of \a other quantum
     * circuit description is inserted before the \a pos_dit classical dit
     * index of the current quantum circuit description (in the classical
     * dits array), the rest following in order. If absent (default),
     * insertion is performed at the end.
     * \return Reference to the current instance
     */
    QCircuit& couple_circuit_right(QCircuit other,
                                   const std::vector<idx>& target,
                                   std::optional<idx> pos_dit = std::nullopt) {
        // EXCEPTION CHECKS
        // check equal dimensions
        if (other.d_ != d_) {
            throw exception::DimsNotEqual(
                "qpp::QCircuit::couple_circuit_right()", "other");
        }
        // check classical dits
        if (!pos_dit.has_value()) {
            pos_dit = nc_;
        } else if (internal::is_negative(pos_dit.value()) ||
                   pos_dit.value() > nc_) {
            throw exception::OutOfRange("qpp::QCircuit::couple_circuit_right()",
                                        "pos_dit");
        }
        // check valid target
        if (static_cast<idx>(target.size()) != other.nq_) {
            throw exception::OutOfRange("qpp::QCircuit::couple_circuit_right()",
                                        "target");
        }
        if (static_cast<idx>(target.size()) > nq_) {
            throw exception::OutOfRange("qpp::QCircuit::couple_circuit_right()",
                                        "target");
        }
        if (!internal::check_no_duplicates(target)) {
            throw exception::Duplicates("qpp::QCircuit::couple_circuit_right()",
                                        "target");
        }
        for (idx elem : target) {
            if (elem >= nq_) {
                throw exception::OutOfRange(
                    "qpp::QCircuit::couple_circuit_right()", "target");
            }
        }
        // check matching qudits (in the current instance) were not already
        // measured destructively
        for (idx elem : target) {
            if (was_measured_d(elem)) {
                throw exception::QuditAlreadyMeasured(
                    "qpp::QCircuit::couple_circuit_right()", "target");
            }
        }
        // END EXCEPTION CHECKS

        // STEP 0: insert classical dits from the to-be-coupled circuit
        add_dit(other.nc_, pos_dit.value());

        // STEP 1
        //
        // update [c]ctrl, target, and measurement indexes of other
        // update gate indexes
        auto gate_step_visitor = [&](internal::QCircuitGateStep& gate_step) {
            // update the cctrl indexes
            if (is_cCTRL(gate_step)) {
                for (idx& dit : gate_step.ctrl_.value()) {
                    dit += pos_dit.value();
                }
            }
            // update the ctrl indexes
            if (is_CTRL(gate_step)) {
                for (idx& pos : gate_step.ctrl_.value()) {
                    pos = target[pos];
                }
            }
            // update the target indexes
            for (idx& pos : gate_step.target_) {
                pos = target[pos];
            }
        };
        // update measurement indexes
        auto measurement_step_visitor =
            [&](internal::QCircuitMeasurementStep& measurement_step) {
                measurement_step.c_reg_ += pos_dit.value();
                for (idx& pos : measurement_step.target_) {
                    pos = target[pos];
                }
            };
        auto nop_step_visitor = [&](internal::QCircuitNOPStep&) {};
        auto conditional_step_visitor =
            [&](const internal::QCircuitConditionalStep&) {};

        visit_circuit_(other.circuit_, conditional_step_visitor,
                       gate_step_visitor, measurement_step_visitor,
                       nop_step_visitor);

        // TODO check this

        // replace the corresponding elements of measured_, measured_nd_,
        // clean_qudits_, clean_dits_, and measurement_dits_ with the ones
        // of other
        for (idx i = 0; i < static_cast<idx>(other.measured_d_.size()); ++i) {
            if (other.measured_d_[i]) {
                measured_d_[target[i]] = true;
            }
        }
        for (idx i = 0; i < static_cast<idx>(other.measured_nd_.size()); ++i) {
            if (other.measured_nd_[i]) {
                measured_nd_[target[i]] = true;
            }
        }
        for (idx i = 0; i < static_cast<idx>(other.clean_qudits_.size()); ++i) {
            if (!other.clean_qudits_[i]) {
                clean_qudits_[target[i]] = false;
            }
        }

        for (idx i = 0; i < static_cast<idx>(other.clean_dits_.size()); ++i) {
            if (!other.clean_dits_[i]) {
                clean_dits_[target[i]] = false;
            }
        }
        for (idx i = 0; i < static_cast<idx>(other.measurement_dits_.size());
             ++i) {
            if (other.measurement_dits_[i]) {
                measurement_dits_[target[i]] = true;
            }
        }

        // STEP 2: append the copy of other to the current instance
        circuit_.insert(circuit_.end(), other.circuit_.begin(),
                        other.circuit_.end());

        // STEP 3: modify gate counts, hash tables etc. accordingly
        // update matrix hash table
        for (auto& elem : other.cmat_hash_tbl_) {
            cmat_hash_tbl_[elem.first] = elem.second;
        }
        // update gate counts
        for (auto& elem : other.gate_count_) {
            gate_count_[elem.first] += elem.second;
        }
        // update measurement counts
        for (auto& elem : other.measurement_count_) {
            measurement_count_[elem.first] += elem.second;
        }

        return *this;
    }

    /**
     * \brief Composes (appends) a controlled quantum circuit description to
     * the end of the current one, with the current instance acting as the
     * control
     * \see qpp::QCircuit::compose_circuit()
     *
     * \note For each gate or controlled-gate in the target circuit
     * description, a control is being added from the control circuit
     * description control qudits. Classical control gates and
     * measurements/resets/discards in the target circuit description are
     * left unchanged.
     *
     * \note If the qudit indexes of the added quantum circuit description
     * do not totally overlap with the indexes of the current quantum
     * circuit description, then the required number of additional qudits
     * are automatically added to the current quantum circuit description
     *
     * \param ctrl Control qubits
     * \param qc_target Quantum circuit description
     * \param pos_qudit The index of the first/top qudit of \a qc_target
     * quantum circuit description relative to the index of the first/top
     * qudit of the current quantum circuit description, with the rest
     * following in order. If negative or greater than the total number of
     * qudits of the current quantum circuit description, then the required
     * number of additional qudits are automatically added to the current
     * quantum circuit description.
     * \param shift Optional, performs the control as if the \a ctrl qudit state
     * was \f$X\f$-incremented by \a shift
     * \param pos_dit Optional, the first classical dit of \a
     * qc_target quantum circuit description is inserted before the \a
     * pos_dit classical dit index of the current quantum circuit
     * description (in the classical dits array), the rest following in
     * order. If absent (default), insertion is performed at the end.
     * \return Reference to the current instance
     */
    QCircuit&
    compose_CTRL_circuit(const std::vector<idx>& ctrl, QCircuit qc_target,
                         bigint pos_qudit,
                         std::optional<std::vector<idx>> shift = std::nullopt,
                         std::optional<idx> pos_dit = std::nullopt) {
        // EXCEPTION CHECKS
        // check non-empty control circuit
        if (this->get_nq() == 0) {
            throw exception::ZeroSize("qpp::QCircuit::compose_CTRL_circuit()",
                                      "Current qpp::QCircuit instance");
        }

        // check non-empty qc_target
        if (qc_target.get_nq() == 0) {
            throw exception::ZeroSize("qpp::QCircuit::compose_CTRL_circuit()",
                                      "qc_target");
        }

        // check equal dimensions
        if (qc_target.d_ != d_) {
            throw exception::DimsNotEqual(
                "qpp::QCircuit::compose_CTRL_circuit()", "qc_target");
        }

        // check classical dits
        if (!pos_dit.has_value()) {
            pos_dit = nc_;
        } else {
            if (internal::is_negative(pos_dit.value()) ||
                pos_dit.value() > nc_) {
                throw exception::OutOfRange(
                    "qpp::QCircuit::compose_CTRL_circuit()", "pos_dit");
            }
        }

        // check that overlapping qudits (in the current instance) were not
        // already destructively measured
        if (pos_qudit < 0 &&
            (pos_qudit + static_cast<bigint>(qc_target.nq_)) >= 0) {
            for (idx i = 0;
                 i < std::min(static_cast<idx>(pos_qudit + static_cast<bigint>(
                                                               qc_target.nq_)),
                              nq_);
                 ++i) {
                if (was_measured_d(i)) {
                    throw exception::QuditAlreadyMeasured(
                        "qpp::QCircuit::compose_CTRL_circuit()",
                        "Current qpp::QCircuit instance");
                }
            }
        }
        if (pos_qudit >= 0 && static_cast<idx>(pos_qudit) < nq_) {
            for (idx i = 0;
                 i < std::min(static_cast<idx>(nq_ - pos_qudit), qc_target.nq_);
                 ++i) {
                if (was_measured_d(pos_qudit + i)) {
                    throw exception::QuditAlreadyMeasured(
                        "qpp::QCircuit::compose_CTRL_circuit()",
                        "Current qpp::QCircuit instance");
                }
            }
        }

        // check valid ctrl
        if (ctrl.empty()) {
            throw exception::ZeroSize("qpp::QCircuit::compose_CTRL_circuit()",
                                      "ctrl");
        }
        for (idx elem : ctrl) {
            if (elem >= nq_) {
                throw exception::OutOfRange(
                    "qpp::QCircuit::compose_CTRL_circuit()", "ctrl");
            }
            // check ctrl was not measured before
            if (was_measured_d(elem)) {
                throw exception::QuditAlreadyMeasured(
                    "qpp::QCircuit::compose_CTRL_circuit()", "ctrl");
            }
        }
        // check no duplicates ctrl
        if (!internal::check_no_duplicates(ctrl)) {
            throw exception::Duplicates("qpp::QCircuit::compose_CTRL_circuit()",
                                        "ctrl");
        }
        // check that ctrl and target do not overlap
        idx target_begin = pos_qudit;                    // including it
        idx target_end = pos_qudit + qc_target.get_nq(); // excluding it
        for (idx elem : ctrl) {
            if (elem >= target_begin && elem < target_end) {
                throw exception::OutOfRange(
                    "qpp::QCircuit::compose_CTRL_circuit()", "ctrl/qc_target");
            }
        }

        // check shift
        if (shift.has_value() && (shift.value().size() != ctrl.size())) {
            throw exception::SizeMismatch(
                "qpp::QCircuit::compose_CTRL_circuit()", "ctrl/shift");
        }
        if (shift.has_value()) {
            for (idx elem : shift.value()) {
                if (elem >= d_) {
                    throw exception::OutOfRange(
                        "qpp::QCircuit::compose_CTRL_circuit()", "shift");
                }
            }
        }
        // END EXCEPTION CHECKS

        // compose the current instance with qc_target
        idx end_ctrl_circuit = this->get_step_count();
        this->compose_circuit(qc_target, pos_qudit, pos_dit);
        idx end_composed_circuit = this->get_step_count();

        // modify all gate steps to control gate steps in the added
        // qc_target then enlarge ctrl and shift accordingly
        auto gate_step_visitor = [&](internal::QCircuitGateStep& gate_step) {
            switch (gate_step.gate_type_) {
                // don't do anything to these gates and return
                case internal::QCircuitGateStep::Type::NONE:
                case internal::QCircuitGateStep::Type::cCTRL:
                case internal::QCircuitGateStep::Type::cCTRL_FAN:
                    return;
                // change the type of these gates to CTRL
                case internal::QCircuitGateStep::Type::SINGLE:
                case internal::QCircuitGateStep::Type::TWO:
                case internal::QCircuitGateStep::Type::THREE:
                case internal::QCircuitGateStep::Type::JOINT:
                    gate_step.gate_type_ =
                        internal::QCircuitGateStep::Type::CTRL;
                    break;
                // change the type of FAN gate to CTRL_FAN
                case internal::QCircuitGateStep::Type::FAN:
                    gate_step.gate_type_ =
                        internal::QCircuitGateStep::Type::CTRL_FAN;
                    break;
                    // don't do anything to these gates
                case internal::QCircuitGateStep::Type::CTRL:
                case internal::QCircuitGateStep::Type::CTRL_FAN:
                    break;
            }

            // ctrl gate
            if (gate_step.ctrl_.has_value()) {
                // enlarge ctrl
                gate_step.ctrl_.value().insert(gate_step.ctrl_.value().begin(),
                                               ctrl.begin(), ctrl.end());
                // enlarge shifts
                // both shifts present
                if (shift.has_value() && gate_step.shift_.has_value()) {
                    gate_step.shift_.value().insert(
                        gate_step.shift_.value().begin(), shift.value().begin(),
                        shift.value().end());
                }
                // current shift present (ctrl circuit)
                else if (shift.has_value() && !gate_step.shift_.has_value()) {
                    // fill with zero at the end
                    gate_step.shift_ =
                        std::vector<idx>(gate_step.ctrl_.value().size(), 0);
                    // copy shift at the beginning
                    std::copy(shift.value().begin(), shift.value().end(),
                              gate_step.shift_.value().begin());
                }
                // old shift present (qc_target)
                else if (!shift.has_value() && gate_step.shift_.has_value()) {
                    gate_step.shift_.value().insert(
                        gate_step.shift_.value().begin(),
                        ctrl.size() - gate_step.shift_.value().size(), 0);
                }
            }
            // regular gate
            else {
                gate_step.ctrl_ = ctrl;
                if (shift.has_value()) {
                    gate_step.shift_ = shift;
                }
            }
        };
        // measurements are left unchanged
        auto measurement_step_visitor =
            [&](const internal::QCircuitMeasurementStep&) {};
        // NOPs are left unchanged
        auto nop_step_visitor = [&](const internal::QCircuitNOPStep&) {};
        auto conditional_step_visitor =
            [&](const internal::QCircuitConditionalStep&) {};
        for (idx i = end_ctrl_circuit; i < end_composed_circuit; ++i) {
            visit_circuit_step_(circuit_[i], conditional_step_visitor,
                                gate_step_visitor, measurement_step_visitor,
                                nop_step_visitor);
        }

        return *this;
    }

    /**
     * \brief Kronecker product with another quantum circuit description, in
     * place
     *
     * \param qc Quantum circuit description
     * \return Reference to the current instance
     */
    QCircuit& kron(QCircuit qc) {
        compose_circuit(std::move(qc), static_cast<bigint>(nq_));

        return *this;
    }

    /**
     * \brief Returns true if the quantum circuit description contains any
     * measurements, false otherwise
     *
     * \return True if the quantum circuit description contains any
     * measurements, false otherwise
     */
    bool has_measurements() const noexcept {
        return std::find_if(circuit_.begin(), circuit_.end(), [](auto&& arg) {
                   return std::holds_alternative<
                       internal::QCircuitMeasurementStep>(arg);
               }) != circuit_.end();
    }

    /**
     * \brief Returns true if the quantum circuit description
     * contains any operations/measurements that remove qudits,
     * false otherwise
     *
     * \return True if the quantum circuit description contains any
     * operations/measurements that remove qudits, false otherwise
     */
    bool removes_qudits() const noexcept {
        return std::find_if(circuit_.begin(), circuit_.end(), [](auto&& arg) {
                   if (std::holds_alternative<
                           internal::QCircuitMeasurementStep>(arg)) {
                       switch (std::get<internal::QCircuitMeasurementStep>(arg)
                                   .measurement_type_) {
                           case internal::QCircuitMeasurementStep::Type::
                               MEASURE:
                           case internal::QCircuitMeasurementStep::Type::
                               MEASURE_MANY:
                           case internal::QCircuitMeasurementStep::Type::
                               MEASURE_V:
                           case internal::QCircuitMeasurementStep::Type::
                               MEASURE_V_JOINT:
                           case internal::QCircuitMeasurementStep::Type::
                               DISCARD:
                           case internal::QCircuitMeasurementStep::Type::
                               DISCARD_MANY:
                           case internal::QCircuitMeasurementStep::Type::
                               POST_SELECT:
                           case internal::QCircuitMeasurementStep::Type::
                               POST_SELECT_MANY:
                           case internal::QCircuitMeasurementStep::Type::
                               POST_SELECT_V:
                           case internal::QCircuitMeasurementStep::Type::
                               POST_SELECT_V_JOINT:
                               return true;
                           default:
                               return false;
                       }
                   }
                   return false;
               }) != circuit_.end();
    }

    /**
     * \brief Adjoint quantum circuit description, in place
     * \note The circuit should not contain any measurement steps
     *
     * \return Reference to the current instance
     */
    QCircuit& adjoint() {
        // EXCEPTION CHECKS
        if (this->has_measurements()) {
            throw exception::QuditAlreadyMeasured(
                "qpp::QCircuit::adjoint()", "Current qpp::QCircuit instance");
        }
        // END EXCEPTION CHECKS

        auto htbl = cmat_hash_tbl_; // copy the gate hash table of other
        cmat_hash_tbl_.clear();

        std::reverse(circuit_.begin(), circuit_.end());

        auto gate_step_visitor = [&](internal::QCircuitGateStep& gate_step) {
            // get the gate and its corresponding
            // hash
            std::size_t hashU = gate_step.gate_hash_;
            cmat U = htbl[hashU];

            // compute the adjoints
            cmat Udagger = qpp::adjoint(U);
            std::size_t hashUdagger = hash_eigen(Udagger);

            // modify and add hash
            gate_step.gate_hash_ = hashUdagger;
            if (gate_step.name_.has_value()) {
                gate_step.name_.value() += "+";
            }
            add_hash_(hashUdagger, Udagger);
        };
        auto measurement_step_visitor =
            [&](const internal::QCircuitMeasurementStep&) {};
        auto nop_step_visitor = [&](const internal::QCircuitNOPStep&) {};
        auto conditional_step_visitor =
            [&](const internal::QCircuitConditionalStep&) {};

        visit_circuit_(circuit_, conditional_step_visitor, gate_step_visitor,
                       measurement_step_visitor, nop_step_visitor);

        return *this;
    }

    /**
     * \brief Checks whether a qudit in the circuit was used before
     * or not
     * \see qpp::QCircuit::get_clean_qudits(),
     * qpp::QCircuit::get_dirty_qudits()
     *
     * \param i Qudit index
     * \return True if the qudit \a i was used before (by a gate
     * and/or measurement, either destructive or non-destructive),
     * false otherwise
     */
    bool is_clean_qudit(idx i) const {
        // EXCEPTION CHECKS
        // check valid target
        if (i >= nq_) {
            throw exception::OutOfRange("qpp::QCircuit::is_clean_qudit()", "i");
        }
        // END EXCEPTION CHECKS

        return clean_qudits_[i];
    }

    /**
     * \brief Checks whether a classical dit in the circuit was used before
     * or not
     * \see qpp::QCircuit::get_clean_dits(), qpp::QCircuit::get_dirty_dits()
     *
     * \param i Classical dit index
     * \return True if the classical dit \a i was used before (by a cCTRL
     * gate and/or measurement, either destructive or non-destructive),
     * false otherwise
     */
    bool is_clean_dit(idx i) const {
        // EXCEPTION CHECKS
        // check valid target
        if (i >= nc_) {
            throw exception::OutOfRange("qpp::QCircuit::is_clean_dit()", "i");
        }
        // END EXCEPTION CHECKS

        return clean_dits_[i];
    }

    /**
     * \brief Checks whether a classical dit in the circuit was used to
     * store the result of a measurement (either destructive or
     * non-destructive)
     * \see qpp::QCircuit::get_measurement_dits()
     *
     * \param i Classical dit index
     * \return True if the classical dit \a i was used before to
     * store the result of a measurement, false otherwise
     */
    bool is_measurement_dit(idx i) const {
        // EXCEPTION CHECKS
        // check valid target
        if (i >= nc_) {
            throw exception::OutOfRange("qpp::QCircuit::is_measurement_dit()",
                                        "i");
        }
        // END EXCEPTION CHECKS

        return measurement_dits_[i];
    }

    /**
     * \brief Vector of clean qudits
     * \see qpp::QCircuit::is_clean_qudit(),
     * qpp::QCircuit::get_dirty_qudits()
     *
     * \return Vector of clean qudits
     */
    std::vector<idx> get_clean_qudits() const {
        std::vector<idx> result;
        for (idx i = 0; i < nq_; ++i) {
            if (is_clean_qudit(i)) {
                result.emplace_back(i);
            }
        }

        return result;
    }

    /**
     * \brief Vector of dirty qudits
     * \see qpp::QCircuit::is_clean_qudit(),
     * qpp::QCircuit::get_clean_qudits()
     *
     * \return Vector of dirty qudits
     */
    std::vector<idx> get_dirty_qudits() const {
        return complement(get_clean_qudits(), get_nq());
    }

    /**
     * \brief Vector of clean classical dits
     * \see qpp::QCircuit::is_clean_dit(),
     * qpp::QCircuit::get_dirty_dits()
     *
     * \return Vector of clean classical dits
     */
    std::vector<idx> get_clean_dits() const {
        std::vector<idx> result;
        for (idx i = 0; i < nc_; ++i) {
            if (is_clean_dit(i)) {
                result.emplace_back(i);
            }
        }

        return result;
    }

    /**
     * \brief Vector of dirty classical dits
     * \see qpp::QCircuit::is_clean_dit(),
     * qpp::QCircuit::get_clean_dits()
     *
     * \return Vector of dirty classical dits
     */
    std::vector<idx> get_dirty_dits() const {
        return complement(get_clean_dits(), get_nc());
    }

    /**
     * \brief Removes clean qudit from the quantum circuit
     * description and relabels the rest of the qudits accordingly
     * \see qpp::QCircuit::is_clean_qudit(),
     * qpp::QCircuit::compress()
     *
     * \param target Target clean qudit index that is removed
     * \return Reference to the current instance
     */
    QCircuit& remove_clean_qudit(idx target) {
        // EXCEPTION CHECKS
        // check valid target and clean qudit
        if (target >= nq_ || !is_clean_qudit(target)) {
            throw exception::OutOfRange("qpp::QCircuit::remove_clean_qudit()",
                                        "target");
        }
        // END EXCEPTION CHECKS

        auto gate_step_visitor = [&](internal::QCircuitGateStep& gate_step) {
            if (is_CTRL(gate_step)) {
                for (idx& pos : gate_step.ctrl_.value()) {
                    if (pos > target) {
                        --pos;
                    }
                }
            }
            for (idx& pos : gate_step.target_) {
                if (pos > target) {
                    --pos;
                }
            }
        };
        auto measurement_step_visitor =
            [&](internal::QCircuitMeasurementStep& measurement_step) {
                for (idx& pos : measurement_step.target_) {
                    if (pos > target) {
                        --pos;
                    }
                }
            };
        auto nop_step_visitor = [&](const internal::QCircuitNOPStep&) {};
        auto conditional_step_visitor =
            [&](const internal::QCircuitConditionalStep&) {};

        visit_circuit_(circuit_, conditional_step_visitor, gate_step_visitor,
                       measurement_step_visitor, nop_step_visitor);

        clean_qudits_.erase(std::next(clean_qudits_.begin(),
                                      static_cast<std::ptrdiff_t>(target)));
        --nq_;

        return *this;
    }

    /**
     * \brief Removes clean classical dit from the quantum circuit
     * description and relabels the rest of the classical dits accordingly
     * \see qpp::QCircuit::is_clean_dit(), qpp::QCircuit::compress()
     *
     * \param target Target clean classical dit index that is
     * removed
     * \return Reference to the current instance
     */
    QCircuit& remove_clean_dit(idx target) {
        // EXCEPTION CHECKS
        // check valid target and clean dit
        if (target >= nc_ || !is_clean_dit(target)) {
            throw exception::OutOfRange("qpp::QCircuit::remove_clean_dit()",
                                        "target");
        }
        // END EXCEPTION CHECKS

        auto gate_step_visitor = [&](internal::QCircuitGateStep& gate_step) {
            if (is_cCTRL(gate_step)) {
                for (idx& pos : gate_step.ctrl_.value()) {
                    if (pos > target) {
                        --pos;
                    }
                }
            }
        };
        auto measurement_step_visitor =
            [&](internal::QCircuitMeasurementStep& measurement_step) {
                if (measurement_step.c_reg_ > target) {
                    --measurement_step.c_reg_;
                }
            };
        auto nop_step_visitor = [&](const internal::QCircuitNOPStep&) {};
        auto conditional_step_visitor =
            [&](const internal::QCircuitConditionalStep&) {};

        visit_circuit_(circuit_, conditional_step_visitor, gate_step_visitor,
                       measurement_step_visitor, nop_step_visitor);

        clean_dits_.erase(std::next(clean_dits_.begin(),
                                    static_cast<std::ptrdiff_t>(target)));
        --nc_;

        return *this;
    }

    /**
     * \brief Removes list of clean qudits from the quantum circuit
     * description and relabels the rest of the qudits accordingly
     * \see qpp::QCircuit::is_clean_qudit(), qpp::QCircuit::compress()
     *
     * \param target Target clean qudit indexes that are removed
     * \return Reference to the current instance
     */
    QCircuit& remove_clean_qudits(std::vector<idx> target) {
        // EXCEPTION CHECKS
        // check valid target
        for (idx elem : target) {
            // removing non-existing or non-clean qudit
            if (elem >= nq_ || !is_clean_qudit(elem)) {
                throw exception::OutOfRange(
                    "qpp::QCircuit::remove_clean_qudits()", "target");
            }
        }
        // END EXCEPTION CHECKS

        // sort the target
        std::sort(target.begin(), target.end());
        idx dirty = 0;

        for (idx elem : target) {
            remove_clean_qudit(elem - dirty++);
        }

        return *this;
    }

    /**
     * \brief Removes list of clean classical dits from the quantum circuit
     * description and relabels the rest of the classical dits accordingly
     *
     * \see qpp::QCircuit::is_clean_dit(), qpp::QCircuit::compress()
     *
     * \param target Target clean classical dit indexes that are
     * removed
     * \return Reference to the current instance
     */
    QCircuit& remove_clean_dits(std::vector<idx> target) {
        // EXCEPTION CHECKS
        // check valid target
        for (idx elem : target) {
            // removing non-existing or non-clean dit
            if (elem >= nc_ || !is_clean_dit(elem)) {
                throw exception::OutOfRange(
                    "qpp::QCircuit::remove_clean_dits()", "target");
            }
        }
        // END EXCEPTION CHECKS

        // sort the target
        std::sort(target.begin(), target.end());
        idx dirty = 0;

        for (idx elem : target) {
            remove_clean_dit(elem - dirty++);
        }

        return *this;
    }

    /**
     * \brief Removes all clean qudits form the quantum circuit
     * description and relabels the rest of the qudits accordingly
     * \see qpp::QCircuit::remove_clean_qudits(),
     * qpp::QCircuit::remove_clean_dits()
     *
     * \param compress_dits If true, removes clean classical dits, false by
     * default.
     * \return Reference to the current instance
     */
    QCircuit& compress(bool compress_dits = false) {
        remove_clean_qudits(get_clean_qudits());
        if (compress_dits) {
            remove_clean_dits(get_clean_dits());
        }

        return *this;
    }

    /**
     * \brief Equality operator
     * \note Ignores names (e.g., circuit names, gate names etc.)
     * and does not perform any circuit simplifications, in other
     * words the circuits have to have the exact same number of
     * qubits/classical dits and the exact same gates/measurements
     * placed in the exact same order. For example, the circuit
     * \f$X_1 Z_2\f$ is considered different from \f$Z_2 X_1\f$,
     * although logically they are the same.
     *
     * \param rhs Quantum circuit description against which the
     * equality is being tested
     * \return True if the quantum circuit descriptions are equal, false
     * otherwise
     */
    bool operator==(const QCircuit& rhs) const noexcept {
        return std::tie(rhs.circuit_, rhs.nq_, rhs.nc_, rhs.measured_d_,
                        rhs.measured_nd_, rhs.d_, rhs.cmat_hash_tbl_) ==
               std::tie(circuit_, nq_, nc_, measured_d_, measured_nd_, d_,
                        cmat_hash_tbl_);
    }

    /**
     * \brief Inequality operator
     *
     * \param rhs Quantum circuit description against which the
     * inequality is being tested
     * \return True if the quantum circuit descriptions are not equal, false
     * otherwise
     */
    bool operator!=(const QCircuit& rhs) const noexcept {
        return !(*this == rhs);
    }

    /**
     * \brief qpp::IJSON::to_JSON() override
     *
     * Displays the quantum circuit description in JSON format
     *
     * \param enclosed_in_curly_brackets If true, encloses the
     * result in curly brackets
     * \return String containing the JSON representation of the quantum
     * circuit description
     */
    std::string to_JSON(bool enclosed_in_curly_brackets = true) const override;

  private:
    /**
     * \brief qpp::IDisplay::display() override
     *
     * Writes to the output stream a textual representation of the
     * quantum circuit
     *
     * \param os Output stream passed by reference
     * \return Reference to the output stream
     */
    std::ostream& display(std::ostream& os) const override;
}; /* class QCircuit */

template <class... Fns>
void QCircuit::visit_circuit_step_(QCircuit::VarStep& circuit_step,
                                   Fns&&... fns) {
    std::visit(overloaded{std::forward<Fns>(fns)...}, circuit_step);
}

template <class... Fns>
void QCircuit::visit_circuit_step_(const QCircuit::VarStep& circuit_step,
                                   Fns&&... fns) const {
    std::visit(overloaded{std::forward<Fns>(fns)...}, circuit_step);
}

template <class... Fns>
void QCircuit::visit_circuit_(std::vector<VarStep>& circuit, Fns&&... fns) {
    for (auto& step : circuit) {
        std::visit(overloaded{std::forward<Fns>(fns)...}, step);
    }
}

template <class... Fns>
void QCircuit::visit_circuit_(const std::vector<VarStep>& circuit,
                              Fns&&... fns) const {
    for (auto& step : circuit) {
        std::visit(overloaded{std::forward<Fns>(fns)...}, step);
    }
}

namespace internal {
/**
 * \class qpp::internal::QCircuitIterator
 * \brief Quantum circuit description bound-checking (safe) forward iterator
 *
 * \note The iterator is a const_iterator by default
 */
class QCircuitIterator {
    /**
     * \class qpp::internal::QCircuitIterator::value_type_
     * \brief Value type class for qpp::internal::QCircuitIterator
     */
    class value_type_ : public IDisplay /*, IJSON*/ {
        ///< non-owning pointer to the grand-parent const quantum circuit
        ///< description
        const QCircuit* qc_ptr_;
        idx ip_;                 ///< instruction pointer
        QCircuit::VarStep step_; ///< current circuit step
      public:
        /**
         * \brief Constructor
         *
         * \param qc_ptr Pointer to constant quantum circuit description
         */
        explicit value_type_(const QCircuit* qc_ptr, idx ip,
                             QCircuit::VarStep step)
            : qc_ptr_{qc_ptr}, ip_{ip}, step_{std::move(step)} {}

        // silence -Weffc++ class has pointer data members
        /**
         * \brief Default copy constructor
         */
        value_type_(const value_type_&) = default;

        // silence -Weffc++ class has pointer data members
        /**
         * \brief Default copy assignment operator
         *
         * \return Reference to the current instance
         */
        value_type_& operator=(const value_type_&) = default;

        /**
         * \brief Equality operator
         *
         * \param rhs Instance against which the equality is being tested
         * \return True if the instances not are equal (bitwise), false
         * otherwise
         */
        bool operator==(const value_type_& rhs) {
            return std::tie(qc_ptr_, ip_, step_) ==
                   std::tie(rhs.qc_ptr_, rhs.ip_, rhs.step_);
        }

        /**
         * \brief Inequality operator
         *
         * \param rhs Instance against which the inequality is being tested
         * \return True if the instances are not equal (bit by bit), false
         * otherwise
         */
        bool operator!=(const value_type_& rhs) { return !(*this == rhs); }

        // getters
        /**
         * \brief Pointer to underlying quantum circuit description
         * \return Pointer to underlying quantum circuit description
         */
        const QCircuit* get_qc_ptr() const { return qc_ptr_; }

        /**
         * \brief Current quantum circuit description instruction
         * pointer
         *
         * \return Current quantum circuit description instruction
         * pointer
         */
        idx get_ip() const { return ip_; }

        /**
         * \brief Current quantum circuit description step
         *
         * \return Current quantum circuit description step
         */
        QCircuit::VarStep get_step() const { return step_; }

      private:
        std::ostream& display(std::ostream& os) const override {
            // field spacing for the step number
            idx text_width = std::to_string(qc_ptr_->get_step_count()).size();

            os << std::left;
            os << std::setw(static_cast<int>(text_width)) << ip_ << ": ";
            os << std::right;

            auto gate_step_visitor =
                [&](const internal::QCircuitGateStep& gate_step) {
                    os << gate_step;
                };
            auto measurement_step_visitor = [&](const internal::
                                                    QCircuitMeasurementStep&
                                                        measurement_step) {
                switch (measurement_step.measurement_type_) {
                    case internal::QCircuitMeasurementStep::Type::NONE:
                        break;
                    case internal::QCircuitMeasurementStep::Type::MEASURE:
                    case internal::QCircuitMeasurementStep::Type::MEASURE_MANY:
                    case internal::QCircuitMeasurementStep::Type::MEASURE_V:
                    case internal::QCircuitMeasurementStep::Type::
                        MEASURE_V_JOINT:
                        os << "|> ";
                        break;
                    case internal::QCircuitMeasurementStep::Type::MEASURE_ND:
                    case internal::QCircuitMeasurementStep::Type::
                        MEASURE_MANY_ND:
                    case internal::QCircuitMeasurementStep::Type::MEASURE_V_ND:
                    case internal::QCircuitMeasurementStep::Type::
                        MEASURE_V_JOINT_ND:
                        os << "|] ";
                        break;
                    case internal::QCircuitMeasurementStep::Type::RESET:
                    case internal::QCircuitMeasurementStep::Type::RESET_MANY:
                        os << "|* ";
                        break;
                    case internal::QCircuitMeasurementStep::Type::DISCARD:
                    case internal::QCircuitMeasurementStep::Type::DISCARD_MANY:
                        os << "|x ";
                        break;
                    case internal::QCircuitMeasurementStep::Type::POST_SELECT:
                    case internal::QCircuitMeasurementStep::Type::
                        POST_SELECT_MANY:
                    case internal::QCircuitMeasurementStep::Type::POST_SELECT_V:
                    case internal::QCircuitMeasurementStep::Type::
                        POST_SELECT_V_JOINT:
                        os << "<| ";
                        break;
                    case internal::QCircuitMeasurementStep::Type::
                        POST_SELECT_ND:
                    case internal::QCircuitMeasurementStep::Type::
                        POST_SELECT_MANY_ND:
                    case internal::QCircuitMeasurementStep::Type::
                        POST_SELECT_V_ND:
                    case internal::QCircuitMeasurementStep::Type::
                        POST_SELECT_V_JOINT_ND:
                        os << "[| ";
                        break;
                } /* end switch */
                os << measurement_step;
            };
            auto nop_step_visitor =
                [&](const internal::QCircuitNOPStep& nop_step) {
                    os << nop_step;
                };
            auto conditional_step_visitor =
                [&](const internal::QCircuitConditionalStep& conditional_step) {
                    os << conditional_step;
                };

            qc_ptr_->visit_circuit_step_(
                qc_ptr_->circuit_[ip_], conditional_step_visitor,
                gate_step_visitor, measurement_step_visitor, nop_step_visitor);

            return os;
        }
    }; /* class QCircuit::QCircuitIterator::value_type_ */

  private:
    ///< non-owning pointer to the parent const quantum circuit
    ///< description
    const QCircuit* qc_ptr_{nullptr};
    idx ip_{static_cast<idx>(-1)}; ///< instruction pointer

  public:
    // iterator traits
    using value_type = value_type_;                      ///< iterator trait
    using pointer = const value_type*;                   ///< iterator trait
    using reference = const value_type&;                 ///< iterator trait
    using difference_type = std::ptrdiff_t;              ///< iterator trait
    using iterator_category = std::forward_iterator_tag; ///< iterator trait

    /**
     * \brief Default constructor
     */
    QCircuitIterator() = default;

    /**
     * \brief Explicit constructor
     *
     * \param qc_ptr Pointer to underlying quantum circuit description
     * \param ip Quantum circuit description instruction pointer
     */
    explicit QCircuitIterator(const QCircuit* qc_ptr, idx ip)
        : qc_ptr_{qc_ptr}, ip_{ip} {}

    // silence -Weffc++ class has pointer data members
    /**
     * \brief Default copy constructor
     */
    QCircuitIterator(const QCircuitIterator&) = default;

    // silence -Weffc++ class has pointer data members
    /**
     * \brief Default copy assignment operator
     *
     * \return Reference to the current instance
     */
    QCircuitIterator& operator=(const QCircuitIterator&) = default;

    /**
     * \brief Prefix increment operator
     *
     * \return Reference to the current instance
     */
    QCircuitIterator& operator++() {
        // EXCEPTION CHECKS
        // protects against incrementing invalid iterators
        if (qc_ptr_ == nullptr) {
            throw exception::InvalidIterator(
                "qpp::internal::QCircuitIterator::operator++()",
                "No qpp::QCircuit assigned");
        }

        idx num_steps = qc_ptr_->get_step_count();
        // protects against incrementing an empty circuit iterator
        if (num_steps == 0) {
            throw exception::InvalidIterator(
                "qpp::internal::QCircuitIterator::operator++()",
                "Zero-sized qpp::QCircuit");
        }

        // protects against incrementing past the end
        if (ip_ == num_steps) {
            throw exception::InvalidIterator(
                "qpp::internal::QCircuitIterator::operator++()",
                "Incrementing past the end");
        }
        // END EXCEPTION CHECKS

        ++ip_;

        return *this;
    }

    /**
     * \brief Postfix increment operator
     *
     * \return Copy of the current instance before the increment
     */
    QCircuitIterator operator++(int) {
        QCircuitIterator retval{*this};
        this->operator++();

        return retval;
    }

    /**
     * \brief Equality operator
     *
     * \param rhs Iterator against which the equality is being tested
     * \return True if the iterators are equal, false otherwise
     */
    bool operator==(const QCircuitIterator& rhs) const noexcept {
        return std::tie(ip_, qc_ptr_) == std::tie(rhs.ip_, rhs.qc_ptr_);
    }

    /**
     * \brief Inequality operator
     *
     * \param rhs Iterator against which the inequality is being tested
     * \return True if the iterators are not equal (bitwise), false
     * otherwise
     */
    bool operator!=(const QCircuitIterator& rhs) const noexcept {
        return !(*this == rhs);
    }

    /**
     * \brief Safe de-referencing operator
     *
     * \return De-referenced iterator element
     */
    value_type operator*() const {
        // EXCEPTION CHECKS
        // protects against de-referencing past the last element or
        // against de-referencing invalid iterators
        idx num_steps = qc_ptr_->get_step_count();
        if (qc_ptr_ == nullptr) {
            throw exception::InvalidIterator(
                "qpp::internal::QCircuitIterator::operator*()",
                "No qpp::QCircuit assigned");
        }
        if (num_steps == 0) {
            throw exception::InvalidIterator(
                "qpp::internal::QCircuitIterator::operator*()",
                "Zero-sized qpp::QCircuit");
        }
        if (ip_ == num_steps) {
            throw exception::InvalidIterator(
                "qpp::internal::QCircuitIterator::operator*()",
                "Dereferencing past the end");
        }
        // END EXCEPTION CHECKS

        return value_type{qc_ptr_, ip_, qc_ptr_->circuit_[ip_]};
    }
}; /* class internal::QCircuitIterator */
} /* namespace internal */

inline QCircuit::iterator QCircuit::begin() noexcept {
    return std::as_const(*this).begin();
}

inline QCircuit::iterator QCircuit::end() noexcept {
    return std::as_const(*this).end();
}

inline QCircuit::const_iterator QCircuit::begin() const noexcept {
    idx ip = get_step_count() != 0 ? 0 : static_cast<idx>(-1);

    return const_iterator{this, ip};
}

inline QCircuit::const_iterator QCircuit::end() const noexcept {
    idx step_count = this->get_step_count();
    idx ip = std::numeric_limits<idx>::max();
    if (step_count != 0) {
        ip = step_count;
    }

    return const_iterator{this, ip};
}

inline QCircuit::const_iterator QCircuit::cbegin() const noexcept {
    return std::as_const(*this).begin();
}

inline QCircuit::const_iterator QCircuit::cend() const noexcept {
    return std::as_const(*this).end();
}

inline std::ostream& QCircuit::display(std::ostream& os) const {
    os << "[QCircuit ";
    os << "nq: " << nq_ << ", nc: " << nc_ << ", d: " << d_;
    if (name_.has_value()) {
        os << ", name: ";
        os << std::quoted(name_.value());
    }
    os << "]\n";

    std::string sep{};
    for (auto&& elem : *this) {
        os << sep << elem;
        sep = '\n';
    }

    return os;
}

inline std::string QCircuit::to_JSON(bool enclosed_in_curly_brackets) const {
    std::string result;

    if (enclosed_in_curly_brackets) {
        result += "{";
    }

    if (name_.has_value()) {
        result += "\"name\": \"" + name_.value() + "\", ";
    }

    std::string sep;
    std::ostringstream ss;
    result += "\"steps\": [";
    for (auto&& elem : *this) {
        result += sep;
        sep = ", ";
        result += "{\"step\": " + std::to_string(elem.get_ip()) + ", ";
        result += "\"type\": ";

        auto gate_step_visitor =
            [&](const internal::QCircuitGateStep& gate_step) {
                ss.str("");
                ss.clear();
                ss << gate_step.gate_type_;
                result += '\"' + ss.str() + "\", ";
                if (gate_step.ctrl_.has_value() &&
                    !gate_step.ctrl_.value().empty()) {
                    ss.str("");
                    ss.clear();
                    ss << disp(gate_step.ctrl_.value(),
                               IOManipContainerOpts{}.set_sep(", "));
                    switch (gate_step.gate_type_) {
                        case internal::QCircuitGateStep::Type::CTRL:
                        case internal::QCircuitGateStep::Type::CTRL_FAN:
                            result += "\"ctrl\": " + ss.str() + ", ";
                            break;
                        case internal::QCircuitGateStep::Type::cCTRL:
                        case internal::QCircuitGateStep::Type::cCTRL_FAN:
                            result += "\"c_ctrl\": " + ss.str() + ", ";
                            break;
                        default:
                            break;
                    }
                }
                ss.str("");
                ss.clear();
                ss << disp(gate_step.target_,
                           IOManipContainerOpts{}.set_sep(", "));
                result += "\"target\": " + ss.str();

                if (gate_step.shift_.has_value()) {
                    result += ", \"shift\": ";
                    ss.str("");
                    ss.clear();
                    ss << disp(gate_step.shift_.value(),
                               IOManipContainerOpts{}.set_sep(", "));
                    result += ss.str();
                }

                if (gate_step.name_.has_value()) {
                    result += ", \"name\": \"" + gate_step.name_.value() + "\"";
                }
                result += "}";
            };
        auto measurement_step_visitor =
            [&](const internal::QCircuitMeasurementStep& measurement_step) {
                ss.str("");
                ss.clear();
                ss << measurement_step.measurement_type_;
                result += '\"' + ss.str() + "\", ";
                ss.str("");
                ss.clear();
                ss << disp(measurement_step.target_,
                           IOManipContainerOpts{}.set_sep(", "));
                result += "\"target\": " + ss.str();

                if (measurement_step.ps_vals_.has_value()) {
                    ss.str("");
                    ss.clear();
                    ss << disp(measurement_step.ps_vals_.value(),
                               IOManipContainerOpts{}.set_sep(", "));
                    result += ", \"ps_vals\": " + ss.str();
                }

                if (measurement_step.measurement_type_ !=
                        internal::QCircuitMeasurementStep::Type::RESET &&
                    measurement_step.measurement_type_ !=
                        internal::QCircuitMeasurementStep::Type::RESET_MANY &&
                    measurement_step.measurement_type_ !=
                        internal::QCircuitMeasurementStep::Type::DISCARD &&
                    measurement_step.measurement_type_ !=
                        internal::QCircuitMeasurementStep::Type::DISCARD_MANY) {
                    result += ", \"c_reg\": " +
                              std::to_string(measurement_step.c_reg_);
                }

                if (measurement_step.name_.has_value()) {
                    result += ", \"name\": \"" +
                              measurement_step.name_.value() + "\"";
                }
                result += "}";
            };
        auto nop_step_visitor = [&](const internal::QCircuitNOPStep&) {
            result += std::string{"\"NOP\""} + "}";
        };
        auto conditional_step_visitor =
            [&](const internal::QCircuitConditionalStep&) {};

        visit_circuit_step_(elem.get_step(), conditional_step_visitor,
                            gate_step_visitor, measurement_step_visitor,
                            nop_step_visitor);
    } // end for
    result += "], "; // end steps

    result += get_resources().to_JSON(false) + ", ";

    ss.str("");
    ss.clear();
    ss << disp(get_measured_d(), IOManipContainerOpts{}.set_sep(", "));
    result += "\"measured/discarded (destructive)\": " + ss.str() + ", ";

    ss.str("");
    ss.clear();
    ss << disp(get_measured_nd(), IOManipContainerOpts{}.set_sep(", "));
    result += "\"measured (non-destructive)\": " + ss.str() + ", ";

    ss.str("");
    ss.clear();
    ss << disp(get_measurement_dits(), IOManipContainerOpts{}.set_sep(", "));
    result += "\"measurement dits\": " + ss.str();

    if (enclosed_in_curly_brackets) {
        result += "}";
    }

    return result;
} /* QCircuit::to_JSON() */

/**
 * \class qpp::QCircuitTraits
 * \brief Generic type traits for classes that describe quantum circuits
 */
template <typename T>
struct QCircuitTraits;

/**
 * \class qpp::QCircuitTraits<QCircuit>
 * \brief Specialization type traits for qpp::QCircuit
 */
template <>
struct QCircuitTraits<QCircuit> {
    using iterator_type = QCircuit::iterator;
    using value_type = QCircuit::iterator::value_type;
};

// QCircuit free functions

/**
 * \brief Composes (appends) a quantum circuit description to the end of
 * another one
 * \see qpp::couple_circuit_left() and qpp::couple_circuit_right()
 *
 * \note If qudit indexes of the second quantum circuit description do
 * not totally overlap with the indexes of the first quantum circuit
 * description, then the required number of additional qudits are
 * automatically added to the output quantum circuit description
 *
 * \param qc1 Quantum circuit description
 * \param qc2 Quantum circuit description
 * \param pos_qudit The index of the first/top qudit of \a qc2 quantum
 * circuit description relative to the index of the first/top qudit of the
 * \a qc1 quantum circuit description, with the rest following in order.
 * If negative or greater than the total number of qudits of \a qc1,
 * then the required number of additional qudits are automatically added
 * to the output quantum circuit description.
 * \param name Optional result's circuit name
 * \param pos_dit Optional, the first classical dit of \a qc2 quantum
 * circuit description is inserted before the \a pos_dit classical dit index
 * of the \a qc1 quantum circuit description (in the classical dits array),
 * the rest following in order. If absent (default), insertion is performed
 * at the end.
 * \return Combined quantum circuit description, with \a qc2 added at the
 * end of \a qc1
 */
inline QCircuit compose_circuit(QCircuit qc1, const QCircuit& qc2,
                                bigint pos_qudit,
                                std::optional<std::string> name = std::nullopt,
                                std::optional<idx> pos_dit = std::nullopt) {
    // EXCEPTION CHECKS
    // check equal dimensions
    if (qc1.get_d() != qc2.get_d()) {
        throw exception::DimsNotEqual("qpp::compose_circuit()", "other");
    }
    // check classical dits
    if (!pos_dit.has_value()) {
        pos_dit = qc1.get_nc();
    } else {
        if (internal::is_negative(pos_dit.value()) ||
            pos_dit.value() > qc1.get_nc()) {
            throw exception::OutOfRange("qpp::compose_circuit()", "pos_dit");
        }
    }
    // check that overlapping qudits (in the current instance) were not
    // already destructively measured
    if (pos_qudit < 0 && (pos_qudit + static_cast<bigint>(qc2.get_nq())) >= 0) {
        for (idx i = 0;
             i < std::min(static_cast<idx>(pos_qudit +
                                           static_cast<bigint>(qc2.get_nq())),
                          qc1.get_nq());
             ++i) {
            if (qc1.was_measured_d(i)) {
                throw exception::QuditAlreadyMeasured("qpp::compose_circuit()",
                                                      "qc1");
            }
        }
    }
    if (pos_qudit >= 0 && static_cast<idx>(pos_qudit) < qc1.get_nq()) {
        for (idx i = 0; i < std::min(static_cast<idx>(qc1.get_nq() - pos_qudit),
                                     qc2.get_nq());
             ++i) {
            if (qc1.was_measured_d(pos_qudit + i)) {
                throw exception::QuditAlreadyMeasured("qpp::compose_circuit()",
                                                      "qc1");
            }
        }
    }
    // END EXCEPTION CHECKS

    if (name.has_value()) {
        qc1.set_name(name.value());
    }

    return qc1.compose_circuit(qc2, pos_qudit, pos_dit);
}

/**
 * \brief Couples in-place a quantum circuit description \a qc2 to another
 * quantum circuit description \a qc1, with the \a qc2 quantum circuit
 * description placed at the left (beginning) of the first quantum circuit
 * description
 * \see qpp::couple_circuit_right() and qpp::compose_circuit()
 *
 * \note The added quantum circuit description \a qc2 cannot be larger than
 * the \a qc1 quantum circuit description, i.e., all qudit indexes of the
 * added quantum circuit description must match with qudits from the \a qc1
 * quantum circuit description (and those coupled of the latter must contain
 * no measurements)
 *
 * \note The classical dits are not relabeled
 *
 * \param qc1 Quantum circuit description
 * \param qc2 Quantum circuit description
 * \param target Qudit indexes of the \a qc1 circuit description where the
 * qudits of \a qc2 are being coupled, i.e., the first/top qudit of
 * \a qc2 quantum circuit description is coupled with the target[0] qudit
 * of the \a qc1 circuit description, and so on
 * \param name Optional result's circuit name
 * \param pos_dit Optional, the first classical dit of \a qc2 quantum
 * circuit description is inserted before the \a pos_dit classical dit index
 * of the \a qc1 quantum circuit description (in the classical dits array),
 * the rest following in order. If absent (default), insertion is performed
 * at the end.
 * \return Combined quantum circuit description
 */
inline QCircuit
couple_circuit_left(QCircuit qc1, const QCircuit& qc2,
                    const std::vector<idx>& target,
                    std::optional<std::string> name = std::nullopt,
                    std::optional<idx> pos_dit = std::nullopt) {
    // EXCEPTION CHECKS
    // check equal dimensions
    if (qc1.get_d() != qc2.get_d()) {
        throw exception::DimsNotEqual("qpp::couple_circuit_left()", "qc1/qc2");
    }
    // check classical dits
    if (!pos_dit.has_value()) {
        pos_dit = qc1.get_nc();
    } else if (internal::is_negative(pos_dit.value()) ||
               pos_dit.value() > qc1.get_nc()) {
        throw exception::OutOfRange("qpp::couple_circuit_left()", "pos_dit");
    }
    // check no measurement for the coupled circuit
    if (!qc2.get_measured_d().empty()) {
        throw exception::QuditAlreadyMeasured("qpp::couple_circuit_left()",
                                              "qc2");
    }
    // check valid target
    if (static_cast<idx>(target.size()) != qc1.get_nq()) {
        throw exception::OutOfRange("qpp::couple_circuit_left()", "target");
    }
    if (static_cast<idx>(target.size()) > qc2.get_nq()) {
        throw exception::OutOfRange("qpp::couple_circuit_left()", "target");
    }
    if (!internal::check_no_duplicates(target)) {
        throw exception::Duplicates("qpp::couple_circuit_left()", "target");
    }
    for (idx elem : target) {
        if (elem >= qc1.get_nq()) {
            throw exception::OutOfRange("qpp::couple_circuit_left()", "target");
        }
    }
    // check matching qudits (in the current instance) were not already
    // measured destructively
    for (idx elem : target) {
        if (qc1.was_measured_d(elem)) {
            throw exception::QuditAlreadyMeasured("qpp::couple_circuit_left()",
                                                  "target");
        }
    }
    // END EXCEPTION CHECKS

    if (name.has_value()) {
        qc1.set_name(name.value());
    }

    return qc1.couple_circuit_left(qc2, target, pos_dit);
}

/**
 * \brief Couples in-place a quantum circuit description \a qc2 to another
 * quantum circuit description \a qc1, with the \a qc2 quantum circuit
 * description placed at the right (end) of the first quantum circuit
 * description
 * \see qpp::couple_circuit_left() and qpp::compose_circuit()
 *
 * \note The added quantum circuit description \a qc2 cannot be larger than
 * the \a qc1 quantum circuit description, i.e., all qudit indexes of the
 * added quantum circuit description must match with qudits from the \a qc1
 * quantum circuit description (and those coupled of the latter must contain
 * no measurements)
 *
 * \note The classical dits are not relabeled
 *
 * \param qc1 Quantum circuit description
 * \param qc2 Quantum circuit description
 * \param target Qudit indexes of the \a qc1 circuit description where the
 * qudits of \a qc2 are being coupled, i.e., the first/top qudit of \a qc2
 * quantum circuit description is coupled with the target[0] qudit of the
 * \a qc1 circuit description, and so on
 * \param name Optional result's circuit name
 * \param pos_dit Optional, the first classical dit of \a qc2 quantum
 * circuit description is inserted before the \a pos_dit classical dit index
 * of the \a qc1 quantum circuit description (in the classical dits array),
 * the rest following in order. If absent (default), insertion is performed
 * at the end.
 * \return Combined quantum circuit description
 */
inline QCircuit
couple_circuit_right(QCircuit qc1, const QCircuit& qc2,
                     const std::vector<idx>& target,
                     std::optional<std::string> name = std::nullopt,
                     std::optional<idx> pos_dit = std::nullopt) {

    // EXCEPTION CHECKS
    // check equal dimensions
    if (qc1.get_d() != qc2.get_d()) {
        throw exception::DimsNotEqual("qpp::couple_circuit_right()", "qc1/qc2");
    }
    // check classical dits
    if (!pos_dit.has_value()) {
        pos_dit = qc1.get_nc();
    } else if (internal::is_negative(pos_dit.value()) ||
               pos_dit.value() > qc1.get_nc()) {
        throw exception::OutOfRange("qpp::couple_circuit_right()", "pos_dit");
    }
    // check valid target
    if (static_cast<idx>(target.size()) != qc1.get_nq()) {
        throw exception::OutOfRange("qpp::couple_circuit_right()", "target");
    }
    if (static_cast<idx>(target.size()) > qc2.get_nq()) {
        throw exception::OutOfRange("qpp::couple_circuit_right()", "target");
    }
    if (!internal::check_no_duplicates(target)) {
        throw exception::Duplicates("qpp::couple_circuit_right()", "target");
    }
    for (idx elem : target) {
        if (elem >= qc1.get_nq()) {
            throw exception::OutOfRange("qpp::couple_circuit_right()",
                                        "target");
        }
    }
    // check matching qudits (in the current instance) were not already
    // measured destructively
    for (idx elem : target) {
        if (qc1.was_measured_d(elem)) {
            throw exception::QuditAlreadyMeasured("qpp::couple_circuit_right()",
                                                  "target");
        }
    }
    // END EXCEPTION CHECKS

    if (name.has_value()) {
        qc1.set_name(name.value());
    }

    return qc1.couple_circuit_right(qc2, target, pos_dit);
}

/**
 * \brief Composes (appends) the \a qc_target controlled quantum circuit
 * description to the end of the \a qc_ctrl quantum circuit description;
 * \a qc_ctrl controls the \a qc_target.
 * \see qpp::compose_circuit()
 *
 * \note For each gate or controlled-gate in the target circuit description,
 * a control is being added from the control circuit description control
 * qudits. Classical control gates and measurements/resets/discards in the
 * target circuit description are left unchanged.
 *
 * \note If qudit indexes of the second quantum circuit description do
 * not totally overlap with the indexes of the first quantum circuit
 * description, then the required number of additional qudits are
 * automatically added to the output quantum circuit description
 *
 * \param qc_ctrl Control quantum circuit description
 * \param ctrl Control qubits
 * \param qc_target Target quantum circuit description
 * \param pos_qudit The index of the first/top qudit of \a qc_target quantum
 * circuit description relative to the index of the first/top qudit of the
 * \a qc_ctrl quantum circuit description, with the rest following in order.
 * If negative or greater than the total number of qudits of \a qc_ctrl,
 * then the required number of additional qudits are automatically added
 * to the output quantum circuit description.
 * \param shift Optional, performs the control as if the \a ctrl qudit state
 * was \f$X\f$-incremented by \a shift
 * \param pos_dit Optional, the first classical dit of \a qc_target quantum
 * circuit description is inserted before the \a pos_dit classical dit index
 * of the \a qc_ctrl quantum circuit description (in the classical dits
 * array), the rest following in order. If absent (default), insertion is
 * performed at the end.
 * \param name Optional result's circuit name
 * \return Combined quantum circuit description, with \a qc_target added at
 * the end of \a qc_ctrl
 */
inline QCircuit
compose_CTRL_circuit(QCircuit qc_ctrl, const std::vector<idx>& ctrl,
                     const QCircuit& qc_target, bigint pos_qudit,
                     std::optional<std::vector<idx>> shift = std::nullopt,
                     std::optional<idx> pos_dit = std::nullopt,
                     std::optional<std::string> name = std::nullopt) {

    // EXCEPTION CHECKS
    // check non-empty qc_ctrl and qc_target
    if (qc_ctrl.get_nq() == 0) {
        throw exception::ZeroSize("qpp::compose_CTRL_circuit()", "qc_ctrl");
    }
    if (qc_target.get_nq() == 0) {
        throw exception::ZeroSize("qpp::compose_CTRL_circuit()", "qc_target");
    }

    // check equal dimensions
    if (qc_ctrl.get_d() != qc_target.get_d()) {
        throw exception::DimsNotEqual("qpp::compose_CTRL_circuit()",
                                      "qc_ctrl/qc_target");
    }

    // check classical dits
    if (!pos_dit.has_value()) {
        pos_dit = qc_ctrl.get_nc();
    } else {
        if (internal::is_negative(pos_dit.value()) ||
            pos_dit.value() > qc_ctrl.get_nc()) {
            throw exception::OutOfRange("qpp::compose_CTRL_circuit()",
                                        "pos_dit");
        }
    }

    // check that overlapping qudits (in the current instance) were not
    // already destructively measured
    if (pos_qudit < 0 &&
        (pos_qudit + static_cast<bigint>(qc_target.get_nq())) >= 0) {
        for (idx i = 0;
             i < std::min(static_cast<idx>(pos_qudit + static_cast<bigint>(
                                                           qc_target.get_nq())),
                          qc_ctrl.get_nq());
             ++i) {
            if (qc_ctrl.was_measured_d(i)) {
                throw exception::QuditAlreadyMeasured(
                    "qpp::compose_CTRL_circuit()", "qc_ctrl");
            }
        }
    }
    if (pos_qudit >= 0 && static_cast<idx>(pos_qudit) < qc_ctrl.get_nq()) {
        for (idx i = 0;
             i < std::min(static_cast<idx>(qc_ctrl.get_nq() - pos_qudit),
                          qc_target.get_nq());
             ++i) {
            if (qc_ctrl.was_measured_d(pos_qudit + i)) {
                throw exception::QuditAlreadyMeasured(
                    "qpp::compose_CTRL_circuit()", "qc_ctrl");
            }
        }
    }

    // check valid ctrl
    if (ctrl.empty()) {
        throw exception::ZeroSize("qpp::compose_CTRL_circuit()", "ctrl");
    }
    for (idx elem : ctrl) {
        if (elem >= qc_ctrl.get_nq()) {
            throw exception::OutOfRange("qpp::compose_CTRL_circuit()", "ctrl");
        }
        // check ctrl was not measured before
        if (qc_ctrl.was_measured_d(elem)) {
            throw exception::QuditAlreadyMeasured("qpp::compose_CTRL_circuit()",
                                                  "ctrl");
        }
    }
    // check no duplicates ctrl
    if (!internal::check_no_duplicates(ctrl)) {
        throw exception::Duplicates("qpp::compose_CTRL_circuit()", "ctrl");
    }
    // check that ctrl and target do not overlap
    idx target_begin = pos_qudit;                    // including it
    idx target_end = pos_qudit + qc_target.get_nq(); // excluding it
    for (idx elem : ctrl) {
        if (elem >= target_begin && elem < target_end) {
            throw exception::OutOfRange("qpp::compose_CTRL_circuit()",
                                        "ctrl/qc_target");
        }
    }

    // check shift
    if (shift.has_value() && (shift.value().size() != ctrl.size())) {
        throw exception::SizeMismatch("qpp::compose_CTRL_circuit()",
                                      "ctrl/shift");
    }
    if (shift.has_value()) {
        for (idx elem : shift.value()) {
            if (elem >= qc_ctrl.get_d()) {
                throw exception::OutOfRange("qpp::compose_CTRL_circuit()",
                                            "shift");
            }
        }
    }
    // END EXCEPTION CHECKS

    if (name.has_value()) {
        qc_ctrl.set_name(name.value());
    }

    return qc_ctrl.compose_CTRL_circuit(ctrl, qc_target, pos_qudit, shift,
                                        pos_dit);
}

//  TODO: check for reset/measured non-destructively etc.
/**
 * \brief Adjoint quantum circuit description
 *
 * \param qc Quantum circuit description
 * \param name Optional result's circuit name
 * \return Adjoint quantum circuit description
 */
inline QCircuit adjoint(QCircuit qc,
                        std::optional<std::string> name = nullptr) {
    // EXCEPTION CHECKS
    if (qc.has_measurements()) {
        throw exception::QuditAlreadyMeasured("qpp::adjoint()", "qc");
    }
    // END EXCEPTION CHECKS

    if (name.has_value()) {
        qc.set_name(name.value());
    }

    return qc.adjoint();
}

/**
 * \brief Kronecker product between two quantum circuit descriptions
 *
 * \param qc1 Quantum circuit description
 * \param qc2 Quantum circuit description
 * \param name Optional result's circuit name
 * \return Quantum circuit description of the Kronecker product of \a qc1
 * with \a qc2
 */
inline QCircuit kron(QCircuit qc1, const QCircuit& qc2,
                     std::optional<std::string> name = std::nullopt) {
    if (name.has_value()) {
        qc1.set_name(name.value());
    }

    return qc1.kron(qc2);
}

/**
 * \brief Replicates a quantum circuit description
 * \note The circuit should not contain any steps that remove qudits
 *
 * \param qc Quantum circuit description
 * \param n Number of repetitions. If \a n == 1, returns the original
 * circuit.
 * \param name Optional result's circuit name
 * \return Replicated quantum circuit description
 */
inline QCircuit replicate(QCircuit qc, idx n,
                          std::optional<std::string> name = std::nullopt) {
    if (name.has_value()) {
        qc.set_name(name.value());
    }

    // EXCEPTION CHECKS
    if (n == 0) {
        throw exception::OutOfRange("qpp::replicate()", "n");
    }
    if (qc.removes_qudits()) {
        throw exception::QuditAlreadyMeasured("qpp::replicate()", "qc");
    }
    if (n == 1) {
        return qc;
    }
    // END EXCEPTION CHECKS

    return qc.replicate(n);
}

/**
 * \brief Random quantum circuit description generator for fixed gate count
 *
 * \param nq Number of qudits
 * \param d Subsystem dimensions
 * \param gate_count Circuit gate count
 * \param p_two Optional, probability of applying a two qudit gate, must belong
 * to the interval (0,1). If the two qudit gate set has more than one element,
 * then the gate is chosen at random from the set.
 * \param with_respect_to_gate Optional, if present, gate count is calculated
 * with respect to this particular gate (absent by default, so by default gate
 * count is calculated with respect to all gates in the circuit)
 * \param one_qudit_gate_set Optional set of one qudit gates, must be specified
 * for \a d > 2
 * \param two_qudit_gate_set Optional set of two qudit gates, must be specified
 * for \a d > 2
 * \param one_qudit_gate_names Optional one qudit gate names
 * \param two_qudit_gate_names Optional two qudit gate names
 * \return Instance of random qpp::QCircuit for fixed gate count
 */
inline QCircuit random_circuit_count(
    idx nq, idx d, idx gate_count, std::optional<realT> p_two,
    std::optional<cmat> with_respect_to_gate = std::nullopt,
    std::optional<std::vector<cmat>> one_qudit_gate_set = std::nullopt,
    std::optional<std::vector<cmat>> two_qudit_gate_set = std::nullopt,
    std::optional<std::vector<std::string>> one_qudit_gate_names = std::nullopt,
    std::optional<std::vector<std::string>> two_qudit_gate_names =
        std::nullopt) {
    // EXCEPTION CHECKS
    // check valid dimension
    if (d < 2) {
        throw exception::DimsInvalid("qpp::random_circuit_count()", "d");
    }
    // check valid probabilities
    if (p_two.has_value()) {
        if (!(p_two.value() > 0 && p_two.value() < 1)) {
            throw exception::OutOfRange("qpp::random_circuit_count()", "p_two");
        }
    }
    if (nq < 1 || (p_two.has_value() && nq == 1)) {
        throw exception::OutOfRange("qpp::random_circuit_count()", "nq/p_two");
    }
    // pre-fill gate sets for qubits (in case they're empty)
    if (d == 2 && !one_qudit_gate_set.has_value()) {
        one_qudit_gate_set = std::vector{
            Gates::get_instance().X, Gates::get_instance().Y,
            Gates::get_instance().Z, Gates::get_instance().H,
            Gates::get_instance().S, adjoint(Gates::get_instance().S),
            Gates::get_instance().T, adjoint(Gates::get_instance().T)};
    }
    if (d == 2 && !two_qudit_gate_set.has_value()) {
        two_qudit_gate_set = std::vector{Gates::get_instance().CNOT};
    }
    // check gate sets are not empty for d > 2
    if (d > 2) {
        if (!one_qudit_gate_set.has_value()) {
            throw exception::ZeroSize("qpp::random_circuit_count()",
                                      "one_qudit_gate_set");
        }
        if (p_two.has_value() && !two_qudit_gate_set.has_value()) {
            throw exception::ZeroSize("qpp::random_circuit_count()",
                                      "two_qudit_gate_set");
        }
    }
    // check one qudit gate name sizes
    if (one_qudit_gate_names.has_value()) {
        if (!one_qudit_gate_set.has_value() ||
            one_qudit_gate_names.value().size() !=
                one_qudit_gate_set.value().size()) {
            throw exception::SizeMismatch("qpp::random_circuit_count()",
                                          "one_qudit_gate_names");
        }
    }
    // check two qudit gate name sizes
    if (two_qudit_gate_names.has_value()) {
        if (!two_qudit_gate_set.has_value() ||
            two_qudit_gate_names.value().size() !=
                two_qudit_gate_set.value().size()) {
            throw exception::SizeMismatch("qpp::random_circuit_count()",
                                          "two_qudit_gate_names");
        }
    }
    // check one qudit gate sizes
    if (one_qudit_gate_set.has_value()) {
        for (auto&& gate : one_qudit_gate_set.value()) {
            if (!internal::check_square_mat(gate)) {
                throw exception::MatrixNotSquare("qpp::random_circuit_count()",
                                                 "one_qudit_gate_set");
            }
            if (static_cast<idx>(gate.rows()) != d) {
                throw exception::MatrixMismatchSubsys(
                    "qpp::random_circuit_count()", "one_qudit_gate_set");
            }
        }
    }
    // check two qudit gate sizes
    if (two_qudit_gate_set.has_value()) {
        for (auto&& gate : two_qudit_gate_set.value()) {
            if (!internal::check_square_mat(gate)) {
                throw exception::MatrixNotSquare("qpp::random_circuit_count()",
                                                 "two_qudit_gate_set");
            }
            if (static_cast<idx>(gate.rows()) != d * d) {
                throw exception::MatrixMismatchSubsys(
                    "qpp::random_circuit_count()", "two_qudit_gate_set");
            }
        }
    }
    // check with_respect_to_gate
    if (with_respect_to_gate.has_value()) {
        cmat wrt_gate = with_respect_to_gate.value();
        auto rows = static_cast<idx>(wrt_gate.rows());
        // check square matrix
        if (!internal::check_square_mat(wrt_gate)) {
            throw exception::MatrixNotSquare("qpp::random_circuit_count()",
                                             "with_respect_to_gate");
        }
        // check size (either 1 qubit or 2 qubit gate)
        bool is_one_qubit_gate = (rows == d);
        bool is_two_qubit_gate = (rows == d * d);
        if (!(is_one_qubit_gate || is_two_qubit_gate)) {
            throw exception::MatrixMismatchSubsys("qpp::random_circuit_count()",
                                                  "with_respect_to_gate");
        }
        // check with_respect_to_gate is present in the gate sets
        if (is_one_qubit_gate) {
            bool found =
                std::find(one_qudit_gate_set.value().begin(),
                          one_qudit_gate_set.value().end(),
                          wrt_gate) != one_qudit_gate_set.value().end();
            if (!found) {
                throw exception::NotFound(
                    "qpp::random_circuit_count()",
                    "with_respect_to_gate/one_qudit_gate_set");
            }
        }
        if (is_two_qubit_gate) {
            bool found =
                std::find(two_qudit_gate_set.value().begin(),
                          two_qudit_gate_set.value().end(),
                          wrt_gate) != two_qudit_gate_set.value().end();
            if (!found) {
                throw exception::NotFound(
                    "qpp::random_circuit_count()",
                    "with_respect_to_gate/two_qudit_gate_set");
            }
        }
    } // end if (with_respect_to_gate.has_value())
    // END EXCEPTION CHECKS

    realT p_one = p_two.has_value() ? 1 - p_two.value() : 1;
    QCircuit qc{nq, 0, d};
    idx current_count = 0;
    while (current_count < gate_count) {
        bool is_one_qudit_gate = bernoulli(p_one);
        if (is_one_qudit_gate) {
            idx q = randidx(0, nq - 1);
            idx gate = randidx(0, one_qudit_gate_set.value().size() - 1);
            if (!one_qudit_gate_names.has_value()) {
                qc.gate(one_qudit_gate_set.value()[gate], q);
            } else {
                qc.gate(one_qudit_gate_set.value()[gate], q,
                        one_qudit_gate_names.value()[gate]);
            }
        } else {
            idx ctrl = randidx(0, nq - 1);
            idx target = randidx(0, nq - 1);
            while (ctrl == target) {
                target = randidx(0, nq - 1);
            }
            idx gate = randidx(0, two_qudit_gate_set.value().size() - 1);
            if (!two_qudit_gate_names.has_value()) {
                qc.gate(two_qudit_gate_set.value()[gate], ctrl, target);
            } else {
                qc.gate(two_qudit_gate_set.value()[gate], ctrl, target,
                        two_qudit_gate_names.value()[gate]);
            }
        }
        current_count = with_respect_to_gate.has_value()
                            ? qc.get_gate_count(with_respect_to_gate.value())
                            : qc.get_gate_count();
    }

    return qc;
}

/**
 * \brief Random quantum circuit description generator for fixed gate depth
 *
 * \param nq Number of qudits
 * \param d Subsystem dimensions
 * \param gate_depth Circuit gate depth
 * \param p_two Optional, probability of applying a two qudit gate, must belong
 * to the interval (0,1). If the two qudit gate set has more than one element,
 * then the gate is chosen at random from the set.
 * \param with_respect_to_gate Optional, if present, gate depth is calculated
 * with respect to this particular gate (absent by default, so by default gate
 * depth is calculated with respect to all gates in the circuit)
 * \param one_qudit_gate_set Set of one qudit gates (optional, must be specified
 * for \a d > 2)
 * \param two_qudit_gate_set Set of two qudit gates (optional, must be specified
 * for \a d > 2);
 * \param one_qudit_gate_names One qudit gate names (optional)
 * \param two_qudit_gate_names Two qudit gate names (optional)
 * \return Instance of random qpp::QCircuit for fixed circuit gate depth
 */
inline QCircuit random_circuit_depth(
    idx nq, idx d, idx gate_depth, std::optional<realT> p_two,
    std::optional<cmat> with_respect_to_gate = std::nullopt,
    std::optional<std::vector<cmat>> one_qudit_gate_set = std::nullopt,
    std::optional<std::vector<cmat>> two_qudit_gate_set = std::nullopt,
    std::optional<std::vector<std::string>> one_qudit_gate_names = std::nullopt,
    std::optional<std::vector<std::string>> two_qudit_gate_names =
        std::nullopt) {
    // EXCEPTION CHECKS
    // check valid dimension
    if (d < 2) {
        throw exception::DimsInvalid("qpp::random_circuit_depth()", "d");
    }
    // check valid probabilities
    if (p_two.has_value()) {
        if (!(p_two.value() > 0 && p_two.value() < 1)) {
            throw exception::OutOfRange("qpp::random_circuit_depth()", "p_two");
        }
    }
    if (nq < 1 || (p_two.has_value() && nq == 1)) {
        throw exception::OutOfRange("qpp::random_circuit_depth()", "nq/p_two");
    }
    // pre-fill gate sets for qubits (in case they're empty)
    if (d == 2 && !one_qudit_gate_set.has_value()) {
        one_qudit_gate_set = std::vector{
            Gates::get_instance().X, Gates::get_instance().Y,
            Gates::get_instance().Z, Gates::get_instance().H,
            Gates::get_instance().S, adjoint(Gates::get_instance().S),
            Gates::get_instance().T, adjoint(Gates::get_instance().T)};
    }
    if (d == 2 && !two_qudit_gate_set.has_value()) {
        two_qudit_gate_set = std::vector{Gates::get_instance().CNOT};
    }
    // check gate sets are not empty for d > 2
    if (d > 2) {
        if (!one_qudit_gate_set.has_value()) {
            throw exception::ZeroSize("qpp::random_circuit_depth()",
                                      "one_qudit_gate_set");
        }
        if (p_two.has_value() && !two_qudit_gate_set.has_value()) {
            throw exception::ZeroSize("qpp::random_circuit_depth()",
                                      "two_qudit_gate_set");
        }
    }
    // check one qudit gate name sizes
    if (one_qudit_gate_names.has_value()) {
        if (!one_qudit_gate_set.has_value() ||
            one_qudit_gate_names.value().size() !=
                one_qudit_gate_set.value().size()) {
            throw exception::SizeMismatch("qpp::random_circuit_depth()",
                                          "one_qudit_gate_names");
        }
    }
    // check two qudit gate name sizes
    if (two_qudit_gate_names.has_value()) {
        if (!two_qudit_gate_set.has_value() ||
            two_qudit_gate_names.value().size() !=
                two_qudit_gate_set.value().size()) {
            throw exception::SizeMismatch("qpp::random_circuit_depth()",
                                          "two_qudit_gate_names");
        }
    }
    // check one qudit gate sizes
    if (one_qudit_gate_set.has_value()) {
        for (auto&& gate : one_qudit_gate_set.value()) {
            if (!internal::check_square_mat(gate)) {
                throw exception::MatrixNotSquare("qpp::random_circuit_depth()",
                                                 "one_qudit_gate_set");
            }
            if (static_cast<idx>(gate.rows()) != d) {
                throw exception::MatrixMismatchSubsys(
                    "qpp::random_circuit_depth()", "one_qudit_gate_set");
            }
        }
    }
    // check two qudit gate sizes
    if (two_qudit_gate_set.has_value()) {
        for (auto&& gate : two_qudit_gate_set.value()) {
            if (!internal::check_square_mat(gate)) {
                throw exception::MatrixNotSquare("qpp::random_circuit_depth()",
                                                 "two_qudit_gate_set");
            }
            if (static_cast<idx>(gate.rows()) != d * d) {
                throw exception::MatrixMismatchSubsys(
                    "qpp::random_circuit_depth()", "two_qudit_gate_set");
            }
        }
    }
    // check with_respect_to_gate
    if (with_respect_to_gate.has_value()) {
        cmat wrt_gate = with_respect_to_gate.value();
        auto rows = static_cast<idx>(wrt_gate.rows());
        // check square matrix
        if (!internal::check_square_mat(wrt_gate)) {
            throw exception::MatrixNotSquare("qpp::random_circuit_depth()",
                                             "with_respect_to_gate");
        }
        // check size (either 1 qubit or 2 qubit gate)
        bool is_one_qubit_gate = (rows == d);
        bool is_two_qubit_gate = (rows == d * d);
        if (!(is_one_qubit_gate || is_two_qubit_gate)) {
            throw exception::MatrixMismatchSubsys("qpp::random_circuit_depth()",
                                                  "with_respect_to_gate");
        }
        // check with_respect_to_gate is present in the gate sets
        if (is_one_qubit_gate) {
            bool found =
                std::find(one_qudit_gate_set.value().begin(),
                          one_qudit_gate_set.value().end(),
                          wrt_gate) != one_qudit_gate_set.value().end();
            if (!found) {
                throw exception::NotFound(
                    "qpp::random_circuit_depth()",
                    "with_respect_to_gate/one_qudit_gate_set");
            }
        }
        if (is_two_qubit_gate) {
            bool found =
                std::find(two_qudit_gate_set.value().begin(),
                          two_qudit_gate_set.value().end(),
                          wrt_gate) != two_qudit_gate_set.value().end();
            if (!found) {
                throw exception::NotFound(
                    "qpp::random_circuit_depth()",
                    "with_respect_to_gate/two_qudit_gate_set");
            }
        }
    } // end if (with_respect_to_gate.has_value())
    // END EXCEPTION CHECKS

    realT p_one = p_two.has_value() ? 1 - p_two.value() : 1;
    QCircuit qc{nq, 0, d};
    idx current_depth = 0;
    while (current_depth < gate_depth) {
        bool is_one_qudit_gate = bernoulli(p_one);
        if (is_one_qudit_gate) {
            idx q = randidx(0, nq - 1);
            idx gate = randidx(0, one_qudit_gate_set.value().size() - 1);
            if (!one_qudit_gate_names.has_value()) {
                qc.gate(one_qudit_gate_set.value()[gate], q);
            } else {
                qc.gate(one_qudit_gate_set.value()[gate], q,
                        one_qudit_gate_names.value()[gate]);
            }
        } else {
            idx ctrl = randidx(0, nq - 1);
            idx target = randidx(0, nq - 1);
            while (ctrl == target) {
                target = randidx(0, nq - 1);
            }
            idx gate = randidx(0, two_qudit_gate_set.value().size() - 1);
            if (!two_qudit_gate_names.has_value()) {
                qc.gate(two_qudit_gate_set.value()[gate], ctrl, target);
            } else {
                qc.gate(two_qudit_gate_set.value()[gate], ctrl, target,
                        two_qudit_gate_names.value()[gate]);
            }
        }
        current_depth = with_respect_to_gate.has_value()
                            ? qc.get_gate_depth(with_respect_to_gate.value())
                            : qc.get_gate_depth();
    }

    return qc;
}

/**
 * \brief Constructs a quantum phase estimation circuit with a \n bits of
 * precision
 * \note The number of classical bits will be equal with the number of
 * ancillas
 *
 * \param U Unitary matrix
 * \param n Number of ancilla qubits (translates to number of bits of
 * precision)
 * \param omit_measurements Optional (true by default); if true, omits
 * measurements of the ancilla qubits
 * \param d Subsystem dimensions (default is qubit, i.e., \a d =2)
 * \param name Optional ("qpe" by default) circuit name
 * \return Quantum phase estimation circuit with \a n bits of precision
 */
inline QCircuit qpe_circuit(cmat U, qpp::idx n, bool omit_measurements = true,
                            idx d = 2,
                            std::optional<std::string> name = "qpe") {
    // EXCEPTION CHECKS
    // check square matrix for the gate
    if (!internal::check_square_mat(U)) {
        throw exception::MatrixNotSquare("qpp::qpe_circuit()", "U");
    }

    auto D = static_cast<idx>(U.rows());
    idx m = internal::get_num_subsys(D, d);
    if (static_cast<idx>(std::pow(d, m)) != D) {
        throw exception::MatrixMismatchSubsys("qpp::qpe_circuit()", "U");
    }
    // END EXCEPTION CHECK

    QCircuit qc{static_cast<idx>(n + m), n, d, name};
    std::vector<idx> counting_qubits(n);
    std::iota(counting_qubits.begin(), counting_qubits.end(), 0);
    std::vector<idx> ancilla(m);
    std::iota(ancilla.begin(), ancilla.end(), n);

    if (d == 2) {
        qc.gate_fan(Gates::get_instance().H, counting_qubits);
    } else {
        qc.gate_fan(Gates::get_instance().Fd(d), counting_qubits, "Fd");
    }
    for (idx i = n; i-- > 0;) {
        qc.CTRL(U, i, ancilla);
        U = powm(U, 2);
    }
    qc.TFQ(counting_qubits); // inverse Fourier transform
    // measure many qubits at once, store starting with the 0
    // classical dit
    if (!omit_measurements) {
        qc.measure(counting_qubits);
    }

    return qc;
}

namespace internal {
/**
 * \brief True if the qpp::internal::QCircuitMeasurementStep is a projective
 * measurement step (including post-selection), false otherwise
 *
 * \param measurement_step Instance of qpp::internal::QCircuitMeasurementStep
 * \return True if the qpp::internal::QCircuitMeasurementStep is a projective
 * measurement step (including post-selection), false otherwise
 */
inline bool is_projective_measurement(
    const internal::QCircuitMeasurementStep& measurement_step) {
    switch (measurement_step.measurement_type_) {
        case internal::QCircuitMeasurementStep::Type::MEASURE:
        case internal::QCircuitMeasurementStep::Type::MEASURE_ND:
        case internal::QCircuitMeasurementStep::Type::MEASURE_MANY:
        case internal::QCircuitMeasurementStep::Type::MEASURE_MANY_ND:
        case internal::QCircuitMeasurementStep::Type::POST_SELECT:
        case internal::QCircuitMeasurementStep::Type::POST_SELECT_ND:
        case internal::QCircuitMeasurementStep::Type::POST_SELECT_MANY:
        case internal::QCircuitMeasurementStep::Type::POST_SELECT_MANY_ND:
            return true;
        default:
            return false;
    }
}

/**
 * \brief True if the quantum circuit step is a projective measurement step
 * (including post-selection), false otherwise
 *
 * \param elem Quantum circuit step
 * \return True if the quantum circuit step is a projective measurement
 * step (including post-selection), false otherwise
 */
inline bool
is_projective_measurement(const QCircuit::iterator::value_type& elem) {
    auto step = elem.get_step();
    if (std::holds_alternative<internal::QCircuitMeasurementStep>(step)) {
        auto measurement_step =
            std::get<internal::QCircuitMeasurementStep>(step);
        return is_projective_measurement(measurement_step);
    }

    return false;
}

/**
 * \brief True if the quantum circuit iterator points to a projective
 * measurement step (including post-selection), false otherwise
 *
 * \param it Quantum circuit iterator
 * \return True if the quantum circuit iterator points to a projective
 * measurement step (including post-selection), false otherwise
 */
inline bool is_projective_measurement(QCircuit::iterator it) {
    return is_projective_measurement(*it);
}

/**
 * \brief True if the qpp::internal::QCircuitMeasurementStep is a measurement
 * step (projective or not, including post-selection), false otherwise
 *
 * \param measurement_step Instance of qpp::internal::QCircuitMeasurementStep
 * \return True if the qpp::internal::QCircuitMeasurementStep is a measurement
 * step (projective or not, including post-selection), false otherwise
 */
inline bool
is_measurement(const internal::QCircuitMeasurementStep& measurement_step) {
    switch (measurement_step.measurement_type_) {
        case internal::QCircuitMeasurementStep::Type::MEASURE_V:
        case internal::QCircuitMeasurementStep::Type::MEASURE_V_ND:
        case internal::QCircuitMeasurementStep::Type::MEASURE_V_JOINT:
        case internal::QCircuitMeasurementStep::Type::MEASURE_V_JOINT_ND:
        case internal::QCircuitMeasurementStep::Type::POST_SELECT_V:
        case internal::QCircuitMeasurementStep::Type::POST_SELECT_V_ND:
        case internal::QCircuitMeasurementStep::Type::POST_SELECT_V_JOINT:
        case internal::QCircuitMeasurementStep::Type::POST_SELECT_V_JOINT_ND:
            return true;
        default:
            return is_projective_measurement(measurement_step);
    }
}

/**
 * \brief True if the quantum circuit step is a measurement step (projective
 * or not, including post-selection), false otherwise
 *
 * \param elem Quantum circuit step
 * \return True if the quantum circuit step is a measurement step
 * (projective or not, including post-selection), false otherwise
 */
inline bool is_measurement(const QCircuit::iterator::value_type& elem) {
    auto step = elem.get_step();
    if (std::holds_alternative<internal::QCircuitMeasurementStep>(step)) {
        auto measurement_step =
            std::get<internal::QCircuitMeasurementStep>(step);
        return is_measurement(measurement_step);
    }

    return false;
}

/**
 * \brief True if the quantum circuit iterator points to a measurement step
 * (projective or not, including post-selection), false otherwise
 *
 * \param it Quantum circuit iterator
 * \return True if the quantum circuit iterator points to a measurement step
 * (projective or not, including post-selection), false otherwise
 */
inline bool is_measurement(QCircuit::iterator it) {
    return is_measurement(*it);
}

/**
 * \brief True if the qpp::internal::QCircuitMeasurementStep performs a
 * destructive measurement of any kind (including post-selection), false
 * otherwise
 *
 * \param measurement_step Instance of qpp::internal::QCircuitMeasurementStep
 * \return True if the qpp::internal::QCircuitMeasurementStep performs a
 * destructive measurement of any kind (including post-selection), false
 * otherwise
 */
inline bool is_destructive_measurement(
    const internal::QCircuitMeasurementStep& measurement_step) {
    switch (measurement_step.measurement_type_) {
        case internal::QCircuitMeasurementStep::Type::MEASURE_ND:
        case internal::QCircuitMeasurementStep::Type::MEASURE_MANY_ND:
        case internal::QCircuitMeasurementStep::Type::MEASURE_V_ND:
        case internal::QCircuitMeasurementStep::Type::MEASURE_V_JOINT_ND:
        case internal::QCircuitMeasurementStep::Type::POST_SELECT_ND:
        case internal::QCircuitMeasurementStep::Type::POST_SELECT_MANY_ND:
        case internal::QCircuitMeasurementStep::Type::POST_SELECT_V_ND:
        case internal::QCircuitMeasurementStep::Type::POST_SELECT_V_JOINT_ND:
            return false;
        default:
            return true;
    }
}

/**
 * \brief True if the quantum circuit step performs a destructive measurement
 * of any kind, false otherwise
 *
 * \param elem Quantum circuit step
 * \return True if the quantum circuit step performs a destructive measurement
 * of any kind, false otherwise
 */
inline bool
is_destructive_measurement(const QCircuit::iterator::value_type& elem) {
    auto step = elem.get_step();
    if (std::holds_alternative<internal::QCircuitMeasurementStep>(step)) {
        auto measurement_step =
            std::get<internal::QCircuitMeasurementStep>(step);
        return is_destructive_measurement(measurement_step);
    }
    return false;
}

/**
 * \brief True if the quantum circuit iterator points to a destructive
 * measurement step of any kind, false otherwise
 *
 * \param it Quantum circuit iterator
 * \return True if the quantum circuit iterator points to a destructive
 * measurement step of any kind, false otherwise
 */
inline bool is_destructive_measurement(QCircuit::iterator it) {
    return is_destructive_measurement(*it);
}

/**
 * \brief True if the qpp::internal::QCircuitMeasurementStep is a projective
 * post-selection step, false otherwise
 *
 * \param measurement_step Instance of qpp::internal::QCircuitMeasurementStep
 * \return True if the qpp::internal::QCircuitMeasurementStep is a projective
 * post-selection step, false otherwise
 */
inline bool is_projective_post_selection(
    const internal::QCircuitMeasurementStep& measurement_step) {
    switch (measurement_step.measurement_type_) {
        case internal::QCircuitMeasurementStep::Type::POST_SELECT:
        case internal::QCircuitMeasurementStep::Type::POST_SELECT_ND:
        case internal::QCircuitMeasurementStep::Type::POST_SELECT_MANY:
        case internal::QCircuitMeasurementStep::Type::POST_SELECT_MANY_ND:
            return true;
        default:
            return false;
    }
}

/**
 * \brief True if the quantum circuit step is a projective post-selection step,
 * false otherwise
 *
 * \param elem Quantum circuit step
 * \return True if the quantum circuit step is a projective post-selection
 * step, false otherwise
 */
inline bool
is_projective_post_selection(const QCircuit::iterator::value_type& elem) {
    auto step = elem.get_step();
    if (std::holds_alternative<internal::QCircuitMeasurementStep>(step)) {
        auto measurement_step =
            std::get<internal::QCircuitMeasurementStep>(step);
        return is_projective_post_selection(measurement_step);
    }

    return false;
}

/**
 * \brief True if the quantum circuit iterator points to a projective
 * post-selection step, false otherwise
 *
 * \param it Quantum circuit iterator
 * \return True if the quantum circuit iterator points to a projective
 * post-selection step, false otherwise
 */
inline bool is_projective_post_selection(QCircuit::iterator it) {
    return is_projective_post_selection(*it);
}

/**
 * \brief True if the qpp::internal::QCircuitMeasurementStep is a
 * post-selection step (projective or not), false otherwise
 *
 * \param measurement_step Instance of qpp::internal::QCircuitMeasurementStep
 * \return True if the qpp::internal::QCircuitMeasurementStep is a
 * post-selection step (projective or not), false otherwise
 */
inline bool
is_post_selection(const internal::QCircuitMeasurementStep& measurement_step) {
    switch (measurement_step.measurement_type_) {
        case internal::QCircuitMeasurementStep::Type::POST_SELECT_V:
        case internal::QCircuitMeasurementStep::Type::POST_SELECT_V_ND:
        case internal::QCircuitMeasurementStep::Type::POST_SELECT_V_JOINT:
        case internal::QCircuitMeasurementStep::Type::POST_SELECT_V_JOINT_ND:
            return true;
        default:
            return is_projective_post_selection(measurement_step);
    }
}

/**
 * \brief True if the quantum circuit step is a post-selection step (projective
 * or not), false otherwise
 *
 * \param elem Quantum circuit step
 * \return True if the quantum circuit step is a post-selection step
 * (projective or not), false otherwise
 */
inline bool is_post_selection(const QCircuit::iterator::value_type& elem) {
    auto step = elem.get_step();
    if (std::holds_alternative<internal::QCircuitMeasurementStep>(step)) {
        auto measurement_step =
            std::get<internal::QCircuitMeasurementStep>(step);
        return is_post_selection(measurement_step);
    }

    return false;
}

/**
 * \brief True if the quantum circuit iterator points to a post-selection step
 * (projective or not), false otherwise
 *
 * \param it Quantum circuit iterator
 * \return True if the quantum circuit iterator points to a post-selection step
 * (projective or not), false otherwise
 */
inline bool is_post_selection(QCircuit::iterator it) {
    return is_post_selection(*it);
}

/**
 * \brief True if the qpp::internal::QCircuitMeasurementStep is a discard step,
 * false otherwise
 *
 * \param measurement_step Instance of qpp::internal::QCircuitMeasurementStep
 * \return True if the qpp::internal::QCircuitMeasurementStep is a discard
 * step, false otherwise
 */
inline bool
is_discard(const internal::QCircuitMeasurementStep& measurement_step) {
    switch (measurement_step.measurement_type_) {
        case internal::QCircuitMeasurementStep::Type::DISCARD:
        case internal::QCircuitMeasurementStep::Type::DISCARD_MANY:
            return true;
        default:
            return false;
    }
}

/**
 * \brief True if the quantum circuit step is a discard step, false
 * otherwise
 *
 * \param elem Quantum circuit step
 * \return True if the quantum circuit step is a discard step, false
 * otherwise
 */
inline bool is_discard(const QCircuit::iterator::value_type& elem) {
    auto step = elem.get_step();
    if (std::holds_alternative<internal::QCircuitMeasurementStep>(step)) {
        auto measurement_step =
            std::get<internal::QCircuitMeasurementStep>(step);
        return is_discard(measurement_step);
    }

    return false;
}

/**
 * \brief True if the quantum circuit iterator points to a discard step,
 * false otherwise
 *
 * \param it Quantum circuit iterator
 * \return True if the quantum circuit iterator points to a discard step,
 * false otherwise
 */
inline bool is_discard(QCircuit::iterator it) { return is_discard(*it); }

/**
 * \brief True if the qpp::internal::QCircuitMeasurementStep is a reset step,
 * false otherwise
 *
 * \param measurement_step Instance of qpp::internal::QCircuitMeasurementStep
 * \return True if the qpp::internal::QCircuitMeasurementStep is a reset step,
 * false otherwise
 */
inline bool
is_reset(const internal::QCircuitMeasurementStep& measurement_step) {
    switch (measurement_step.measurement_type_) {
        case internal::QCircuitMeasurementStep::Type::RESET:
        case internal::QCircuitMeasurementStep::Type::RESET_MANY:
            return true;
        default:
            return false;
    }
}

/**
 * \brief True if the quantum circuit step is a reset step,
 * false otherwise
 *
 * \param elem Quantum circuit step
 * \return True if the quantum circuit step is a reset step,
 * false otherwise
 */
inline bool is_reset(const QCircuit::iterator::value_type& elem) {
    auto step = elem.get_step();
    if (std::holds_alternative<internal::QCircuitMeasurementStep>(step)) {
        auto measurement_step =
            std::get<internal::QCircuitMeasurementStep>(step);
        return is_reset(measurement_step);
    }

    return false;
}

/**
 * \brief True if the quantum circuit iterator points to a reset step, false
 * otherwise
 *
 * \param it Quantum circuit iterator
 * \return True if the quantum circuit iterator points to a reset step,
 * false otherwise
 */
inline bool is_reset(QCircuit::iterator it) { return is_reset(*it); }

/** \brief True if the qpp::internal::QCircuitGateStep is a classical CTRL
 * step, false otherwise
 *
 * \param gate_step Instance of qpp::internal::QCircuitGateStep
 * \return True if the qpp::internal::QCircuitGateStep is a classical CTRL
 * step, false otherwise
 */
inline bool is_cCTRL(const internal::QCircuitGateStep& gate_step) {
    return QCircuit::is_cCTRL(gate_step);
}

/** \brief True if the quantum circuit step is a classical CTRL step, false
 * otherwise
 *
 * \param elem Quantum circuit step
 * \return True if the quantum circuit step is a classical CTRL step, false
 * otherwise
 */
inline bool is_cCTRL(const QCircuit::iterator::value_type& elem) {
    auto step = elem.get_step();
    if (std::holds_alternative<internal::QCircuitGateStep>(step)) {
        auto gate_step = std::get<internal::QCircuitGateStep>(step);
        return is_cCTRL(gate_step);
    }

    return false;
}

/**
 * \brief True if the quantum circuit iterator points to a classical CTRL
 * step, false otherwise
 *
 * \param it Quantum circuit iterator
 * \return True if the quantum circuit iterator points to a classical CTRL
 * step, false otherwise
 */
inline bool is_cCTRL(QCircuit::iterator it) { return is_cCTRL(*it); }

/**
 * \brief Extracts ctrl, target, ps_vals, and c_reg vectors (in this order)
 * from a quantum circuit step, as a tuple
 *
 * \param elem Quantum circuit step
 * \return Tuple with vectors representing (in this order) ctrl, target, and
 * c_reg, respectively
 */
inline std::tuple<std::vector<idx>, std::vector<idx>, std::vector<idx>,
                  std::vector<idx>>
extract_ctrl_target_ps_vals_c_reg(const QCircuit::iterator::value_type& elem) {
    std::vector<idx> ctrl, target, ps_vals, c_reg;

    // measurement
    if (is_measurement(elem)) {
        auto current_measurement_step =
            std::get<internal::QCircuitMeasurementStep>(elem.get_step());

        ctrl = {};
        target = current_measurement_step.target_;
        if (current_measurement_step.ps_vals_.has_value()) {
            ps_vals = current_measurement_step.ps_vals_.value();
        }

        switch (current_measurement_step.measurement_type_) {
            // measure()/post-select()
            case internal::QCircuitMeasurementStep::Type::MEASURE:
            case internal::QCircuitMeasurementStep::Type::MEASURE_ND:
            case internal::QCircuitMeasurementStep::Type::MEASURE_MANY:
            case internal::QCircuitMeasurementStep::Type::MEASURE_MANY_ND:
            case internal::QCircuitMeasurementStep::Type::POST_SELECT:
            case internal::QCircuitMeasurementStep::Type::POST_SELECT_ND:
            case internal::QCircuitMeasurementStep::Type::POST_SELECT_MANY:
            case internal::QCircuitMeasurementStep::Type::POST_SELECT_MANY_ND:
                c_reg.resize(target.size());
                std::iota(c_reg.begin(), c_reg.end(),
                          current_measurement_step.c_reg_);
                break;
            // measureV()/post-selectV()
            case internal::QCircuitMeasurementStep::Type::MEASURE_V:
            case internal::QCircuitMeasurementStep::Type::MEASURE_V_ND:
            case internal::QCircuitMeasurementStep::Type::MEASURE_V_JOINT:
            case internal::QCircuitMeasurementStep::Type::MEASURE_V_JOINT_ND:
            case internal::QCircuitMeasurementStep::Type::POST_SELECT_V:
            case internal::QCircuitMeasurementStep::Type::POST_SELECT_V_ND:
            case internal::QCircuitMeasurementStep::Type::POST_SELECT_V_JOINT:
            case internal::QCircuitMeasurementStep::Type::
                POST_SELECT_V_JOINT_ND:
                c_reg.resize(1);
                c_reg[0] = current_measurement_step.c_reg_;
                break;
            // discard/reset/none
            default:
                c_reg = {};
        }
    }
    // cCTRL
    else if (is_cCTRL(elem)) {
        auto current_gate_step =
            std::get<internal::QCircuitGateStep>(elem.get_step());

        ctrl = {};
        target = current_gate_step.target_;
        c_reg = current_gate_step.ctrl_.value();

    }
    // otherwise
    else {
        auto step = elem.get_step();
        if (std::holds_alternative<internal::QCircuitGateStep>(step)) {
            auto current_gate_step = std::get<internal::QCircuitGateStep>(step);

            ctrl = current_gate_step.ctrl_.value_or(std::vector<idx>{});
            target = current_gate_step.target_;
            c_reg = {};
        }
    }

    return {ctrl, target, ps_vals, c_reg};
}

/**
 * \brief Extracts ctrl, target, and c_reg vectors (in this order) from a
 * quantum circuit iterator, as a tuple
 *
 * \param it Quantum circuit iterator
 * \return Tuple with vectors representing (in this order) ctrl, target, and
 * c_reg, respectively
 */
inline std::tuple<std::vector<idx>, std::vector<idx>, std::vector<idx>,
                  std::vector<idx>>
extract_ctrl_target_ps_vals_c_reg(QCircuit::iterator it) {
    return extract_ctrl_target_ps_vals_c_reg(*it);
}

/**
 * \brief True if two quantum circuit steps (assumed added via the
 * qpp::QCircuit API) can be swapped, false otherwise
 *
 * \param elem1 Quantum circuit step
 * \param elem2 Quantum circuit step
 * \return True if two quantum circuit steps can be swapped,
 * false otherwise
 */
inline bool can_swap(const QCircuit::iterator::value_type& elem1,
                     const QCircuit::iterator::value_type& elem2) {

    auto [ctrl1, target1, ps_vals1, c_reg1] =
        extract_ctrl_target_ps_vals_c_reg(elem1);
    auto [ctrl2, target2, ps_vals2, c_reg2] =
        extract_ctrl_target_ps_vals_c_reg(elem2);

    // if ctrls overlap, return false
    if (std::find_first_of(ctrl1.begin(), ctrl1.end(), ctrl2.begin(),
                           ctrl2.end()) != ctrl1.end()) {
        return false;
    }

    // if targets overlap, return false
    if (std::find_first_of(target1.begin(), target1.end(), target2.begin(),
                           target2.end()) != target1.end()) {
        return false;
    }

    // if c_reg overlap, return false
    if (std::find_first_of(c_reg1.begin(), c_reg1.end(), c_reg2.begin(),
                           c_reg2.end()) != c_reg1.end()) {
        return false;
    }

    return true;
}

/**
 * \brief True if two quantum circuit steps (assumed added via the
 * qpp::QCircuit API) can be swapped, false otherwise
 *
 * \param it1 Quantum circuit iterator pointing to a quantum circuit step
 * \param it2 Quantum circuit iterator pointing to a quantum circuit step
 * \return True if two quantum circuit steps can be swapped, false otherwise
 */
inline bool can_swap(QCircuit::iterator it1, QCircuit::iterator it2) {
    return can_swap(*it1, *it2);
}

/**
 * \brief Converts a quantum (sub)-circuit description to a vector of
 * quantum circuit iterators
 *
 * \param start Quantum circuit iterator pointing to the first element
 * \param finish Quantum circuit iterator pointing to the last element (not
 * included)
 * \return Vector of quantum circuit iterators
 */
inline std::vector<QCircuit::iterator>
circuit_as_iterators(QCircuit::iterator start, QCircuit::iterator finish) {
    std::vector<QCircuit::iterator> steps; // circuit steps (as iterators)
    for (auto it = start; it != finish; ++it) {
        steps.emplace_back(it);
    }

    return steps;
}

/**
 * \brief Converts a quantum (sub)-circuit description to a vector of
 * quantum circuit iterators
 *
 * \param qc Quantum circuit description
 * \return Vector of quantum circuit iterators
 */
inline std::vector<QCircuit::iterator>
circuit_as_iterators(const QCircuit& qc) {
    return circuit_as_iterators(qc.begin(), qc.end());
}

/**
 * \brief Puts a quantum (sub)-circuit description in the canonical form,
 * i.e., starting with the first measurement step in the circuit range
 * [start, finish), pushes all cCTRLs and measurements to the end of the
 * circuit, so, if possible, the circuit will be of the form
 * [Gates, cCTRLs, Measurements]
 *
 * \note This function does not interchange measurements, i.e., the re-ordering
 * is stable
 *
 * \param start Quantum circuit iterator pointing to the first element
 * \param finish Quantum circuit iterator pointing to the last element (not
 * included)
 * \return Quantum circuit canonical form represented as a vector of quantum
 * circuit iterators
 */
inline std::vector<QCircuit::iterator>
canonical_form(QCircuit::iterator start, QCircuit::iterator finish) {

    auto first_measurement_it = std::find_if(
        start, finish, [](auto&& elem) { return is_measurement(elem); });

    std::vector<QCircuit::iterator> steps_before_measurement =
        circuit_as_iterators(start, first_measurement_it);

    std::vector<QCircuit::iterator> steps =
        circuit_as_iterators(first_measurement_it, finish);

    // push all cCTRL steps to the end
    for (auto&& rit = steps.rbegin(); rit != steps.rend(); ++rit) {
        if (internal::is_cCTRL(*rit)) {
            // try to push it to the end
            for (auto&& after_it = rit.base(); after_it != steps.end();
                 ++after_it) {
                auto cur_it = std::prev(after_it);
                if (internal::can_swap(*cur_it, *after_it) &&
                    !internal::is_cCTRL(*after_it)) {
                    std::swap(*cur_it, *after_it);
                } else {
                    break;
                }
            }
        }
    }

    // next push all measurement steps to the end
    for (auto&& rit = steps.rbegin(); rit != steps.rend(); ++rit) {
        if (internal::is_measurement(*rit)) {
            // try to push it to the end
            for (auto&& after_it = rit.base(); after_it != steps.end();
                 ++after_it) {
                auto cur_it = std::prev(after_it);
                if (internal::can_swap(*cur_it, *after_it) &&
                    !internal::is_measurement(*after_it)) {
                    std::swap(*cur_it, *after_it);
                } else {
                    break;
                }
            }
        }
    }

    steps_before_measurement.insert(steps_before_measurement.end(),
                                    steps.begin(), steps.end());

    return steps_before_measurement;
}

/**
 * \brief Puts a quantum (sub)-circuit description in the canonical form,
 * i.e., starting with the first measurement step in the circuit, pushes all
 * cCTRLs and measurements to the end of the circuit, so, if possible, the
 * circuit will be of the form [Gates, cCTRLs, Measurements]
 *
 * \note This function does not interchange measurements, i.e., the re-ordering
 * is stable
 *
 * \param qc Quantum circuit description
 * \return Quantum circuit canonical form represented as a vector of quantum
 * circuit iterators
 */
inline std::vector<QCircuit::iterator> canonical_form(const QCircuit& qc) {
    return canonical_form(qc.begin(), qc.end());
}

} /* namespace internal */
} /* namespace qpp */

#endif /* QPP_CLASSES_QCIRCUIT_HPP_ */
