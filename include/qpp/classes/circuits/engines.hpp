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
 * \file classes/circuits/engines.hpp
 * \brief Qudit quantum engines
 */

#ifndef QPP_CLASSES_CIRCUITS_ENGINES_HPP_
#define QPP_CLASSES_CIRCUITS_ENGINES_HPP_

#include <algorithm>
#include <cmath>
#include <iterator>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <optional>
#include <ostream>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <variant>
#include <vector>

#include "qpp/input_output.hpp"
#include "qpp/instruments.hpp"
#include "qpp/operations.hpp"
#include "qpp/types.hpp"

#include "qpp/classes/circuits/circuits.hpp"
#include "qpp/classes/idisplay.hpp"
#include "qpp/classes/states.hpp"

namespace qpp {
/**
 * \class qpp::QEngine
 * \brief Quantum circuit engine, executes qpp::QCircuit
 * \see qpp::QCircuit
 *
 * \note When interfacing with OpenQASM, the classical dits/registers are
 * evaluated in little-endian order, with  the least significant bit being
 * stored first. For example, [1,0,0] is interpreted as 1 (and not 4).
 * See https://github.com/softwareQinc/qpp/issues/75 for more details.
 */
class QEngine : public IDisplay, public IJSON {
  public:
    /**
     * \brief Sampling/measurement statistics
     */
    class Statistics : public IDisplay, public IJSON {
        /**
         * \brief Measurement/sampling statistics
         */
        using stats_t_ = std::map<std::vector<idx>, idx>;
        mutable stats_t_ stats_data_{}; ///< statistics data

      public:
        /**
         * \brief Default constructor
         */
        Statistics() = default;

        /**
         * \brief Constructor
         *
         * \param stats Instance of qpp::QEngine::Statistics
         */
        explicit Statistics(stats_t_ stats) : stats_data_{std::move(stats)} {}

        /**
         * \brief Number of samples
         *
         * \return Number of samples
         */
        idx get_num_reps() const {
            idx result = 0;
            for (auto&& [_, val] : stats_data_) {
                result += val;
            }
            return result;
        }

        /**
         * \brief Number of distinct outcomes
         *
         * \return Number of distinct outcomes
         */
        idx get_num_outcomes() const { return stats_data_.size(); }

        /**
         * \brief Raw data structure representing the statistics
         *
         * \return Raw data structure representing the statistics
         */
        stats_t_& data() const& { return stats_data_; }

        /**
         * \brief qpp::IJSON::to_JSON() override
         *
         * Displays the statistics in JSON format
         *
         * \param enclosed_in_curly_brackets If true, encloses the result in
         * curly brackets
         * \return String containing the JSON representation of the statistics
         */
        std::string
        to_JSON(bool enclosed_in_curly_brackets = true) const override {
            std::string result;

            if (enclosed_in_curly_brackets) {
                result += "{";
            }

            std::ostringstream ss;
            ss << "\"num_reps\": " << get_num_reps() << ", ";
            ss << "\"num_outcomes\": " << get_num_outcomes() << ", ";

            bool is_first = true;
            std::string sep{};
            ss << "\"outcomes\": ";
            ss << "{";
            for (auto&& [key, val] : stats_data_) {
                ss << sep << "\""
                   << disp(key, IOManipContainerOpts{}.set_sep(", "))
                   << "\": " << val;
                if (is_first) {
                    is_first = false;
                    sep = ", ";
                }
            }
            ss << "}";
            result += ss.str();

            if (enclosed_in_curly_brackets) {
                result += "}";
            }

            return result;
        }

      private:
        /**
         * \brief qpp::IDisplay::display() override
         *
         * Writes to the output stream a textual representation of the
         * statistics
         *
         * \param os Output stream passed by reference
         * \return Reference to the output stream
         */
        std::ostream& display(std::ostream& os) const override {
            os << "[Statistics]\n";
            os << "\tnum_reps: " << get_num_reps() << '\n';
            os << "\tnum_outcomes: " << get_num_outcomes() << '\n';
            bool is_first = true;
            std::string sep{};
            for (auto&& [key, val] : stats_data_) {
                os << sep << '\t'
                   << disp(key, IOManipContainerOpts{}.set_sep(" ")) << ": "
                   << val;
                if (is_first) {
                    is_first = false;
                    sep = '\n';
                }
            }

            return os;
        };
    }; /* class QEngine::Statistics */

  protected:
    const QCircuit*
        qc_ptr_; ///< pointer to constant quantum circuit description

    /**
     * \class qpp::QEngine::state_
     * \brief Current state of the engine
     */
    struct state_ {
        const QCircuit* qc_ptr_;     ///< non-owning pointer to the parent
                                     ///< const quantum circuit description
        ket psi_{};                  ///< state vector
        std::vector<realT> probs_{}; ///< measurement probabilities
        std::vector<idx> dits_{};    ///< classical dits (where measurement
                                     ///< results are usually stored)
        std::vector<idx> subsys_{};  ///< keeps track of the measured
                                     ///< subsystems, re-label them after
                                     ///< measurements

        /**
         * \brief Constructor
         *
         * \param qc_ptr Non-owning pointer to the parent const quantum circuit
         * description
         */
        explicit state_(const QCircuit* qc_ptr) : qc_ptr_{qc_ptr} {
            // EXCEPTION CHECKS

            if (qc_ptr->get_nq() == 0) {
                throw exception::ZeroSize("qpp::QEngine::state_::reset()",
                                          "nq");
            }
            // END EXCEPTION CHECKS
            reset();
        }

        // silence -Weffc++ class has pointer data members
        /**
         * \brief Default copy constructor
         */
        state_(const state_&) = default;

        // silence -Weffc++ class has pointer data members
        /**
         * \brief Default copy assignment operator
         *
         * \return Reference to the current instance
         */
        state_& operator=(const state_&) = default;

        /**
         * \brief Resets the engine state
         *
         * \param psi Optional engine's initial quantum state
         */
        void reset(std::optional<ket> psi = std::nullopt) {
            if (psi.has_value()) {
                idx D = static_cast<idx>(std::llround(
                    std::pow(qc_ptr_->get_d(), qc_ptr_->get_nq())));
                if (static_cast<idx>(psi.value().rows()) != D) {
                    if (static_cast<idx>(psi.value().rows()) != D) {
                        throw exception::DimsNotEqual(
                            "qpp::QEngine::state_::reset()", "psi");
                    }
                }
                psi_ = psi.value();
            } else {
                psi_ = States::get_no_thread_local_instance().zero(
                    qc_ptr_->get_nq(), qc_ptr_->get_d());
            }
            probs_ = std::vector<realT>(qc_ptr_->get_nc(), 0);
            dits_ = std::vector<idx>(qc_ptr_->get_nc(), 0);
            subsys_ = std::vector<idx>(qc_ptr_->get_nq(), 0);
            std::iota(subsys_.begin(), subsys_.end(), 0);
        }
    } st_;               ///< current state of the engine
    Statistics stats_{}; ///< measurement statistics for multiple runs

    /**
     * \brief Marks qudit \a i as measured then re-label accordingly the
     * remaining non-measured qudits
     * \param i Qudit index
     */
    void set_measured_(idx i) {
        // EXCEPTION CHECKS

        if (was_measured(i)) {
            throw exception::QuditAlreadyMeasured(
                "qpp::QEngine::set_measured_()", "i");
        }
        // END EXCEPTION CHECKS
        st_.subsys_[i] =
            std::numeric_limits<idx>::max(); // set qudit i to measured state
        for (idx m = i; m < qc_ptr_->get_nq(); ++m) {
            if (!was_measured(m)) {
                --st_.subsys_[m];
            }
        }
    }

    // giving a vector of non-measured qudits, get their relative position wrt
    // the measured qudits
    /**
     * \brief Giving a vector \a v of non-measured qudits, gets their relative
     * position with respect to the measured qudits
     *
     * \param v Vector of non-measured qudit indexes
     * \return Vector of qudit indexes
     */
    std::vector<idx> get_relative_pos_(std::vector<idx> v) {
        idx vsize = v.size();
        for (idx i = 0; i < vsize; ++i) {
            // EXCEPTION CHECKS

            if (was_measured(v[i])) {
                throw exception::QuditAlreadyMeasured(
                    "qpp::QEngine::get_relative_pos_()", "v[i]");
            }
            // END EXCEPTION CHECKS
            v[i] = st_.subsys_[v[i]];
        }
        return v;
    }

  private:
    bool can_sample; ///< can sample when executing with multiple repetitions

  public:
    /**
     * \brief Constructs a quantum engine out of a quantum circuit description
     *
     * \note The quantum circuit description must be an lvalue
     * \see qpp::QEngine(QCircuit&&)
     *
     * \note The initial underlying quantum state is set to
     * \f$|0\rangle^{\otimes n}\f$
     *
     * \param qc Quantum circuit description
     */
    explicit QEngine(const QCircuit& qc)
        : qc_ptr_{std::addressof(qc)}, st_{qc_ptr_}, stats_{},
          can_sample{false} {}

    // silence -Weffc++ class has pointer data members
    /**
     * \brief Default copy constructor
     */
    QEngine(const QEngine&) = default;

    // silence -Weffc++ class has pointer data members
    /**
     * \brief Default copy assignment operator
     *
     * \return Reference to the current instance
     */
    QEngine& operator=(const QEngine&) = default;

    /**
     * \brief Disables rvalue QCircuit
     */
    QEngine(QCircuit&&) = delete;

    /**
     * \brief Default virtual destructor
     */
    ~QEngine() override = default;

    // getters
    /**
     * \brief Underlying quantum state
     *
     * \return Underlying quantum state
     */
    ket get_psi() const { return st_.psi_; }

    /**
     * \brief Vector with the values of the underlying classical dits
     * \see qpp::QEngine::set_dits()
     *
     * \note When interfacing with OpenQASM, the classical dits/registers
     * are evaluated in little-endian order, with  the least significant bit
     * being stored first. For example, [1,0,0] is interpreted as 1 (and not
     * 4).
     *
     * \return Vector of underlying classical dits
     */
    std::vector<idx> get_dits() const { return st_.dits_; }

    /**
     * \brief Value of the classical dit at position \a i
     * \see qpp::QEngine::set_dit()
     *
     * \note When interfacing with OpenQASM, the classical dits/registers
     * are evaluated in little-endian order, with  the least significant bit
     * being stored first. For example, [1,0,0] is interpreted as 1 (and not
     * 4).
     *
     * \param i Classical dit index
     * \return Value of the classical dit at position \a i
     */
    idx get_dit(idx i) const {
        // EXCEPTION CHECKS

        if (i >= qc_ptr_->get_nc()) {
            throw exception::OutOfRange("qpp::QEngine::get_dit()", "i");
        }
        // END EXCEPTION CHECKS

        return st_.dits_[i];
    }

    /**
     * \brief Vector of underlying measurement outcome probabilities
     *
     * \note Those should be interpreted as conditional probabilities based
     * on the temporal order of the measurements, i.e., if we measure qubit
     * 0, then measure qubit 1, and finally qubit 2, the resulting vector of
     * outcome probabilities probs[2] should be interpreted as the
     * conditional probability of qubit 2 having the outcome it had given
     * that qubit 1 and qubit 0 had their given outcomes, respectively. As
     * an example, if we measure the qubit 0 followed by the qubit 1 of a
     * maximally entangled state \f$(|00\rangle + |11\rangle)/\sqrt{2}\f$,
     * then the vector of outcome probabilities will be [0.5, 1].
     *
     * \note The probability vector has the same length as the vector of
     * classical dits. If the measurement result is stored at the index
     * \a c_reg, then the outcome probability is automatically stored at
     * the same index \a c_reg in the probability vector.
     *
     * \return Vector of underlying measurement outcome probabilities
     */
    std::vector<realT> get_probs() const { return st_.probs_; }

    /**
     * \brief Check whether qudit \a i was already measured (destructively)
     * at the current engine state
     *
     * \param i Qudit index
     * \return True if qudit \a i was already measured, false otherwise
     */
    bool was_measured(idx i) const {
        return st_.subsys_[i] == std::numeric_limits<idx>::max();
    }

    /**
     * \brief Vector of already measured (destructively) qudit indexes at
     * the current engine state
     *
     * \return Vector of already measured qudit (destructively) indexes at
     * the current engine state
     */
    std::vector<idx> get_measured() const {
        std::vector<idx> result;
        for (idx i = 0; i < qc_ptr_->get_nq(); ++i) {
            if (was_measured(i)) {
                result.emplace_back(i);
            }
        }

        return result;
    }

    /**
     * \brief Vector of non-measured qudit indexes
     *
     * \return Vector of non-measured qudit indexes
     */
    std::vector<idx> get_non_measured() const {
        std::vector<idx> result;
        for (idx i = 0; i < qc_ptr_->get_nq(); ++i) {
            if (!was_measured(i)) {
                result.emplace_back(i);
            }
        }

        return result;
    }

    /**
     * \brief Quantum circuit description, lvalue ref qualifier
     *
     * \return Const reference to the underlying quantum circuit description
     */
    const QCircuit& get_circuit() const& noexcept { return *qc_ptr_; }

    /**
     * \brief Quantum circuit description, rvalue ref qualifier
     *
     * \return Copy of the underlying quantum circuit description
     */
    QCircuit get_circuit() const&& noexcept { return *qc_ptr_; }

    /**
     * \brief Measurement statistics for multiple runs
     *
     * \return Hash table with collected measurement statistics for multiple
     * runs, with hash key being the string representation of the vector of
     * measurement results and value being the number of occurrences (of the
     * vector of measurement results), with the most significant bit located
     * at index 0 (i.e., top/left).
     */
    QEngine::Statistics get_stats() const { return stats_; }

    /**
     * \brief Determines if engines derived from \a qpp::QEngine are noisy or
     * not at runtime
     *
     * \return True if the engine is noisy, false otherwise
     */
    virtual bool is_noisy() const { return false; }
    // end getters

    // setters
    /**
     * \brief Sets the classical dit at position \a i
     * \see qpp::QEngine::get_dit()
     *
     * \note When interfacing with OpenQASM, the classical dits/registers
     * are evaluated in little-endian order, with  the least significant bit
     * being stored first. For example, [1,0,0] is interpreted as 1 (and not
     * 4).
     *
     * \param i Classical dit index
     * \param value Classical dit value
     * \return Reference to the current instance
     */
    QEngine& set_dit(idx i, idx value) {
        // EXCEPTION CHECKS

        if (i >= qc_ptr_->get_nc()) {
            throw exception::OutOfRange("qpp::QEngine::set_dit()", "i");
        }
        // END EXCEPTION CHECKS
        st_.dits_[i] = value;

        return *this;
    }

    /**
     * \brief Set the classical dits to \a dits
     * \see qpp::QEngine::get_dits()
     *
     * \note When interfacing with OpenQASM, the classical dits/registers
     * are evaluated in little-endian order, with  the least significant bit
     * being stored first. For example, [1,0,0] is interpreted as 1 (and not 4).
     *
     * \param dits Vector of classical dits, must have the same size as the
     * internal vector of classical dits returned by qpp::QEngine::get_dits()
     * \return Reference to the current instance
     */
    QEngine& set_dits(std::vector<idx> dits) {
        // EXCEPTION CHECKS

        if (dits.size() != st_.dits_.size()) {
            throw exception::SizeMismatch("qpp::QEngine::set_dits()", "dits");
        }
        // END EXCEPTION CHECKS
        st_.dits_ = std::move(dits);

        return *this;
    }

    /**
     * \brief Sets the underlying quantum state to \a psi
     *
     * \note The order is lexicographical with respect to the remaining
     * non-measured qudits
     *
     * \param psi State vector
     * \return Reference to the current instance
     */
    QEngine& set_psi(const ket& psi) {
        // EXCEPTION CHECKS

        idx n = get_non_measured().size();
        idx D = static_cast<idx>(std::llround(std::pow(qc_ptr_->get_d(), n)));
        if (static_cast<idx>(psi.rows()) != D) {
            throw exception::DimsNotEqual("qpp::QEngine::set_psi()", "psi");
        }
        // END EXCEPTION CHECKS

        st_.psi_ = psi;

        return *this;
    }
    // end setters

    /**
     * \brief Resets the collected measurement statistics hash table
     *
     * \return Reference to the current instance
     */
    QEngine& reset_stats() {
        stats_ = {};

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
    virtual QEngine& reset(bool reset_stats = true) {
        this->can_sample = false;
        st_.reset(st_.psi_);
        if (reset_stats) {
            this->reset_stats();
        }

        return *this;
    }

    /**
     * \brief Executes one step in the quantum circuit description
     *
     * \note Override only this QEngine::execute() member function in every
     * derived class to achieve the desired behaviour
     *
     * \param elem Step to be executed
     * \return Reference to the current instance
     */
    virtual QEngine& execute(const QCircuit::iterator::value_type& elem) {
        // EXCEPTION CHECKS

        // iterator must point to the same quantum circuit description
        if (elem.get_qc_ptr() != qc_ptr_) {
            throw exception::InvalidIterator(
                "qpp::QEngine::execute()",
                "Iterator does not point to the same circuit description");
        }
        // the rest of exceptions are caught by the iterator::operator*()
        // END EXCEPTION CHECKS

        auto h_tbl = qc_ptr_->get_cmat_hash_tbl();
        idx d = qc_ptr_->get_d();

        std::visit(
            overloaded{
                [&](const QCircuit::GateStep& gate_step) {
                    std::vector<idx> ctrl_rel_pos;
                    std::vector<idx> target_rel_pos =
                        get_relative_pos_(gate_step.target_);

                    // regular gate
                    switch (gate_step.gate_type_) {
                        case QCircuit::GateStep::Type::NONE:
                            break;
                        case QCircuit::GateStep::Type::SINGLE:
                        case QCircuit::GateStep::Type::TWO:
                        case QCircuit::GateStep::Type::THREE:
                        case QCircuit::GateStep::Type::JOINT:
                            st_.psi_ =
                                apply(st_.psi_, h_tbl[gate_step.gate_hash_],
                                      target_rel_pos, d);
                            break;
                        case QCircuit::GateStep::Type::FAN:
                            for (idx m = 0;
                                 m < static_cast<idx>(gate_step.target_.size());
                                 ++m) {
                                st_.psi_ =
                                    apply(st_.psi_, h_tbl[gate_step.gate_hash_],
                                          {target_rel_pos[m]}, d);
                            }
                            break;
                        default:
                            break;
                    }

                    // controlled gate
                    if (QCircuit::is_CTRL(gate_step)) {
                        ctrl_rel_pos =
                            get_relative_pos_(gate_step.ctrl_.value());
                        bool is_fan = (gate_step.gate_type_ ==
                                       QCircuit::GateStep::Type::CTRL_FAN);
                        st_.psi_ =
                            is_fan ? applyCTRL_fan(st_.psi_,
                                                   h_tbl[gate_step.gate_hash_],
                                                   ctrl_rel_pos, target_rel_pos,
                                                   d, gate_step.shift_)
                                   : applyCTRL(st_.psi_,
                                               h_tbl[gate_step.gate_hash_],
                                               ctrl_rel_pos, target_rel_pos, d,
                                               gate_step.shift_);
                    }

                    // classically-controlled gate
                    if (QCircuit::is_cCTRL(gate_step)) {
                        bool is_fan = (gate_step.gate_type_ ==
                                       QCircuit::GateStep::Type::cCTRL_FAN);
                        if (!st_.dits_.empty()) {
                            {
                                bool should_apply = true;
                                idx first_dit;
                                // we have a shift
                                if (gate_step.shift_.has_value()) {
                                    first_dit =
                                        (st_.dits_[(
                                             gate_step.ctrl_.value())[0]] +
                                         gate_step.shift_.value()[0]) %
                                        d;
                                    for (idx m = 1;
                                         m <
                                         static_cast<idx>(
                                             gate_step.ctrl_.value().size());
                                         ++m) {
                                        if ((st_.dits_[(
                                                 gate_step.ctrl_.value())[m]] +
                                             gate_step.shift_.value()[m]) %
                                                d !=
                                            first_dit) {
                                            should_apply = false;
                                            break;
                                        }
                                    }
                                }
                                // no shift
                                else {
                                    first_dit =
                                        st_.dits_[(gate_step.ctrl_.value())[0]];
                                    for (idx m = 1;
                                         m <
                                         static_cast<idx>(
                                             gate_step.ctrl_.value().size());
                                         ++m) {
                                        if (st_.dits_[(
                                                gate_step.ctrl_.value())[m]] !=
                                            first_dit) {
                                            should_apply = false;
                                            break;
                                        }
                                    }
                                }
                                if (should_apply) {
                                    cmat U = powm(h_tbl[gate_step.gate_hash_],
                                                  first_dit);
                                    if (is_fan) {
                                        for (idx qudit : target_rel_pos) {
                                            st_.psi_ =
                                                apply(st_.psi_, U, {qudit}, d);
                                        }
                                    } else {
                                        st_.psi_ = apply(
                                            st_.psi_,
                                            powm(h_tbl[gate_step.gate_hash_],
                                                 first_dit),
                                            target_rel_pos, d);
                                    }
                                }
                            }
                        }
                        // TODO check if this can happen (st_.dits_.empty())
                        else {
                            if (is_fan) {
                                for (idx qudit : target_rel_pos) {
                                    st_.psi_ = apply(
                                        st_.psi_, h_tbl[gate_step.gate_hash_],
                                        {qudit}, d);
                                }
                            } else {
                                st_.psi_ =
                                    apply(st_.psi_, h_tbl[gate_step.gate_hash_],
                                          target_rel_pos, d);
                            }
                        }
                    } // end if classically-controlled gate
                },
                [&](const QCircuit::MeasurementStep& measure_step) {
                    std::vector<idx> target_rel_pos =
                        get_relative_pos_(measure_step.target_);

                    idx mres = 0;
                    std::vector<idx> results;
                    std::vector<realT> probs;
                    std::vector<ket> states;

                    switch (measure_step.measurement_type_) {
                        case QCircuit::MeasurementStep::Type::NONE:
                            break;
                        case QCircuit::MeasurementStep::Type::MEASURE:
                            std::tie(results, probs, st_.psi_) =
                                measure_seq(st_.psi_, target_rel_pos, d);
                            st_.dits_[measure_step.c_reg_] = results[0];
                            st_.probs_[measure_step.c_reg_] = probs[0];
                            set_measured_(measure_step.target_[0]);
                            break;
                        case QCircuit::MeasurementStep::Type::MEASURE_MANY:
                            std::tie(results, probs, st_.psi_) =
                                measure_seq(st_.psi_, target_rel_pos, d);
                            std::copy(results.begin(), results.end(),
                                      std::next(st_.dits_.begin(),
                                                measure_step.c_reg_));
                            std::copy(probs.begin(), probs.end(),
                                      std::next(st_.probs_.begin(),
                                                measure_step.c_reg_));
                            for (idx target : measure_step.target_) {
                                set_measured_(target);
                            }
                            break;
                        case QCircuit::MeasurementStep::Type::MEASURE_V:
                            std::tie(mres, probs, states) = measure(
                                st_.psi_, h_tbl[measure_step.mats_hash_[0]],
                                target_rel_pos, d);
                            st_.psi_ = states[mres];
                            st_.dits_[measure_step.c_reg_] = mres;
                            st_.probs_[measure_step.c_reg_] = probs[mres];
                            set_measured_(measure_step.target_[0]);
                            break;
                        case QCircuit::MeasurementStep::Type::MEASURE_V_JOINT:
                            std::tie(mres, probs, states) = measure(
                                st_.psi_, h_tbl[measure_step.mats_hash_[0]],
                                target_rel_pos, d);
                            st_.psi_ = states[mres];
                            st_.dits_[measure_step.c_reg_] = mres;
                            st_.probs_[measure_step.c_reg_] = probs[mres];
                            for (idx target : measure_step.target_) {
                                set_measured_(target);
                            }
                            break;
                        case QCircuit::MeasurementStep::Type::MEASURE_ND:
                            std::tie(results, probs, st_.psi_) =
                                measure_seq(st_.psi_, target_rel_pos, d, false);
                            st_.dits_[measure_step.c_reg_] = results[0];
                            st_.probs_[measure_step.c_reg_] = probs[0];
                            break;
                        case QCircuit::MeasurementStep::Type::MEASURE_MANY_ND:
                            std::tie(results, probs, st_.psi_) =
                                measure_seq(st_.psi_, target_rel_pos, d, false);
                            std::copy(results.begin(), results.end(),
                                      std::next(st_.dits_.begin(),
                                                measure_step.c_reg_));
                            std::copy(probs.begin(), probs.end(),
                                      std::next(st_.probs_.begin(),
                                                measure_step.c_reg_));
                            break;
                        case QCircuit::MeasurementStep::Type::MEASURE_V_ND:
                        case QCircuit::MeasurementStep::Type::
                            MEASURE_V_JOINT_ND:
                            std::tie(mres, probs, states) = measure(
                                st_.psi_, h_tbl[measure_step.mats_hash_[0]],
                                target_rel_pos, d, false);
                            st_.psi_ = states[mres];
                            st_.dits_[measure_step.c_reg_] = mres;
                            st_.probs_[measure_step.c_reg_] = probs[mres];
                            break;
                        case QCircuit::MeasurementStep::Type::RESET:
                        case QCircuit::MeasurementStep::Type::RESET_MANY:
                            st_.psi_ = qpp::reset(st_.psi_, target_rel_pos, d);
                            break;
                        case QCircuit::MeasurementStep::Type::DISCARD:
                            std::tie(std::ignore, std::ignore, st_.psi_) =
                                measure_seq(st_.psi_, target_rel_pos, d);
                            set_measured_(measure_step.target_[0]);
                            break;
                        case QCircuit::MeasurementStep::Type::DISCARD_MANY:
                            std::tie(std::ignore, std::ignore, st_.psi_) =
                                measure_seq(st_.psi_, target_rel_pos, d);
                            for (idx target : measure_step.target_) {
                                set_measured_(target);
                            }
                            break;
                    } // end switch on measurement type
                },
                [&](const QCircuit::NOPStep&) {},
            },
            elem.get_step());

        return *this;
    }

    /**
     * \brief Executes one step in the quantum circuit description
     *
     * \note Do not override!
     *
     * \param it Iterator to the step to be executed
     * \return Reference to the current instance
     */
    QEngine& execute(const QCircuit::iterator& it) { return execute(*it); }

    /**
     * \brief Executes the entire quantum circuit description
     *
     * \param reps Number of repetitions
     * \param try_sampling If possible, try to sample from the output (instead
     * of measuring)
     * \return Reference to the current instance
     */
    virtual QEngine& execute(idx reps = 1, bool try_sampling = true) {
        this->reset(false);
        auto steps = (reps > 1 && try_sampling)
                         ? internal::canonical_form(*qc_ptr_)
                         : internal::circuit_as_iterators(*qc_ptr_);
        if (steps.empty()) {
            return *this;
        }

        idx num_steps = steps.size();

        // in the following, we will partition the circuit as
        // [0 ... first_measurement_discard_reset_pos ... end)

        // find the position of the first measurement/reset/discard step
        auto first_measurement_discard_reset_it =
            std::find_if(steps.begin(), steps.end(), [](auto&& elem) {
                return internal::is_measurement(elem) ||
                       internal::is_discard(elem) || internal::is_reset(elem);
            });
        idx first_measurement_discard_reset_pos =
            std::distance(steps.begin(), first_measurement_discard_reset_it);

        // decide if we can sample (every step after first_measurement_pos must
        // be a projective measurement)
        this->can_sample = (reps > 1) && try_sampling;
        for (idx i = first_measurement_discard_reset_pos;
             i < num_steps && this->can_sample; ++i) {
            if (!(internal::is_projective_measurement(steps[i])) ||
                internal::is_discard(steps[i])) {
                this->can_sample = false;
                break;
            }
        }

        // executes everything ONCE in the interval [0, first_measurement_pos)
        for (idx i = 0; i < first_measurement_discard_reset_pos; ++i) {
            execute(steps[i]);
        }
        // saves the state just before the first measurement
        auto current_engine_state = st_;

        // executes repeatedly everything in the remaining interval
        // [sampling_pos, num_steps)

        // can sample (every step must be a projective measurement)
        if (this->can_sample) {
            std::map<idx, idx> used_dits; // records the c <- q map
            bool measured = false;
            for (idx i = first_measurement_discard_reset_pos; i < num_steps;
                 ++i) {
                if (internal::is_projective_measurement(steps[i])) {
                    measured = true;
                    auto [_, target, c_regs] =
                        internal::extract_ctrl_target_c_reg(steps[i]);
                    for (idx q = 0; q < static_cast<idx>(target.size()); ++q) {
                        used_dits[c_regs[q]] = target[q];
                    }
                }
            }
            // at least one qudit was measured
            if (measured) {
                // build the vector of measured qudits that we must sample from
                std::vector<idx> sample_from;
                sample_from.reserve(this->get_dits().size());
                for (auto [dit, qubit] : used_dits) {
                    sample_from.emplace_back(qubit);
                }
                for (idx rep = 0; rep < reps - 1; ++rep) {
                    std::vector<idx> sample_res_restricted_support =
                        sample(current_engine_state.psi_, sample_from,
                               qc_ptr_->get_d());
                    // extend sample_res to full support
                    std::vector<idx> sample_res = this->get_dits();
                    idx i = 0;
                    for (auto [dit, qubit] : used_dits) {
                        sample_res[dit] = sample_res_restricted_support[i++];
                    }
                    ++stats_.data()[sample_res];
                }
                // execute the last repetition, so we can compute the state psi
                for (idx i = first_measurement_discard_reset_pos; i < num_steps;
                     ++i) {
                    execute(steps[i]);
                }
                std::vector<idx> m_res = get_dits();
                ++stats_.data()[m_res];
            }
        }
        // cannot sample
        else {
            for (idx rep = 0; rep < reps; ++rep) {
                // sets the state of the engine to the entry state
                st_ = current_engine_state;
                bool measured = false;
                for (idx i = first_measurement_discard_reset_pos; i < num_steps;
                     ++i) {
                    if (internal::is_measurement(steps[i])) {
                        measured = true;
                    }
                    execute(steps[i]);
                }
                // at least one qudit was measured
                if (measured) {
                    std::vector<idx> m_res = get_dits();
                    ++stats_.data()[m_res];
                }
            }
        }

        return *this;
    }

    /**
     * \brief qpp::IJSON::to_JSON() override
     *
     * Displays the state of the engine in JSON format
     *
     * \param enclosed_in_curly_brackets If true, encloses the result in curly
     * brackets
     * \return String containing the JSON representation of the state of the
     * engine
     */
    std::string to_JSON(bool enclosed_in_curly_brackets = true) const override {
        std::string result;

        if (enclosed_in_curly_brackets) {
            result += "{";
        }

        std::ostringstream ss;
        ss << "\"nq\": " << get_circuit().get_nq() << ", ";
        ss << "\"nc\": " << get_circuit().get_nc() << ", ";
        ss << "\"d\": " << get_circuit().get_d() << ", ";
        ss << "\"name\": ";
        if (get_circuit().get_name().has_value()) {
            ss << "\"" << get_circuit().get_name().value() << "\"";
        } else {
            ss << "null";
        }
        ss << ", ";
        ss << "\"sampling\": " << (this->can_sample ? "true" : "false") << ", ";
        ss << "\"measured/discarded (destructively)\": ";
        ss << disp(get_measured(), IOManipContainerOpts{}.set_sep(", "));
        result += ss.str() + ", ";

        ss.str("");
        ss.clear();
        ss << disp(get_non_measured(), IOManipContainerOpts{}.set_sep(", "));
        result += "\"non-measured (destructively)/non-discarded\": " + ss.str();

        ss.str("");
        ss.clear();
        result += ", \"last probs\": ";
        ss << disp(get_probs(), IOManipContainerOpts{}.set_sep(", "));
        result += ss.str();

        ss.str("");
        ss.clear();
        result += ", \"last dits\": ";
        ss << disp(get_dits(), IOManipContainerOpts{}.set_sep(", "));
        result += ss.str();

        ss.str("");
        ss.clear();

        // compute the statistics
        if (!stats_.data().empty()) {
            result += ", \"stats\": ";
            result += stats_.to_JSON();
        }

        if (enclosed_in_curly_brackets) {
            result += "}";
        }

        return result;
    }

  private:
    /**
     * \brief qpp::IDisplay::display() override
     *
     * Writes to the output stream a textual representation of the state of the
     * engine
     *
     * \param os Output stream passed by reference
     * \return Reference to the output stream
     */
    std::ostream& display(std::ostream& os) const override {
        /*
        os << "measured/discarded (destructive): " << disp(was_measured(),
        ", ")
           << '\n';
        os << "non-measured/non-discarded: " << disp(get_non_measured(), ",
        ")
           << '\n';
        */

        std::string engine_type = is_noisy() ? "[QNoisyEngine]" : "[QEngine]";
        os << engine_type;
        if (this->can_sample) {
            os << " (Sampling)";
        }
        os << '\n';

        os << "<QCircuit nq: " << get_circuit().get_nq()
           << ", nc: " << get_circuit().get_nc()
           << ", d: " << get_circuit().get_d();

        if (get_circuit().get_name().has_value()) {
            os << ", name: ";
            os << "\"" << get_circuit().get_name().value() << "\"";
        }
        os << ">\n";

        os << "last probs: "
           << disp(get_probs(), IOManipContainerOpts{}.set_sep(", ")) << '\n';
        os << "last dits: "
           << disp(get_dits(), IOManipContainerOpts{}.set_sep(", ")) << '\n';

        // compute the statistics
        if (!stats_.data().empty()) {
            os << stats_;
        }

        return os;
    }
}; /* class QEngine */

/**
 * \class qpp::QNoisyEngine
 * \brief Noisy quantum circuit engine, executes qpp::QCircuit
 * \see qpp::QEngine, qpp::QCircuit, qpp::NoiseBase
 *
 * Assumes an uncorrelated noise model that is applied to each non-measured
 * qubit before every non-measurement step in the logical circuit. To add
 * noise before a measurement, insert a no-op via qpp::QCircuit::nop().
 *
 * \tparam NoiseModel Quantum noise model, should be derived from qpp::NoiseBase
 */
template <typename NoiseModel>
class QNoisyEngine : public QEngine {
    NoiseModel noise_;                            ///< quantum noise model
    std::vector<std::vector<idx>> noise_results_; ///< noise results
  public:
    /**
     * \brief Constructs a noisy quantum engine out of a quantum circuit
     * description
     *
     * \param qc Quantum circuit description
     * \param noise Quantum noise model
     */
    explicit QNoisyEngine(const QCircuit& qc, const NoiseModel& noise)
        : QEngine{qc}, noise_{noise}, noise_results_(qc.get_step_count()) {
        // EXCEPTION CHECKS

        // check noise has the correct dimensionality
        if (qc.get_d() != noise.get_d()) {
            throw exception::DimsNotEqual("qpp::QNoisyEngine::QNoisyEngine()",
                                          "noise");
        }
        // END EXCEPTION CHECKS
    }

    using QEngine::execute;

    /**
     * \brief Executes one step in the quantum circuit description
     *
     * \param elem Step to be executed
     */
    QNoisyEngine& execute(const QCircuit::iterator::value_type& elem) override {
        // get the relative position of the target
        std::vector<idx> target_rel_pos = get_relative_pos_(get_non_measured());
        // if (elem.type_ != QCircuit::StepType::MEASUREMENT) {
        // apply the noise
        for (idx i : target_rel_pos) {
            st_.psi_ = noise_(st_.psi_, i);
            // record the Kraus operator that occurred
            noise_results_[elem.get_ip()].emplace_back(noise_.get_last_idx());
        }
        // }
        // execute the circuit step
        (void)QEngine::execute(elem);

        return *this;
    }

    /**
     * \brief Executes the entire quantum circuit description
     *
     * \param reps Number of repetitions
     * \return Reference to the current instance
     */
    QNoisyEngine& execute(idx reps = 1, bool = true) override {
        auto initial_engine_state = st_; // saves the engine entry state

        for (idx i = 0; i < reps; ++i) {
            // sets the state of the engine to the entry state
            st_ = initial_engine_state;

            for (auto&& elem : *qc_ptr_) {
                (void)execute(elem);
            }

            // we measured at least one qudit
            if (qc_ptr_->get_measurement_count() > 0) {
                std::vector<idx> m_res = get_dits();
                ++stats_.data()[m_res];
            }
        }

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
    QNoisyEngine& reset(bool reset_stats = true) override {
        QEngine::reset(reset_stats);
        noise_results_ = {};

        return *this;
    }

    // getters
    /**
     * \brief Vector of noise results obtained before every step in the circuit
     *
     * \note The first vector contains the noise measurement results
     * obtained before applying the first step in the circuit, and so on,
     * ordered by non-measured qudits. That is, the first element in the
     * vector corresponding to noise obtained before a given step in the
     * circuit represents the noise result obtained on the first non-measured
     * qudit etc.
     *
     * \return Vector of noise results
     */
    std::vector<std::vector<idx>> get_noise_results() const {
        return noise_results_;
    }

    /**
     * \brief \a qpp::QEngine::is_noisy() override
     *
     * \return True
     */
    bool is_noisy() const override { return true; }
    // end getters
}; /* class QNoisyEngine */

} /* namespace qpp */

#endif /* QPP_CLASSES_CIRCUITS_ENGINES_HPP_ */
