/*
 * This file is part of Quantum++.
 *
 * MIT License
 *
 * Copyright (c) 2013 - 2020 Vlad Gheorghiu.
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
 * \file classes/circuits/engines.h
 * \brief Qudit quantum engines
 */

#ifndef CLASSES_CIRCUITS_ENGINES_H_
#define CLASSES_CIRCUITS_ENGINES_H_

namespace qpp {
/**
 * \class qpp::QEngine
 * \brief Quantum circuit engine, executes qpp::QCircuit
 * \see qpp::QCircuit
 */
class QEngine : public IDisplay, public IJSON {
  protected:
    const QCircuit* qc_; ///< pointer to constant quantum circuit description

    /**
     * \class qpp::QEngine::state_
     * \brief Current state of the engine
     */
    struct state_ {
        const QCircuit* qc_;          ///< non-owning pointer to the parent
                                      ///< const quantum circuit description
        ket psi_{};                   ///< state vector
        std::vector<double> probs_{}; ///< measurement probabilities
        std::vector<idx> dits_{};     ///< classical dits
        std::vector<idx> subsys_{}; ///< keeps track of the measured subsystems,
                                    ///< re-label them after measurements

        /**
         * \brief Constructor
         *
         * \param qc Non-owning pointer to the parent const quantum circuit
         * description
         */
        explicit state_(const QCircuit* qc) : qc_{qc} {
            // EXECEPTION CHECKS

            if (qc->get_nq() == 0)
                throw exception::ZeroSize("qpp::QEngine::state_::reset()");
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
         */
        void reset() {
            psi_ = States::get_instance().zero(qc_->get_nq(), qc_->get_d());
            probs_ = std::vector<double>(qc_->get_nc(), 0);
            dits_ = std::vector<idx>(qc_->get_nc(), 0);
            subsys_ = std::vector<idx>(qc_->get_nq(), 0);
            std::iota(std::begin(subsys_), std::end(subsys_), 0);
        }
    } st_; ///< current state of the engine
    std::map<std::string, idx, internal::EqualSameSizeStringDits>
        stats_; ///< measurement statistics for multiple runs

    /**
     * \brief Marks qudit \a i as measured then re-label accordingly the
     * remaining non-measured qudits
     * \param i Qudit index
     */
    void set_measured_(idx i) {
        if (get_measured(i))
            throw exception::QuditAlreadyMeasured(
                "qpp::QEngine::set_measured_()");
        st_.subsys_[i] = static_cast<idx>(-1); // set qudit i to measured state
        for (idx m = i; m < qc_->get_nq(); ++m) {
            if (!get_measured(m)) {
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
            if (get_measured(v[i]))
                throw exception::QuditAlreadyMeasured(
                    "qpp::QEngine::get_relative_pos_()");
            v[i] = st_.subsys_[v[i]];
        }
        return v;
    }

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
        : qc_{std::addressof(qc)}, st_{qc_}, stats_{} {}

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
    virtual ~QEngine() = default;

    // getters
    /**
     * \brief Underlying quantum state
     *
     * \return Underlying quantum state
     */
    ket get_psi() const { return st_.psi_; }

    /**
     * \brief Vector with the values of the underlying classical dits
     *
     * \return Vector of underlying classical dits
     */
    std::vector<idx> get_dits() const { return st_.dits_; }

    /**
     * \brief Value of the classical dit at position \a i
     *
     * \param i Classical dit index
     * \return Value of the classical dit at position \a i
     */
    idx get_dit(idx i) const {
        if (i >= qc_->get_nc())
            throw exception::OutOfRange("qpp::QEngine::get_dit()");

        return st_.dits_[i];
    }

    /**
     * \brief Vector of underlying measurement outcome probabilities
     *
     * \note Those should be interpreted as conditional probabilities based on
     * the temporal order of the measurements, i.e. if we measure qubit 0, then
     * measure qubit 1, and finally qubit 2, the resulting vector of outcome
     * probabilities probs[2] should be interpreted as the conditional
     * probability of qubit 2 having the outcome it had given that qubit 1 and
     * qubit 0 had their given outcomes, respectively. As an example, if we
     * measure the qubit 0 followed by the qubit 1 of a maximally entangled
     * state \f$(|00\rangle + |11\rangle)/\sqrt{2}\f$, then the vector of
     * outcome probabilities will be [0.5, 1].
     *
     * \note The probability vector has the same length as the vector of
     * classical dits. If the measurement result is stored at the index
     * \a c_reg, then the outcome probability is automatically stored at
     * the same index \a c_reg in the probability vector.
     *
     * \return Vector of underlying measurement outcome probabilities
     */
    std::vector<double> get_probs() const { return st_.probs_; }

    /**
     * \brief Check whether qudit \a i was already measured
     *
     * \param i Qudit index
     * \return True if qudit \a i was already measured, false othwewise
     */
    bool get_measured(idx i) const {
        return st_.subsys_[i] == static_cast<idx>(-1);
    }

    /**
     * \brief Vector of already measured qudit indexes
     *
     * \return Vector of already measured qudit indexes
     */
    std::vector<idx> get_measured() const {
        std::vector<idx> result;
        for (idx i = 0; i < qc_->get_nq(); ++i)
            if (get_measured(i))
                result.emplace_back(i);

        return result;
    }

    /**
     * \brief Vector of non-measured qudit indexes
     *
     * \return Vector of non-measured qudit indexes
     */
    std::vector<idx> get_non_measured() const {
        std::vector<idx> result;
        for (idx i = 0; i < qc_->get_nq(); ++i)
            if (!get_measured(i))
                result.emplace_back(i);

        return result;
    }

    /**
     * \brief Quantum circuit description, lvalue ref qualifier
     *
     * \return Const reference to the underlying quantum circuit description
     */
    const QCircuit& get_circuit() const& noexcept { return *qc_; }

    /**
     * \brief Quantum circuit description, rvalue ref qualifier
     *
     * \return Copy of the underlying quantum circuit description
     */
    QCircuit get_circuit() const&& noexcept { return *qc_; }

    /**
     * \brief Measurement statistics for multiple runs
     *
     * \return Hash table with collected measurement statistics for multiple
     * runs, with hash key being the string representation of the vector of
     * measurement results and value being the number of occurrences (of the
     * vector of measurement results), with the most significant bit located at
     * index 0 (i.e. top/left).
     */
    std::map<std::string, idx, internal::EqualSameSizeStringDits>
    get_stats() const {
        return stats_;
    }
    // end getters

    // setters
    /**
     * \brief Sets the classical dit at position \a i
     *
     * \param i Classical dit index
     * \param value Classical dit value
     * \return Reference to the current instance
     */
    QEngine& set_dit(idx i, idx value) {
        if (i >= qc_->get_nc())
            throw exception::OutOfRange("qpp::QEngine::set_dit()");
        st_.dits_[i] = value;

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
        idx D = static_cast<idx>(std::llround(std::pow(qc_->get_d(), n)));
        if (static_cast<idx>(psi.rows()) != D)
            throw exception::DimsNotEqual("qpp::QEngine::set_psi()");
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
     * \return Reference to the current instance
     */
    QEngine& reset() {
        st_.reset();

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
        if (elem.value_type_qc_ != qc_)
            throw exception::InvalidIterator("qpp::QEngine::execute()");
        // the rest of exceptions are caught by the iterator::operator*()
        // END EXCEPTION CHECKS

        auto h_tbl = qc_->get_cmat_hash_tbl_();
        idx d = qc_->get_d();

        // gate step
        if (elem.type_ == QCircuit::StepType::GATE) {
            auto gates = qc_->get_gates_();

            idx q_ip =
                std::distance(std::begin(qc_->get_gates_()), elem.gates_ip_);

            std::vector<idx> ctrl_rel_pos;
            std::vector<idx> target_rel_pos =
                get_relative_pos_(gates[q_ip].target_);

            switch (gates[q_ip].gate_type_) {
                case QCircuit::GateType::NONE:
                    break;
                case QCircuit::GateType::SINGLE:
                case QCircuit::GateType::TWO:
                case QCircuit::GateType::THREE:
                case QCircuit::GateType::CUSTOM:
                    st_.psi_ = apply(st_.psi_, h_tbl[gates[q_ip].gate_hash_],
                                     target_rel_pos, d);
                    break;
                case QCircuit::GateType::FAN:
                    for (idx m = 0; m < gates[q_ip].target_.size(); ++m)
                        st_.psi_ =
                            apply(st_.psi_, h_tbl[gates[q_ip].gate_hash_],
                                  {target_rel_pos[m]}, d);
                    break;
                case QCircuit::GateType::SINGLE_CTRL_SINGLE_TARGET:
                case QCircuit::GateType::SINGLE_CTRL_MULTIPLE_TARGET:
                case QCircuit::GateType::MULTIPLE_CTRL_SINGLE_TARGET:
                case QCircuit::GateType::MULTIPLE_CTRL_MULTIPLE_TARGET:
                case QCircuit::GateType::CUSTOM_CTRL:
                    ctrl_rel_pos = get_relative_pos_(gates[q_ip].ctrl_);
                    st_.psi_ = applyCTRL(
                        st_.psi_, h_tbl[gates[q_ip].gate_hash_], ctrl_rel_pos,
                        target_rel_pos, d, gates[q_ip].shift_);
                    break;
                case QCircuit::GateType::SINGLE_cCTRL_SINGLE_TARGET:
                case QCircuit::GateType::SINGLE_cCTRL_MULTIPLE_TARGET:
                case QCircuit::GateType::MULTIPLE_cCTRL_SINGLE_TARGET:
                case QCircuit::GateType::MULTIPLE_cCTRL_MULTIPLE_TARGET:
                case QCircuit::GateType::CUSTOM_cCTRL:
                    if (st_.dits_.empty()) {
                        st_.psi_ =
                            apply(st_.psi_, h_tbl[gates[q_ip].gate_hash_],
                                  target_rel_pos, d);
                    } else {
                        bool should_apply = true;
                        idx first_dit;
                        // we have a shift
                        if (!gates[q_ip].shift_.empty()) {
                            first_dit = (st_.dits_[(gates[q_ip].ctrl_)[0]] +
                                         gates[q_ip].shift_[0]) %
                                        d;
                            for (idx m = 1; m < gates[q_ip].ctrl_.size(); ++m) {
                                if ((st_.dits_[(gates[q_ip].ctrl_)[m]] +
                                     gates[q_ip].shift_[m]) %
                                        d !=
                                    first_dit) {
                                    should_apply = false;
                                    break;
                                }
                            }
                        }
                        // no shift
                        else {
                            first_dit = st_.dits_[(gates[q_ip].ctrl_)[0]];
                            for (idx m = 1; m < gates[q_ip].ctrl_.size(); ++m) {
                                if (st_.dits_[(gates[q_ip].ctrl_)[m]] !=
                                    first_dit) {
                                    should_apply = false;
                                    break;
                                }
                            }
                        }
                        if (should_apply) {
                            st_.psi_ = apply(
                                st_.psi_,
                                powm(h_tbl[gates[q_ip].gate_hash_], first_dit),
                                target_rel_pos, d);
                        }
                    }
                    break;
            } // end switch on gate type
        }     // end if gate step

        // measurement step
        else if (elem.type_ == QCircuit::StepType::MEASUREMENT) {
            auto measurements = qc_->get_measurements_();
            idx m_ip = std::distance(std::begin(qc_->get_measurements_()),
                                     elem.measurements_ip_);

            std::vector<idx> target_rel_pos =
                get_relative_pos_(measurements[m_ip].target_);

            std::vector<idx> resZ;
            double probZ;

            idx mres = 0;
            std::vector<double> probs;
            std::vector<cmat> states;

            switch (measurements[m_ip].measurement_type_) {
                case QCircuit::MeasureType::NONE:
                    break;
                case QCircuit::MeasureType::MEASURE_Z:
                    std::tie(resZ, probZ, st_.psi_) =
                        measure_seq(st_.psi_, target_rel_pos, d);
                    st_.dits_[measurements[m_ip].c_reg_] = resZ[0];
                    st_.probs_[measurements[m_ip].c_reg_] = probZ;
                    set_measured_(measurements[m_ip].target_[0]);
                    break;
                case QCircuit::MeasureType::MEASURE_Z_MANY:
                    std::tie(resZ, probZ, st_.psi_) =
                        measure_seq(st_.psi_, target_rel_pos, d);
                    st_.dits_[measurements[m_ip].c_reg_] = multiidx2n(
                        resZ, std::vector<idx>(target_rel_pos.size(), d));
                    st_.probs_[measurements[m_ip].c_reg_] = probZ;
                    for (auto&& elem : measurements[m_ip].target_)
                        set_measured_(elem);
                    break;
                case QCircuit::MeasureType::MEASURE_V:
                    std::tie(mres, probs, states) = measure(
                        st_.psi_, h_tbl[measurements[m_ip].mats_hash_[0]],
                        target_rel_pos, d);
                    st_.psi_ = states[mres];
                    st_.dits_[measurements[m_ip].c_reg_] = mres;
                    st_.probs_[measurements[m_ip].c_reg_] = probs[mres];
                    set_measured_(measurements[m_ip].target_[0]);
                    break;
                case QCircuit::MeasureType::MEASURE_V_MANY:
                    std::tie(mres, probs, states) = measure(
                        st_.psi_, h_tbl[measurements[m_ip].mats_hash_[0]],
                        target_rel_pos, d);
                    st_.psi_ = states[mres];
                    st_.dits_[measurements[m_ip].c_reg_] = mres;
                    st_.probs_[measurements[m_ip].c_reg_] = probs[mres];
                    for (auto&& elem : measurements[m_ip].target_)
                        set_measured_(elem);
                    break;
                case QCircuit::MeasureType::MEASURE_Z_ND:
                    std::tie(resZ, probZ, st_.psi_) =
                        measure_seq(st_.psi_, target_rel_pos, d, false);
                    st_.dits_[measurements[m_ip].c_reg_] = resZ[0];
                    st_.probs_[measurements[m_ip].c_reg_] = probZ;
                    break;
                case QCircuit::MeasureType::MEASURE_Z_MANY_ND:
                    std::tie(resZ, probZ, st_.psi_) =
                        measure_seq(st_.psi_, target_rel_pos, d, false);
                    st_.dits_[measurements[m_ip].c_reg_] = multiidx2n(
                        resZ, std::vector<idx>(target_rel_pos.size(), d));
                    st_.probs_[measurements[m_ip].c_reg_] = probZ;
                    break;
                case QCircuit::MeasureType::MEASURE_V_ND:
                case QCircuit::MeasureType::MEASURE_V_MANY_ND:
                    std::tie(mres, probs, states) = measure(
                        st_.psi_, h_tbl[measurements[m_ip].mats_hash_[0]],
                        target_rel_pos, d, false);
                    st_.psi_ = states[mres];
                    st_.dits_[measurements[m_ip].c_reg_] = mres;
                    st_.probs_[measurements[m_ip].c_reg_] = probs[mres];
                    break;
                case QCircuit::MeasureType::RESET:
                case QCircuit::MeasureType::RESET_MANY:
                    st_.psi_ = qpp::reset(st_.psi_, target_rel_pos, d);
                    break;
                case QCircuit::MeasureType::DISCARD:
                    std::tie(std::ignore, std::ignore, st_.psi_) =
                        measure_seq(st_.psi_, target_rel_pos, d);
                    set_measured_(measurements[m_ip].target_[0]);
                    break;
                case QCircuit::MeasureType::DISCARD_MANY:
                    std::tie(std::ignore, std::ignore, st_.psi_) =
                        measure_seq(st_.psi_, target_rel_pos, d);
                    for (auto&& elem : measurements[m_ip].target_)
                        set_measured_(elem);
                    break;
            } // end switch on measurement type
        }     // end else if measurement step

        // no-op
        else if (elem.type_ == QCircuit::StepType::NOP) {
        }

        // otherwise
        else {
        }

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
     * \note Do not override!
     *
     * \param reps Number of repetitions
     * \param clear_stats Resets the collected measurement statistics hash table
     * before the run
     * \return Reference to the current instance
     */
    QEngine& execute(idx reps = 1, bool clear_stats = true) {
        auto initial_engine_state = st_; // saves the engine entry state

        if (clear_stats)
            reset_stats();

        for (idx i = 0; i < reps; ++i) {
            // sets the state of the engine to the entry state
            st_ = initial_engine_state;

            for (auto&& elem : *qc_)
                (void) execute(elem);

            // we measured at least one qudit
            if (qc_->get_measurement_count() > 0) {
                std::vector<idx> m_res = get_dits();

                std::stringstream ss;
                ss << disp(m_res, " ", "", "");
                ++stats_[ss.str()];
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

        if (enclosed_in_curly_brackets)
            result += "{";

        std::ostringstream ss;
        ss << disp(get_measured(), ", ");
        result += "\"measured/discarded (destructive)\" : " + ss.str() + ", ";

        ss.str("");
        ss.clear();
        ss << disp(get_non_measured(), ", ");
        result += "\"non-measured/non-discarded\" : " + ss.str();

        ss.str("");
        ss.clear();
        result += ", \"last probs\": ";
        ss << disp(get_probs(), ", ");
        result += ss.str();

        ss.str("");
        ss.clear();
        result += ", \"last dits\": ";
        ss << disp(get_dits(), ", ");
        result += ss.str();

        ss.str("");
        ss.clear();

        // compute the statistics
        if (!stats_.empty()) {
            result += ", \"stats\": {";
            idx reps = 0;
            for (auto&& elem : get_stats())
                reps += elem.second;
            result += "\"reps\": " + std::to_string(reps) + ", ";
            result += "\"outcomes\": " + std::to_string(stats_.size()) + ", ";

            std::vector<idx> dits_dims(qc_->get_nc(), qc_->get_d());
            bool is_first = true;
            for (auto&& elem : get_stats()) {
                if (is_first)
                    is_first = false;
                else
                    ss << ", ";
                ss << "\""
                   << "[" << elem.first << "]"
                   << "\" : " << elem.second;
            }
            ss << '}';
            result += ss.str();
        }

        if (enclosed_in_curly_brackets)
            result += "}";

        return result;
    }

  private:
    /**
     * \brief qpp::IDisplay::display() override
     *
     * Writes to the output stream a textual representation of the state of
     * the engine
     *
     * \param os Output stream passed by reference
     * \return Reference to the output stream
     */
    std::ostream& display(std::ostream& os) const override {
        os << "measured/discarded (destructive): " << disp(get_measured(), ", ")
           << '\n';
        os << "non-measured/non-discarded: " << disp(get_non_measured(), ", ")
           << '\n';
        os << "last probs: " << disp(get_probs(), ", ") << '\n';
        os << "last dits: " << disp(get_dits(), ", ");

        // compute the statistics
        if (!stats_.empty()) {
            idx reps = 0;
            for (auto&& elem : get_stats())
                reps += elem.second;
            os << "\nstats:\n";
            os << '\t' << "reps: " << reps << '\n';
            os << '\t' << "outcomes: " << stats_.size() << '\n';
            std::vector<idx> dits_dims(qc_->get_nc(), qc_->get_d());
            bool is_first = true;
            for (auto&& elem : get_stats()) {
                if (is_first)
                    is_first = false;
                else
                    os << '\n';
                os << '\t' << "[" << elem.first << "]"
                   << ": " << elem.second;
            }
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
 * qubit before every non-measurement step in the logical circuit. To add noise
 * before a measurement, insert a no-op via qpp::QCircuit::nop().
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
        if (qc.get_d() != noise.get_d())
            throw exception::DimsNotEqual("qpp::QNoisyEngine::QNoisyEngine()");
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
        if (elem.type_ != QCircuit::StepType::MEASUREMENT) {
            // apply the noise
            for (auto&& i : target_rel_pos) {
                st_.psi_ = noise_(st_.psi_, i);
                // record the Kraus operator that occurred
                noise_results_[elem.ip_].emplace_back(noise_.get_last_idx());
            }
        }
        // execute the circuit step
        (void) QEngine::execute(elem);

        return *this;
    }

    // getters
    /**
     * \brief Vector of noise results obtained before every step in the
     * circuit
     *
     * \note The first vector contains the noise measurement results obtained
     * before applying the first step in the circuit, and so on, ordered by
     * non-measured qudits. That is, the first element in the vector
     * corresponding to noise obtained before a given step in the circuit
     * represents the noise result obtained on the first non-measured qudit etc.
     *
     * \return Vector of noise results
     */
    std::vector<std::vector<idx>> get_noise_results() const {
        return noise_results_;
    }
    // end getters
}; /* class QNoisyEngine */

} /* namespace qpp */

#endif /* CLASSES_CIRCUITS_ENGINES_H_ */
