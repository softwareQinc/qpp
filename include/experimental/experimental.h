/*
 * This file is part of Quantum++.
 *
 * MIT License
 *
 * Copyright (c) 2013 - 2019 Vlad Gheorghiu (vgheorgh@gmail.com)
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
 * \file experimental/experimental.h
 * \brief Experimental/test functions/classes
 */

#ifndef EXPERIMENTAL_EXPERIMENTAL_H_
#define EXPERIMENTAL_EXPERIMENTAL_H_

namespace qpp {
/**
 * \namespace qpp::experimental
 * \brief Experimental/test functions/classes, do not use or modify
 */
namespace experimental {

// TODO: (in progress) add QFT/TFQ as a "gate" type, what about computing
//  depths?!
//  best is to apply it gate by gate in QCircuitDescription::gate(...)

// TODO: perform exception checking before run() (such as wrong idx on apply or
//  out of range ctrl/targets)
//  i.e. ROBUST EXCEPTION CHECKING, something like a sanitize function! Do this
//  in gate() and not in run()

class QCircuitDescription : public IDisplay {
    friend class QCircuit;
    idx nq_;                               ///< number of qudits
    idx nc_;                               ///< number of classical "dits"
    idx d_;                                ///< dimension
    std::vector<idx> measurement_steps_{}; ///< keeps track of where the
    ///< measurement take place
    std::string name_;           ///< optional circuit name
    std::vector<bool> measured_; ///< keeps track of the measured qudits
    idx step_cnt_;               ///< step counter

    /**
     * \brief Type of operation being executed at one step
     */
    enum class GateType {
        NONE, ///< signals no gate

        SINGLE, ///< unitary gate on a single qudit

        TWO, ///< unitary gate on 2 qudits

        THREE, ///< unitary gate on 3 qudits

        CUSTOM, ///< custom gate on multiple qudits

        FAN, ///< same unitary gate on multiple qudits

        QFT, ///< quantum Fourier transform,

        TFQ, ///< quantum inverse Fourier transform,

        SINGLE_CTRL_SINGLE_TARGET, ///< controlled 1 qudit unitary gate with
        ///< one control and one target

        SINGLE_CTRL_MULTIPLE_TARGET, ///< controlled 1 qudit unitary gate with
        ///< one control and multiple targets

        MULTIPLE_CTRL_SINGLE_TARGET, ///< controlled 1 qudit unitary gate with
        ///< multiple controls and single target

        MULTIPLE_CTRL_MULTIPLE_TARGET, ///< controlled 1 qudit unitary gate with
        ///< multiple controls and multiple
        ///< targets

        CUSTOM_CTRL, ///< custom controlled gate with multiple
        ///< controls and multiple targets

        SINGLE_cCTRL_SINGLE_TARGET, ///< controlled 1 qudit unitary gate with
        ///< one classical control and one target

        SINGLE_cCTRL_MULTIPLE_TARGET, ///< controlled 1 qudit unitary gate with
        ///< one classical control and multiple targets

        MULTIPLE_cCTRL_SINGLE_TARGET, ///< controlled 1 qudit unitary gate with
        ///< multiple classical controls and single target

        MULTIPLE_cCTRL_MULTIPLE_TARGET, ///< controlled 1 qudit unitary gate
        ///< with multiple classical controls and multiple targets

        CUSTOM_cCTRL, ///< custom controlled gate with multiple
        ///< controls and multiple targets
    };
    enum class MeasureType {
        NONE, ///< signals no measurement

        MEASURE_Z, ///< Z measurement of single qudit

        MEASURE_V, ///< measurement of single qudit in the orthonormal basis
        ///< or rank-1 POVM specified by the columns of matrix \a V

        MEASURE_V_MANY, ///< measurement of multiple qudits in the orthonormal
        ///< basis or rank-1 POVM specified by the columns of matrix \a V
    };

    /**
     * \brief One step consisting of unitary gates in the circuit
     */
    struct GateStep {
        GateType gate_type_ = GateType::NONE; ///< gate type
        cmat gate_;                           ///< gate
        std::vector<idx> ctrl_;               ///< control
        std::vector<idx> target_; ///< target where the gate is being applied
        idx step_no_;             ///< step number
        std::string name_;        ///< custom name of the step
        GateStep() = default;
        GateStep(GateType gate_type, const cmat& gate,
                 const std::vector<idx>& ctrl, const std::vector<idx>& target,
                 idx step_no, const std::string& name = "")
            : gate_type_{gate_type}, gate_{gate}, ctrl_{ctrl}, target_{target},
              step_no_{step_no}, name_{name} {}
    };

    struct MeasureStep {
        MeasureType measurement_type_ = MeasureType::NONE; ///< measurement type
        std::vector<cmat> mats_;  ///< matrix/matrices being measured
        std::vector<idx> target_; ///< target where the measurement is applied
        idx c_reg_{}; ///< index of the classical register where the measurement
        ///< result is being stored
        idx step_no_;      ///< step number
        std::string name_; ///< custom name of the step
        MeasureStep() = default;
        MeasureStep(MeasureType measurement_type, const std::vector<cmat>& mats,
                    const std::vector<idx>& target, idx c_reg, idx step_no,
                    const std::string& name = "")
            : measurement_type_{measurement_type}, mats_{mats}, target_{target},
              c_reg_{c_reg}, step_no_{step_no}, name_{name} {}
    };

    friend std::ostream& operator<<(std::ostream& os,
                                    const GateType& gate_type) {
        switch (gate_type) {
        case GateType::NONE:
            os << "GATE NONE";
            break;
        case GateType::SINGLE:
            os << "SINGLE";
            break;
        case GateType::TWO:
            os << "TWO";
            break;
        case GateType::THREE:
            os << "THREE";
            break;
        case GateType::FAN:
            os << "FAN";
            break;
        case GateType::QFT:
            os << "QFT";
            break;
        case GateType::TFQ:
            os << "TFQ";
            break;
        case GateType::CUSTOM:
            os << "CUSTOM";
            break;
        case GateType::SINGLE_CTRL_SINGLE_TARGET:
            os << "SINGLE_CTRL_SINGLE_TARGET";
            break;
        case GateType::SINGLE_CTRL_MULTIPLE_TARGET:
            os << "SINGLE_CTRL_MULTIPLE_TARGET";
            break;
        case GateType::MULTIPLE_CTRL_SINGLE_TARGET:
            os << "MULTIPLE_CTRL_SINGLE_TARGET";
            break;
        case GateType::MULTIPLE_CTRL_MULTIPLE_TARGET:
            os << "MULTIPLE_CTRL_MULTIPLE_TARGET";
            break;
        case GateType::CUSTOM_CTRL:
            os << "CUSTOM_CTRL";
            break;
        case GateType::SINGLE_cCTRL_SINGLE_TARGET:
            os << "SINGLE_cCTRL_SINGLE_TARGET";
            break;
        case GateType::SINGLE_cCTRL_MULTIPLE_TARGET:
            os << "SINGLE_cCTRL_MULTIPLE_TARGET";
            break;
        case GateType::MULTIPLE_cCTRL_SINGLE_TARGET:
            os << "MULTIPLE_cCTRL_SINGLE_TARGET";
            break;
        case GateType::MULTIPLE_cCTRL_MULTIPLE_TARGET:
            os << "MULTIPLE_cCTRL_MULTIPLE_TARGET";
            break;
        case GateType::CUSTOM_cCTRL:
            os << "CUSTOM_cCTRL";
            break;
        }

        return os;
    }

    friend std::ostream& operator<<(std::ostream& os,
                                    const MeasureType& measure_type) {
        switch (measure_type) {
        case MeasureType::NONE:
            os << "MEASURE NONE";
            break;
        case MeasureType::MEASURE_Z:
            os << "MEASURE_Z";
            break;
        case MeasureType::MEASURE_V:
            os << "MEASURE_V";
            break;
        case MeasureType::MEASURE_V_MANY:
            os << "MEASURE_V_MANY";
            break;
        }

        return os;
    }

    friend std::ostream& operator<<(std::ostream& os,
                                    const GateStep& gate_step) {
        // os << "gate type = " << ;
        os << gate_step.gate_type_ << ", ";
        if (gate_step.gate_type_ >= GateType::SINGLE_CTRL_SINGLE_TARGET)
            os << "ctrl = " << disp(gate_step.ctrl_, ",") << ", ";
        os << "target = " << disp(gate_step.target_, ",") << ", ";
        os << "name = " << '\"' << gate_step.name_ << '\"';

        return os;
    }

    friend std::ostream& operator<<(std::ostream& os,
                                    const MeasureStep& measure_step) {
        // os << "measurement type = ";
        os << measure_step.measurement_type_ << ", ";
        os << "target = " << disp(measure_step.target_, ",") << ", ";
        os << "c_reg = " << measure_step.c_reg_ << ", ";
        os << "name = " << '\"' << measure_step.name_ << '\"';
        os << " ";

        return os;
    }

    ///< quantum circuit representation
    std::vector<GateStep> gates_{};
    std::vector<MeasureStep> measurements_{};

  public:
    QCircuitDescription(idx nq, idx nc = 0, idx d = 2,
                        const std::string& name = "")
        : nq_{nq}, nc_{nc}, d_{d}, name_{name},
          measured_(nq, false), step_cnt_{0} {}

    // getters
    // number of qudits
    idx get_nq() const { return nq_; }

    // number of classical dits
    idx get_nc() const { return nc_; }

    // dimension
    idx get_d() const { return d_; }

    // measurement steps
    std::vector<idx> get_measurement_steps() const {
        return measurement_steps_;
    }

    // circuit name
    std::string get_name() const { return name_; }

    // return true if qudit i was measured, false otherwise
    idx was_measured(idx i) const { return measured_[i]; }

    // qudits that were measured
    std::vector<idx> get_measured() const {
        std::vector<idx> result;
        for (idx i = 0; i < nq_; ++i)
            if (was_measured(i))
                result.emplace_back(i);

        return result;
    }

    // qudits that were not (yet) measured
    std::vector<idx> get_non_measured() const {
        std::vector<idx> result;
        for (idx i = 0; i < nq_; ++i)
            if (!was_measured(i))
                result.emplace_back(i);

        return result;
    }

    idx get_gate_count() const { return gates_.size(); }

    idx get_measurement_count() const { return measurements_.size(); }

    idx get_total_count() const { return step_cnt_; }
    // end getters

    // single gate single qudit
    QCircuitDescription& gate(const cmat& U, idx i,
                              const std::string& name = "") {
        gates_.emplace_back(GateType::SINGLE, U, std::vector<idx>{},
                            std::vector<idx>{i}, step_cnt_++, name);

        return *this;
    }

    // single gate 2 qudits
    QCircuitDescription& gate(const cmat& U, idx i, idx j,
                              const std::string& name = "") {
        gates_.emplace_back(GateType::TWO, U, std::vector<idx>{},
                            std::vector<idx>{i, j}, step_cnt_++, name);

        return *this;
    }
    // single gate 3 qudits
    QCircuitDescription& gate(const cmat& U, idx i, idx j, idx k,
                              const std::string& name = "") {
        gates_.emplace_back(GateType::THREE, U, std::vector<idx>{},
                            std::vector<idx>{i, j, k}, step_cnt_++, name);

        return *this;
    }

    // multiple qudits same gate
    QCircuitDescription& gate_fan(const cmat& U, const std::vector<idx>& target,
                                  const std::string& name = "") {
        gates_.emplace_back(GateType::FAN, U, std::vector<idx>{}, target,
                            step_cnt_++, name);

        return *this;
    }

    // multiple qudits same gate on ALL non-measured qudits
    QCircuitDescription& gate_fan(const cmat& U, const std::string& name = "") {
        gates_.emplace_back(GateType::FAN, U, std::vector<idx>{},
                            get_non_measured(), step_cnt_++, name);

        return *this;
    }

    // custom gate
    QCircuitDescription& gate(const cmat& U, const std::vector<idx>& target,
                              const std::string& name = "") {
        gates_.emplace_back(GateType::CUSTOM, U, std::vector<idx>{}, target,
                            step_cnt_++, name);

        return *this;
    }

    // quantum Fourier transform
    QCircuitDescription& QFT(const std::vector<idx>& target) {
        gates_.emplace_back(GateType::QFT, cmat{}, std::vector<idx>{}, target,
                            step_cnt_++, "QFT");

        return *this;
    }

    // quantum inverse Fourier transform
    QCircuitDescription& TFQ(const std::vector<idx>& target) {
        gates_.emplace_back(GateType::TFQ, cmat{}, std::vector<idx>{}, target,
                            step_cnt_++, "TFQ");

        return *this;
    }

    // single ctrl single target
    QCircuitDescription& CTRL(const cmat& U, idx ctrl, idx target,
                              const std::string& name = "") {
        gates_.emplace_back(GateType::SINGLE_CTRL_SINGLE_TARGET, U,
                            std::vector<idx>{ctrl}, std::vector<idx>{target},
                            step_cnt_++, name);

        return *this;
    }

    // single ctrl multiple targets
    QCircuitDescription& CTRL(const cmat& U, idx ctrl,
                              const std::vector<idx>& target,
                              const std::string& name = "") {
        gates_.emplace_back(GateType::SINGLE_CTRL_MULTIPLE_TARGET, U,
                            std::vector<idx>{ctrl}, target, step_cnt_++, name);

        return *this;
    }

    // multiple ctrl single target
    QCircuitDescription& CTRL(const cmat& U, const std::vector<idx>& ctrl,
                              idx target, const std::string& name = "") {
        gates_.emplace_back(GateType::MULTIPLE_CTRL_SINGLE_TARGET, U, ctrl,
                            std::vector<idx>{target}, step_cnt_++, name);

        return *this;
    }

    // multiple ctrl multiple targets
    // FIXME
    QCircuitDescription& CTRL(const cmat& U, const std::vector<idx>& ctrl,
                              const std::vector<idx>& target,
                              const std::string& name = "") {
        gates_.emplace_back(GateType::MULTIPLE_CTRL_MULTIPLE_TARGET, U, ctrl,
                            std::vector<idx>{target}, step_cnt_++, name);

        return *this;
    }

    //  custom controlled gate with multiple controls and multiple targets
    QCircuitDescription& CTRL_custom(const cmat& U,
                                     const std::vector<idx>& ctrl,
                                     const std::vector<idx>& target,
                                     const std::string& name = "") {
        gates_.emplace_back(GateType::CUSTOM_CTRL, U, ctrl, target, step_cnt_++,
                            name);

        return *this;
    }

    // FIXME , use the corresponding dits
    // single ctrl single target
    QCircuitDescription& cCTRL(const cmat& U, idx ctrl, idx target,
                               const std::string& name = "") {
        gates_.emplace_back(GateType::SINGLE_cCTRL_SINGLE_TARGET, U,
                            std::vector<idx>{ctrl}, std::vector<idx>{target},
                            step_cnt_++, name);

        return *this;
    }

    // single ctrl multiple targets
    QCircuitDescription& cCTRL(const cmat& U, idx ctrl,
                               const std::vector<idx>& target,
                               const std::string& name = "") {
        gates_.emplace_back(GateType::SINGLE_cCTRL_MULTIPLE_TARGET, U,
                            std::vector<idx>{ctrl}, target, step_cnt_++, name);

        return *this;
    }

    // multiple ctrl single target
    QCircuitDescription& cCTRL(const cmat& U, const std::vector<idx>& ctrl,
                               idx target, const std::string& name = "") {
        gates_.emplace_back(GateType::MULTIPLE_cCTRL_SINGLE_TARGET, U, ctrl,
                            std::vector<idx>{target}, step_cnt_++, name);

        return *this;
    }

    // multiple ctrl multiple targets
    QCircuitDescription& cCTRL(const cmat& U, const std::vector<idx>& ctrl,
                               const std::vector<idx>& target,
                               const std::string& name = "") {
        gates_.emplace_back(GateType::MULTIPLE_cCTRL_MULTIPLE_TARGET, U, ctrl,
                            std::vector<idx>{target}, step_cnt_++, name);

        return *this;
    }

    //  custom controlled gate with multiple controls and multiple targets
    QCircuitDescription& cCTRL_custom(const cmat& U,
                                      const std::vector<idx>& ctrl,
                                      const std::vector<idx>& target,
                                      const std::string& name = "") {
        gates_.emplace_back(GateType::CUSTOM_cCTRL, U, ctrl, target,
                            step_cnt_++, name);

        return *this;
    }

    // Z measurement of single qudit
    QCircuitDescription& measureZ(idx i, idx c_reg,
                                  const std::string& name = "") {
        // EXCEPTION CHECKS

        // measuring non-existing qudit
        if (i >= nq_)
            throw exception::OutOfRange("qpp::QCircuitDescription::measureZ()");
        // trying to put the result into an non-existing classical slot
        if (c_reg >= nc_)
            throw exception::OutOfRange("qpp::QCircuitDescription::measureZ()");
        // qudit was measured before
        if (measured_[i] == true)
            throw exception::QuditAlreadyMeasured(
                "qpp:QCircuitDescription::measureZ");
        // END EXCEPTION CHECKS

        measured_[i] = true;
        measurements_.emplace_back(MeasureType::MEASURE_Z, std::vector<cmat>{},
                                   std::vector<idx>{i}, c_reg, step_cnt_++,
                                   name);
        measurement_steps_.emplace_back(gates_.size());
        /*try {
        } catch (std::exception& e) {
            std::cerr << "IN qpp::QCircuit::measureZ()\n";
            throw;
        }*/

        return *this;
    }

    // measurement of single qudit in the orthonormal basis or rank-1 POVM
    // specified by the columns of matrix V
    QCircuitDescription& measureV(const cmat& V, idx i, idx c_reg,
                                  const std::string& name = "") {
        // EXCEPTION CHECKS

        // measuring non-existing qudit
        if (i >= nq_)
            throw exception::OutOfRange("qpp::QCircuitDescription::measureV()");
        // trying to put the result into an non-existing classical slot
        if (c_reg >= nc_)
            throw exception::OutOfRange("qpp::QCircuitDescription::measureV()");
        // qudit was measured before
        if (measured_[i] == true)
            throw exception::QuditAlreadyMeasured(
                "qpp:QCircuitDescription::measureV");
        // END EXCEPTION CHECKS

        measured_[i] = true;
        measurements_.emplace_back(MeasureType::MEASURE_V, std::vector<cmat>{V},
                                   std::vector<idx>{i}, c_reg, step_cnt_++,
                                   name);
        measurement_steps_.emplace_back(gates_.size());

        return *this;
    }

    // measurement of multiple qudits in the orthonormal basis or rank-1 POVM
    // specified by the columns of matrix V
    QCircuitDescription& measureV(const cmat& V, const std::vector<idx>& target,
                                  idx c_reg, const std::string& name = "") {
        // EXCEPTION CHECKS

        // measuring non-existing qudit
        for (auto&& i : target)
            if (i >= nq_)
                throw exception::OutOfRange(
                    "qpp::QCircuitDescription::measureV()");
        // trying to put the result into an non-existing classical slot
        if (c_reg >= nc_)
            throw exception::OutOfRange("qpp::QCircuitDescription::measureV()");
        // qudit was measured before
        for (auto&& i : target)
            if (measured_[i] == true)
                throw exception::QuditAlreadyMeasured(
                    "qpp::QCircuitDescription::measureV");
        // END EXCEPTION CHECKS

        for (auto&& i : target)
            measured_[i] = true;
        measurements_.emplace_back(MeasureType::MEASURE_V_MANY,
                                   std::vector<cmat>{V}, target, c_reg,
                                   step_cnt_++, name);
        measurement_steps_.emplace_back(gates_.size());

        return *this;
    }

    std::ostream& display(std::ostream& os) const override {
        os << "nq = " << nq_ << ", nc = " << nc_ << ", d = " << d_;

        if (name_ != "") // if the circuit is named
            os << ", name = \"" << name_ << "\"\n";
        else
            os << ", name = \"\"\n";

        idx gates_size = gates_.size();
        idx measurements_size = measurements_.size();

        idx m_ip = 0; // measurement instruction pointer
        for (idx i = 0; i <= gates_size; ++i) {

            // check for measurements
            if (m_ip < measurements_size) {
                idx m_step = measurement_steps_[m_ip];
                // we have a measurement at step i
                if (m_step == i) {
                    while (measurement_steps_[m_ip] == m_step) {
                        os << std::left;
                        os << std::setw(8) << measurements_[m_ip].step_no_;
                        os << std::right;
                        os << "|> " << measurements_[m_ip] << '\n';
                        if (++m_ip == measurements_size)
                            break;
                    }
                }
            }

            // check for gates
            if (i < gates_size) {
                os << std::left;
                os << std::setw(8) << gates_[i].step_no_;
                os << std::right;
                os << gates_[i] << '\n';
            }
        }

        os << "measurement steps: " << disp(measurement_steps_, ",") << '\n';

        os << std::boolalpha;
        os << "measured qudits: " << disp(measured_, ",") << '\n';
        os << std::noboolalpha;

        os << "measured positions: " << disp(get_measured(), ",") << '\n';
        os << "non-measured positions: " << disp(get_non_measured(), ",")
           << '\n';

        return os;
    }

}; /* class QCircuitDescription */

class QCircuit : public IDisplay {
    const QCircuitDescription qcd_; ///< quantum circuit description
    ket psi_;                       ///< state vector
    std::vector<idx> dits_;         ///< classical dits
    std::vector<double> probs_;     ///< measurement probabilities
    std::vector<idx> subsys_;       ///< keeps track of the measured subsystems,
                                    ///< relabel them after measurements

    idx m_ip_; ///< measurement instruction pointer
    idx q_ip_; ///< quantum gates instruction pointer
    idx ip_;   ///< combined (measurements and gates) instruction pointer

    // mark qudit i as measured then re-label accordingly the remaining
    // non-measured qudits
    void mark_as_measured_(idx i) {
        if (was_measured(i))
            throw exception::QuditAlreadyMeasured(
                "qpp::QCircuit::mark_as_measured_()");
        subsys_[i] = idx_infty; // set qudit i to measured state
        for (idx m = i; m < qcd_.nq_; ++m) {
            if (!was_measured(m)) {
                --subsys_[m];
            }
        }
    }

    // giving a vector of non-measured qudits, get their relative position wrt
    // the measured qudits
    std::vector<idx> get_relative_pos_(std::vector<idx> v) {
        idx vsize = v.size();
        for (idx i = 0; i < vsize; ++i) {
            if (was_measured(v[i]))
                throw exception::QuditAlreadyMeasured(
                    "qpp::QCircuit::get_relative_pos_()");
            v[i] = subsys_[v[i]];
        }
        return v;
    }

  public:
    QCircuit(const QCircuitDescription& qcd)
        : qcd_{qcd}, psi_{st.zero(qcd.nq_, qcd.d_)}, dits_(qcd.nc_, 0),
          probs_(qcd.nc_, 0), subsys_(qcd.nq_, 0), m_ip_{0}, q_ip_{0}, ip_{0} {
        std::iota(std::begin(subsys_), std::end(subsys_), 0);
    }

    // getters
    // QCircuitDescription get_circuit_description() const { return qcd_; }

    ket get_psi() const { return psi_; }

    std::vector<idx> get_dits() const { return dits_; }

    std::vector<double> get_probs() const { return probs_; }

    idx was_measured(idx i) const { return subsys_[i] == idx_infty; }

    std::vector<idx> get_measured() const {
        std::vector<idx> result;
        for (idx i = 0; i < qcd_.nq_; ++i)
            if (was_measured(i))
                result.emplace_back(i);

        return result;
    }

    // qudits that were not (yet) measured
    std::vector<idx> get_not_measured() const {
        std::vector<idx> result;
        for (idx i = 0; i < qcd_.nq_; ++i)
            if (!was_measured(i))
                result.emplace_back(i);

        return result;
    }

    // measurement instruction pointer
    idx get_m_ip() const { return m_ip_; }

    // quantum instruction pointer
    idx get_q_ip() const { return q_ip_; }

    // combined instruction pointer
    idx get_ip() const { return ip_; }

    const QCircuitDescription& get_circuit_description() const { return qcd_; }
    // end getters

    void reset() {
        psi_ = st.zero(qcd_.nq_, qcd_.d_);
        dits_ = std::vector<idx>(qcd_.nc_, 0);
        probs_ = std::vector<double>(qcd_.nc_, 0);
        std::iota(std::begin(subsys_), std::end(subsys_), 0);
        m_ip_ = 0;
        q_ip_ = 0;
        ip_ = 0;
    }

    void run(idx step = idx_infty, bool verbose = false) {
        idx no_steps;
        if (step == idx_infty)
            no_steps = qcd_.get_total_count() - ip_;
        else
            no_steps = step;
        idx saved_ip_ = ip_; // save the instruction pointer at entry

        if (ip_ + no_steps > qcd_.get_total_count())
            throw exception::OutOfRange("qpp::QCircuit::run()");

        for (; q_ip_ <= qcd_.gates_.size(); ++q_ip_) {
            // check for measurements
            if (m_ip_ < qcd_.measurements_.size()) {
                idx m_step = qcd_.measurement_steps_[m_ip_];
                // we have a measurement at step i
                if (m_step == q_ip_) {
                    while (qcd_.measurement_steps_[m_ip_] == m_step) {
                        if (ip_ - saved_ip_ == no_steps)
                            return;

                        if (verbose) {
                            std::cout << std::left;
                            std::cout << std::setw(8)
                                      << qcd_.measurements_[m_ip_].step_no_;
                            std::cout << std::right;
                            std::cout << "|> " << qcd_.measurements_[m_ip_]
                                      << '\n';
                        }

                        std::vector<idx> target_rel_pos = get_relative_pos_(
                            qcd_.measurements_[m_ip_].target_);

                        std::vector<idx> resZ;
                        double probZ;

                        idx mres = 0;
                        std::vector<double> probs;
                        std::vector<cmat> states;

                        switch (qcd_.measurements_[m_ip_].measurement_type_) {
                        case QCircuitDescription::MeasureType::NONE:
                            break;
                        case QCircuitDescription::MeasureType::MEASURE_Z:
                            std::tie(resZ, probZ, psi_) =
                                measure_seq(psi_, target_rel_pos, qcd_.d_);
                            dits_[qcd_.measurements_[m_ip_].c_reg_] = resZ[0];
                            probs_[qcd_.measurements_[m_ip_].c_reg_] = probZ;
                            mark_as_measured_(
                                qcd_.measurements_[m_ip_].target_[0]);
                            break;
                        case QCircuitDescription::MeasureType::MEASURE_V:
                            std::tie(mres, probs, states) = measure(
                                psi_, qcd_.measurements_[m_ip_].mats_[0],
                                target_rel_pos, qcd_.d_);
                            psi_ = states[mres];
                            dits_[qcd_.measurements_[m_ip_].c_reg_] = mres;
                            probs_[qcd_.measurements_[m_ip_].c_reg_] =
                                probs[mres];
                            mark_as_measured_(
                                qcd_.measurements_[m_ip_].target_[0]);
                            break;
                        case QCircuitDescription::MeasureType::MEASURE_V_MANY:
                            std::tie(mres, probs, states) = measure(
                                psi_, qcd_.measurements_[m_ip_].mats_[0],
                                target_rel_pos, qcd_.d_);
                            psi_ = states[mres];
                            dits_[qcd_.measurements_[m_ip_].c_reg_] = mres;
                            probs_[qcd_.measurements_[m_ip_].c_reg_] =
                                probs[mres];
                            for (auto&& i : qcd_.measurements_[m_ip_].target_)
                                mark_as_measured_(i);
                            break;
                        }
                        ++ip_;
                        if (++m_ip_ == qcd_.measurements_.size())
                            break;
                    }
                }
            }

            // check for gates
            if (q_ip_ < qcd_.gates_.size()) {
                if (ip_ - saved_ip_ == no_steps)
                    return;
                if (verbose) {
                    std::cout << std::left;
                    std::cout << std::setw(8) << qcd_.gates_[q_ip_].step_no_;
                    std::cout << std::right;
                    std::cout << qcd_.gates_[q_ip_] << '\n';
                }

                std::vector<idx> ctrl_rel_pos;
                std::vector<idx> target_rel_pos =
                    get_relative_pos_(qcd_.gates_[q_ip_].target_);

                switch (qcd_.gates_[q_ip_].gate_type_) {
                case QCircuitDescription::GateType::NONE:
                    break;
                case QCircuitDescription::GateType::SINGLE:
                case QCircuitDescription::GateType::TWO:
                case QCircuitDescription::GateType::THREE:
                case QCircuitDescription::GateType::CUSTOM:
                    psi_ = apply(psi_, qcd_.gates_[q_ip_].gate_, target_rel_pos,
                                 qcd_.d_);
                    break;
                case QCircuitDescription::GateType::FAN:
                    for (idx m = 0; m < qcd_.gates_[q_ip_].target_.size(); ++m)
                        psi_ = apply(psi_, qcd_.gates_[q_ip_].gate_,
                                     {target_rel_pos[m]}, qcd_.d_);
                    break;
                case QCircuitDescription::GateType::QFT:
                case QCircuitDescription::GateType::TFQ:
                case QCircuitDescription::GateType::SINGLE_CTRL_SINGLE_TARGET:
                case QCircuitDescription::GateType::SINGLE_CTRL_MULTIPLE_TARGET:
                case QCircuitDescription::GateType::MULTIPLE_CTRL_SINGLE_TARGET:
                case QCircuitDescription::GateType::
                    MULTIPLE_CTRL_MULTIPLE_TARGET:
                case QCircuitDescription::GateType::CUSTOM_CTRL:
                    ctrl_rel_pos = get_relative_pos_(qcd_.gates_[q_ip_].ctrl_);
                    psi_ = applyCTRL(psi_, qcd_.gates_[q_ip_].gate_,
                                     ctrl_rel_pos, target_rel_pos, qcd_.d_);
                    break;
                case QCircuitDescription::GateType::SINGLE_cCTRL_SINGLE_TARGET:
                case QCircuitDescription::GateType::
                    SINGLE_cCTRL_MULTIPLE_TARGET:
                case QCircuitDescription::GateType::
                    MULTIPLE_cCTRL_SINGLE_TARGET:
                case QCircuitDescription::GateType::
                    MULTIPLE_cCTRL_MULTIPLE_TARGET:
                case QCircuitDescription::GateType::CUSTOM_cCTRL:
                    if (dits_.size() == 0) {
                        psi_ = apply(psi_, qcd_.gates_[q_ip_].gate_,
                                     target_rel_pos, qcd_.d_);
                    } else {
                        bool should_apply = true;
                        idx first_dit = dits_[(qcd_.gates_[q_ip_].ctrl_)[0]];
                        for (idx m = 0; m < qcd_.gates_[q_ip_].ctrl_.size();
                             ++m) {
                            if (dits_[(qcd_.gates_[q_ip_].ctrl_)[m]] !=
                                first_dit) {
                                should_apply = false;
                                break;
                            }
                        }
                        if (should_apply) {
                            psi_ = apply(
                                psi_, powm(qcd_.gates_[q_ip_].gate_, first_dit),
                                target_rel_pos, qcd_.d_);
                        }
                    }
                    break;
                }
                ++ip_;
            }
        }
        --q_ip_;
    }

    std::ostream& display(std::ostream& os) const override {
        os << "measured: " << disp(get_measured(), ",") << '\n';
        os << "dits: " << disp(dits_, ",") << '\n';
        os << "probs: " << disp(probs_, ",") << '\n';
        //        os << "updated subsystems: ";
        //        std::cout << disp(subsys_, ",") << '\n';

        return os;
    }
};

} // namespace experimental
} /* namespace qpp */

#endif /* EXPERIMENTAL_EXPERIMENTAL_H_ */
