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
 * \file qpp/classes/qengine.hpp
 * \brief Quantum engines
 */

#ifndef QPP_CLASSES_QENGINE_HPP_
#define QPP_CLASSES_QENGINE_HPP_

// #define DEBUGGING
#include <iostream>
// clang-format off
#ifdef DEBUGGING
#define LOG std::cout
#else
#define LOG if (false) std::cout
#endif
// clang-format on

#include <algorithm>
#include <iterator>
#include <limits>
#include <map>
#include <optional>
#include <ostream>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

#include "qpp/input_output.hpp"
#include "qpp/instruments.hpp"
#include "qpp/operations.hpp"
#include "qpp/types.hpp"

#include "qpp/classes/qbase_engine.hpp"
#include "qpp/classes/qcircuit.hpp"

#include "qpp/internal/classes/qcircuit_gate_step.hpp"
#include "qpp/internal/classes/qcircuit_measurement_step.hpp"
#include "qpp/internal/classes/qengine_state.hpp"
#include "qpp/internal/classes/qengine_statistics.hpp"

namespace qpp {
/**
 * \class qpp::QEngineT
 * \brief Quantum engine, executes qpp::QCircuit
 * \see qpp::QNoisyEngineT, qpp::QCircuit
 *
 * \note When interfacing with OpenQASM, the classical dits/registers are
 * evaluated in little-endian order, with  the least significant bit being
 * stored first. For example, [1,0,0] is interpreted as 1 (and not 4).
 * See https://github.com/softwareQinc/qpp/issues/75 for more details.
 *
 * \tparam T Engine's state underlying type, qpp::ket or qpp::cmat
 */
template <typename T>
class QEngineT : public QBaseEngine<T, QCircuit> {
    /**
     * \brief Executes a contiguous series of projective measurement steps
     *
     * \note When \a reps > 1, the last repetition is executed repeatedly until
     * the post-selection succeeds, so one can retrieve the quantum state
     *
     * \param engine_state Instance of qpp::internal::QEngineState<T>
     * \param steps Vector of qpp::QCircuit::iterator
     * \idx pos Index from where the execution starts
     * \idx reps Number of repetitions
     * \return Reference to the current instance
     */
    // IMPORTANT: ALWAYS pass engine_state by value!
    QEngineT& execute_no_sample_(internal::QEngineState<T> engine_state,
                                 const std::vector<QCircuit::iterator>& steps,
                                 idx pos, idx reps) {
        auto worker = [&](bool last_rep) {
            // sets the state of the engine to the entry state
            qeng_st_ = engine_state;
            if (last_rep) {
                qeng_st_.ensure_post_selection_ = true;
            }

            bool measured = false;
            for (idx i = pos; i < steps.size(); ++i) {
                if (internal::is_measurement(steps[i])) {
                    measured = true;
                }
                this->execute(steps[i]);
                // post-selection failed, stop executing
                if (!qeng_st_.post_select_ok_) {
                    return;
                }
            }

            // at least one qudit was measured
            if (measured) {
                std::vector<idx> m_res = get_dits();
                ++stats_.data()[m_res];
            }

            // restores qeng_st_.ensure_post_selection_ flag
            if (last_rep) {
                qeng_st_.ensure_post_selection_ =
                    engine_state.ensure_post_selection_;
            }
        };

        for (idx rep = 0; rep < reps - 1; ++rep) {
            worker(false);
        }
        if (reps > 1) {
            worker(true); // ensure post-selection for >1 reps
        } else {
            worker(qeng_st_.ensure_post_selection_);
        }

        return *this;
    }

    /**
     * \brief Compute the position of the first measurement/discard/reset step
     * \see qpp::internal::canonical_form()
     *
     * \param steps Vector of qpp::QCircuit::iterator
     * \return Position of the first measurement/discard/reset step
     */
    idx compute_first_measurement_discard_reset_pos_(
        const std::vector<QCircuit::iterator>& steps) const {
        // find the position of the first measurement/reset/discard step
        auto first_measurement_discard_reset_it =
            std::find_if(steps.begin(), steps.end(), [](auto&& elem) {
                return internal::is_measurement(elem) ||
                       internal::is_discard(elem) || internal::is_reset(elem);
            });
        idx first_measurement_discard_reset_pos =
            std::distance(steps.begin(), first_measurement_discard_reset_it);

        return first_measurement_discard_reset_pos;
    }

    /**
     * \brief Returns pair of (bool, idx), first true if the canonical form of
     * the circuit can be sampled from, second denoting the position of the
     * first measurement/reset/discard step \see qpp::internal::canonical_form()
     *
     * \note A circuit can be sampled from if and only if it its canonical form
     * have only projective measurements (including post-selection) after the
     * first measurement/reset/discard step. In other words, it's of the form
     * [Gate step(s)] * [Projective measurements step(s)].
     *
     * \param steps Vector of qpp::QCircuit::iterator
     * \return Pair of (bool, idx), first true if the canonical form of the
     * circuit can be sampled from, second denoting the position of the first
     * measurement/reset/discard step
     */
    std::pair<bool, idx>
    can_sample_from_(const std::vector<QCircuit::iterator>& steps) const {

        // in the following, we will partition the circuit as
        // [0 ... first_measurement_discard_reset_pos ... end)

        // find the position of the first measurement/reset/discard step
        idx first_measurement_discard_reset_pos =
            compute_first_measurement_discard_reset_pos_(steps);

        // decide if we can sample (every step after
        // first_measurement_discard_reset_pos must be a projective measurement)
        for (idx i = first_measurement_discard_reset_pos; i < steps.size();
             ++i) {
            if (!(internal::is_projective_measurement(steps[i])) ||
                internal::is_discard(steps[i])) {
                return {false, first_measurement_discard_reset_pos};
            }
        }

        return {true, first_measurement_discard_reset_pos};
    }

    /**
     * \brief Builds internal maps used by qpp::QCircuit::execute_sample_()
     *
     * \param steps Vector of qpp::QCircuit::iterator
     * \param pos Index from where the execution starts; from this index
     * further,
     * all steps are assumed to be a projective measurement (including
     * post-selection)
     * \param[out] c_q Records the c <- (q, [OPTIONAL ps_val]) map
     * \param[out] q_ps_val Records the q <- ps_val map
     * \param[out] q_m Qudits that were measured projectively, including
     * post-selected
     * \return Pair of (bool, bool), first representing whether at least one
     * qudit was measured, and second representing whether post-selected values
     * are compatible (e.g., post-select 0 followed by post-select 1 is
     * incompatible)
     */
    std::pair<bool, bool>
    build_sampling_maps_(const std::vector<QCircuit::iterator>& steps, idx pos,
                         std::map<idx, std::pair<idx, std::optional<idx>>>& c_q,
                         std::map<idx, idx>& q_ps_val, std::set<idx>& q_m) {
        bool measured = false;

        for (idx i = pos; i < steps.size(); ++i) {
            // post-selection
            if (internal::is_projective_post_selection(steps[i])) {
                measured = true;
                auto [_, target, ps_vals, c_regs] =
                    internal::extract_ctrl_target_ps_vals_c_reg(steps[i]);
                for (idx q = 0; q < static_cast<idx>(target.size()); ++q) {
                    q_m.emplace(target[q]);
                    c_q[c_regs[q]] = {target[q], ps_vals[q]};
                    // build the q <- ps_val map
                    // insert if not exists
                    if (auto it = q_ps_val.find(target[q]);
                        it == q_ps_val.end()) {
                        q_ps_val[target[q]] = ps_vals[q];
                    } else
                    // element already exists, make sure that
                    // post-selection is consistent
                    {
                        // post-selection is incompatible
                        if (auto ps_val = it->second; ps_val != ps_vals[q]) {
                            qeng_st_.post_select_ok_ = false;
                            LOG << "fails here" << std::endl;
                            LOG << "sample? " << qeng_st_.can_sample_
                                << std::endl;
                            return {measured, false};
                        }
                    }
                }
                LOG << "1" << std::endl;
            }
            // projective measurement
            else if (internal::is_projective_measurement(steps[i])) {
                measured = true;
                auto [_, target, ps_vals_, c_regs] =
                    internal::extract_ctrl_target_ps_vals_c_reg(steps[i]);
                for (idx q = 0; q < static_cast<idx>(target.size()); ++q) {
                    q_m.emplace(target[q]);
                    c_q[c_regs[q]] = {target[q], {}};
                }
                LOG << "2" << std::endl;
            }
        }

        return {measured, true};
    }

    /**
     * \brief Executes a contiguous series of projective measurement steps via
     * sampling
     *
     * \param engine_state Instance of qpp::internal::QEngineState<T>
     * \param steps Vector of qpp::QCircuit::iterator
     * \param pos Index from where the execution starts; from this index
     * further, all steps are assumed to be a projective measurement (including
     * post-selection)
     * \param reps Number of repetitions
     * \return Reference to the current instance
     */
    QEngineT& execute_sample_(const internal::QEngineState<T>& engine_state,
                              const std::vector<QCircuit::iterator>& steps,
                              idx pos, idx reps) {
        std::map<idx, std::pair<idx, std::optional<idx>>>
            c_q; // records the c <- (q, [OPTIONAL ps_val]) map
        std::map<idx, idx> q_ps_val; // records the q <- ps_val map
        std::set<idx> q_m;           // qudits that were measured projectively,
                                     // including post-selected

        auto [measured, post_select_compatible] =
            build_sampling_maps_(steps, pos, c_q, q_ps_val, q_m);

        if (!post_select_compatible) {
            qeng_st_.post_select_ok_ = false;
            return *this;
        }

        if (!measured) {
            return *this;
        }

        // display sampling maps
        LOG << "\nSAMPLING MAPS etc.\n\n";
        LOG << "c_q:\n";
        for (auto elem : c_q) {
            LOG << elem.first << " - " << elem.second.first << " "
                << (elem.second.second.has_value()
                        ? std::to_string(elem.second.second.value())
                        : "")
                << std::endl;
        }
        LOG << std::endl;

        LOG << "q_ps_val:\n";
        for (auto elem : q_ps_val) {
            LOG << elem.first << " - " << elem.second << std::endl;
        }
        LOG << std::endl;

        LOG << "q_m:" << std::endl;
        for (auto elem : q_m) {
            LOG << elem << " ";
        }
        LOG << std::endl << std::endl;

        // at least one qudit was measured, add it to stats_

        // build the vector of measured qudits that we must sample
        // from
        LOG << "SAMPLE FROM SET 1:" << std::endl;
        std::set<idx> sample_from_set;
        for (auto [dit, qudit_ps_val] : c_q) {
            sample_from_set.emplace(qudit_ps_val.first);
            LOG << qudit_ps_val.first << std::endl;
        }
        LOG << std::endl;

        std::vector<idx> sample_from(sample_from_set.begin(),
                                     sample_from_set.end());

        LOG << "CORRECT 'TILL HERE\n" << std::endl;

        for (idx rep = 0; rep < reps - 1; ++rep) {
        REPEAT_POST_SELECT_SAMPLE:
            // sample from the quantum state
            std::vector<idx> sample_result_restricted_support = sample(
                engine_state.qstate_, sample_from, this->qc_ptr_->get_d());

            // make sure post-selected qudits agree on their values
            bool post_selection_failed = false;
            for (idx q = 0; q < sample_from.size(); ++q) {
                // qudit sample_from[i] is a post-selected qudit
                if (auto it = q_ps_val.find(sample_from[q]);
                    it != q_ps_val.end()) {
                    if (sample_result_restricted_support[q] != it->second) {
                        post_selection_failed = true;
                        break;
                    }
                }
            }

            if (post_selection_failed) {
                if (qeng_st_.ensure_post_selection_) {
                    goto REPEAT_POST_SELECT_SAMPLE;
                } else {
                    continue;
                }
            }

            // extend sample_result to full support
            std::vector<idx> sample_result = this->get_dits();
            idx i = 0;
            for (auto [dit, qudit_ps_val] : c_q) {
                sample_result[dit] = sample_result_restricted_support[i++];
            }

            ++stats_.data()[sample_result];
        }

        // the last state is always computed with ensure_post_selection = true
        execute_no_sample_(qeng_st_, steps, pos, 1);

        return *this;
    }

  protected:
    internal::QEngineState<T> qeng_st_; ///< current state of the engine
    internal::QEngineStatistics
        stats_{}; ///< measurement statistics for multiple runs

    /**
     * \brief Marks qudit \a i as measured then re-label accordingly the
     * remaining non-measured qudits
     *
     * \param i Qudit index
     */
    void set_measured_destructively_(idx i) {
        // EXCEPTION CHECKS
        if (was_measured_destructively(i)) {
            throw exception::QuditAlreadyMeasured(
                "qpp::QEngineT::set_measured_()", "i");
        }
        // END EXCEPTION CHECKS
        qeng_st_.subsys_[i] =
            std::numeric_limits<idx>::max(); // set qudit i to measured state
        for (idx m = i; m < this->qc_ptr_->get_nq(); ++m) {
            if (!was_measured_destructively(m)) {
                --qeng_st_.subsys_[m];
            }
        }
    }

    // giving a vector of non-measured qudits, get their relative position
    // w.r.t. the measured qudits
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
            if (was_measured_destructively(v[i])) {
                throw exception::QuditAlreadyMeasured(
                    "qpp::QEngineT::get_relative_pos_()", "v[i]");
            }
            // END EXCEPTION CHECKS
            v[i] = qeng_st_.subsys_[v[i]];
        }
        return v;
    }

    /**
     * \brief Executes qpp::internal::QCircuitGateStep
     *
     * \param gate_step Instance of qpp::internal::QCircuitGateStep
     */
    virtual void
    execute_gate_step_(const internal::QCircuitGateStep& gate_step) {
        std::vector<idx> ctrl_rel_pos;
        std::vector<idx> target_rel_pos = get_relative_pos_(gate_step.target_);

        auto h_tbl = this->qc_ptr_->get_cmat_hash_tbl();
        idx d = this->qc_ptr_->get_d();

        using GType = internal::QCircuitGateStep::Type;
        GType gstep_type_ = gate_step.gate_type_;

        // regular gate
        switch (gstep_type_) {
            case GType::NONE:
                break;
            case GType::SINGLE:
            case GType::TWO:
            case GType::THREE:
            case GType::JOINT:
                qeng_st_.qstate_ =
                    apply(qeng_st_.qstate_, h_tbl[gate_step.gate_hash_],
                          target_rel_pos, d);
                break;
            case GType::FAN:
                for (idx m = 0; m < static_cast<idx>(gate_step.target_.size());
                     ++m) {
                    qeng_st_.qstate_ =
                        apply(qeng_st_.qstate_, h_tbl[gate_step.gate_hash_],
                              {target_rel_pos[m]}, d);
                }
                break;
            default:
                break;
        }

        // controlled gate
        if (QCircuit::is_CTRL(gate_step)) {
            ctrl_rel_pos = get_relative_pos_(gate_step.ctrl_.value());
            bool is_fan = (gate_step.gate_type_ == GType::CTRL_FAN);
            qeng_st_.qstate_ =
                is_fan
                    ? applyCTRL_fan(qeng_st_.qstate_,
                                    h_tbl[gate_step.gate_hash_], ctrl_rel_pos,
                                    target_rel_pos, d, gate_step.shift_)
                    : applyCTRL(qeng_st_.qstate_, h_tbl[gate_step.gate_hash_],
                                ctrl_rel_pos, target_rel_pos, d,
                                gate_step.shift_);
        }

        // classically-controlled gate
        if (QCircuit::is_cCTRL(gate_step)) {
            bool is_fan = (gate_step.gate_type_ == GType::cCTRL_FAN);
            if (!qeng_st_.dits_.empty()) {
                {
                    bool should_apply = true;
                    idx first_dit;
                    // we have a shift
                    if (gate_step.shift_.has_value()) {
                        first_dit =
                            (qeng_st_.dits_[(gate_step.ctrl_.value())[0]] +
                             gate_step.shift_.value()[0]) %
                            d;
                        for (idx m = 1; m < static_cast<idx>(
                                                gate_step.ctrl_.value().size());
                             ++m) {
                            if ((qeng_st_.dits_[(gate_step.ctrl_.value())[m]] +
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
                            qeng_st_.dits_[(gate_step.ctrl_.value())[0]];
                        for (idx m = 1; m < static_cast<idx>(
                                                gate_step.ctrl_.value().size());
                             ++m) {
                            if (qeng_st_.dits_[(gate_step.ctrl_.value())[m]] !=
                                first_dit) {
                                should_apply = false;
                                break;
                            }
                        }
                    }
                    if (should_apply) {
                        cmat U = powm(h_tbl[gate_step.gate_hash_], first_dit);
                        if (is_fan) {
                            for (idx qudit : target_rel_pos) {
                                qeng_st_.qstate_ =
                                    apply(qeng_st_.qstate_, U, {qudit}, d);
                            }
                        } else {
                            qeng_st_.qstate_ = apply(
                                qeng_st_.qstate_,
                                powm(h_tbl[gate_step.gate_hash_], first_dit),
                                target_rel_pos, d);
                        }
                    }
                }
            }
            // TODO check if this can happen (st_.dits_.empty())
            else {
                if (is_fan) {
                    for (idx qudit : target_rel_pos) {
                        qeng_st_.qstate_ =
                            apply(qeng_st_.qstate_, h_tbl[gate_step.gate_hash_],
                                  {qudit}, d);
                    }
                } else {
                    qeng_st_.qstate_ =
                        apply(qeng_st_.qstate_, h_tbl[gate_step.gate_hash_],
                              target_rel_pos, d);
                }
            }
        } // end if classically-controlled gate
    }

    /**
     * \brief Executes qpp::internal::QCircuitMeasurementStep
     *
     * \param measurement_step Instance of
     * qpp::internal::QCircuitMeasurementStep
     */
    virtual void execute_measurement_step_(
        const internal::QCircuitMeasurementStep& measurement_step) {
        std::vector<idx> target_rel_pos =
            get_relative_pos_(measurement_step.target_);

        idx mres = 0;
        std::vector<idx> results;
        std::vector<realT> probs;
        std::vector<T> states;
        T qstate;

        auto h_tbl = this->qc_ptr_->get_cmat_hash_tbl();
        idx d = this->qc_ptr_->get_d();

        bool destructive =
            internal::is_destructive_measurement(measurement_step);

        using MType = internal::QCircuitMeasurementStep::Type;
        MType mstep_type_ = measurement_step.measurement_type_;

        switch (mstep_type_) {
            case MType::NONE:
                break;

            case MType::MEASURE:
            case MType::MEASURE_ND:
            case MType::POST_SELECT:
            case MType::POST_SELECT_ND:
            REPEAT_POST_SELECT:
                std::tie(results, probs, qstate) = measure_seq(
                    qeng_st_.qstate_, target_rel_pos, d, destructive);

                // POST_SELECT(_ND)
                if (mstep_type_ == MType::POST_SELECT ||
                    mstep_type_ == MType::POST_SELECT_ND) {
                    // post-selection failed
                    if (results[0] != measurement_step.ps_vals_.value()[0]) {
                        if (!qeng_st_.ensure_post_selection_) {
                            qeng_st_.post_select_ok_ = false;
                        } else {
                            goto REPEAT_POST_SELECT;
                        }
                        // std::cerr << "POST_SELECT failed!" << std::endl;
                    }
                }

                qeng_st_.qstate_ = std::move(qstate);
                qeng_st_.dits_[measurement_step.c_reg_] = results[0];
                qeng_st_.probs_[measurement_step.c_reg_] = probs[0];
                if (destructive) {
                    set_measured_destructively_(measurement_step.target_[0]);
                }
                break;

            case MType::MEASURE_MANY:
            case MType::MEASURE_MANY_ND:
            case MType::POST_SELECT_MANY:
            case MType::POST_SELECT_MANY_ND:
            REPEAT_POST_SELECT_MANY:
                std::tie(results, probs, qstate) = measure_seq(
                    qeng_st_.qstate_, target_rel_pos, d, destructive);

                // POST_SELECT_MANY(_ND)
                if (mstep_type_ == MType::POST_SELECT_MANY ||
                    mstep_type_ == MType::POST_SELECT_MANY_ND) {
                    // post-selection failed
                    if (results != measurement_step.ps_vals_.value()) {
                        if (!qeng_st_.ensure_post_selection_) {
                            qeng_st_.post_select_ok_ = false;
                        } else {
                            goto REPEAT_POST_SELECT_MANY;
                        }
                    }
                }

                qeng_st_.qstate_ = std::move(qstate);
                std::copy(
                    results.begin(), results.end(),
                    std::next(qeng_st_.dits_.begin(), measurement_step.c_reg_));
                std::copy(probs.begin(), probs.end(),
                          std::next(qeng_st_.probs_.begin(),
                                    measurement_step.c_reg_));
                if (destructive) {
                    for (idx target : measurement_step.target_) {
                        set_measured_destructively_(target);
                    }
                }
                break;

            case MType::MEASURE_V:
            case MType::MEASURE_V_ND:
            case MType::POST_SELECT_V:
            case MType::POST_SELECT_V_ND:
            REPEAT_POST_SELECT_V:
                std::tie(mres, probs, states) = measure(
                    qeng_st_.qstate_, h_tbl[measurement_step.mats_hash_[0]],
                    target_rel_pos, d, destructive);

                // POST_SELECT_V(_ND)
                if (mstep_type_ == MType::POST_SELECT_V ||
                    mstep_type_ == MType::POST_SELECT_V_ND) {
                    // post-selection failed
                    if (mres != measurement_step.ps_vals_.value()[0]) {
                        if (!qeng_st_.ensure_post_selection_) {
                            qeng_st_.post_select_ok_ = false;
                        } else {
                            goto REPEAT_POST_SELECT_V;
                        }
                    }
                }

                qeng_st_.qstate_ = states[mres];
                qeng_st_.dits_[measurement_step.c_reg_] = mres;
                qeng_st_.probs_[measurement_step.c_reg_] = probs[mres];
                if (destructive) {
                    set_measured_destructively_(measurement_step.target_[0]);
                }
                break;

            case MType::MEASURE_V_JOINT:
            case MType::MEASURE_V_JOINT_ND:
            case MType::POST_SELECT_V_JOINT:
            case MType::POST_SELECT_V_JOINT_ND:
            REPEAT_POST_SELECT_V_JOINT:
                std::tie(mres, probs, states) = measure(
                    qeng_st_.qstate_, h_tbl[measurement_step.mats_hash_[0]],
                    target_rel_pos, d, destructive);

                // POST_SELECT_V(_ND)
                if (mstep_type_ == MType::POST_SELECT_V_JOINT ||
                    mstep_type_ == MType::POST_SELECT_V_JOINT_ND) {
                    // post-selection failed
                    if (mres != measurement_step.ps_vals_.value()[0]) {
                        if (!qeng_st_.ensure_post_selection_) {
                            qeng_st_.post_select_ok_ = false;
                        } else {
                            goto REPEAT_POST_SELECT_V_JOINT;
                        }
                        // std::cerr << "POST_SELECT_V_JOINT failed!" <<
                        // std::endl;
                    }
                }

                qeng_st_.qstate_ = states[mres];
                qeng_st_.dits_[measurement_step.c_reg_] = mres;
                qeng_st_.probs_[measurement_step.c_reg_] = probs[mres];
                if (destructive) {
                    for (idx target : measurement_step.target_) {
                        set_measured_destructively_(target);
                    }
                }
                break;

            case MType::RESET:
            case MType::RESET_MANY:
                qeng_st_.qstate_ =
                    qpp::reset(qeng_st_.qstate_, target_rel_pos, d);
                break;

            case MType::DISCARD:
                std::tie(std::ignore, std::ignore, qeng_st_.qstate_) =
                    measure_seq(qeng_st_.qstate_, target_rel_pos, d);
                set_measured_destructively_(measurement_step.target_[0]);
                break;
            case MType::DISCARD_MANY:
                std::tie(std::ignore, std::ignore, qeng_st_.qstate_) =
                    measure_seq(qeng_st_.qstate_, target_rel_pos, d);
                for (idx target : measurement_step.target_) {
                    set_measured_destructively_(target);
                }
                break;
        } // end switch on measurement type
    };

    /**
     * \brief Executes qpp::internal::QCircuitNOPStep
     *
     * \param nop_step Instance of qpp::internal::QCircuitNOPStep
     */
    virtual void execute_nop_step_(
        [[maybe_unused]] const internal::QCircuitNOPStep& nop_step) {}

  public:
    using QBaseEngine<T, QCircuit>::QBaseEngine;
    using QBaseEngine<T, QCircuit>::execute;

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
     * \param ensure_post_selection If true, repeatedly executes post-selection
     * steps until the post-selection result(s) agree, false by default
     */
    explicit QEngineT(const QCircuit& qc, bool ensure_post_selection = false)
        : QBaseEngine<T, QCircuit>{qc}, qeng_st_{this->qc_ptr_}, stats_{} {
        qeng_st_.ensure_post_selection_ = ensure_post_selection;
    }

    // traits
    /**
     * \brief qpp::IQEngineTraits::traits_get_name() override
     */
    virtual std::string traits_get_name() const override { return "QEngineT"; }

    /**
     * \brief qpp::IQEngineTraits::traits_is_noisy() override
     */
    virtual bool traits_is_noisy() const override { return false; }

    /**
     * \brief qpp::IQEngineTraits::traits_is_pure() override
     */
    virtual bool traits_is_pure() const override {
        if constexpr (std::is_same_v<T, cmat>) {
            return false;
        } else if constexpr (std::is_same_v<T, ket>) {
            return true;
        }
        // default, we assume pure states
        return true;
    }
    // end traits

    // getters
    /**
     * \brief Underlying quantum state
     *
     * \return Underlying quantum state
     */
    T get_state() const override { return qeng_st_.qstate_; }

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
    std::vector<idx> get_dits() const { return qeng_st_.dits_; }

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
        if (i >= this->qc_ptr_->get_nc()) {
            throw exception::OutOfRange("qpp::QEngineT::get_dit()", "i");
        }
        // END EXCEPTION CHECKS

        return qeng_st_.dits_[i];
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
    std::vector<realT> get_probs() const { return qeng_st_.probs_; }

    /**
     * \brief Check whether qudit \a i was already measured (destructively)
     * w.r.t. the current engine state
     *
     * \param i Qudit index
     * \return True if qudit \a i was already measured, false otherwise
     */
    bool was_measured_destructively(idx i) const {
        return qeng_st_.subsys_[i] == std::numeric_limits<idx>::max();
    }

    /**
     * \brief Vector of already measured (destructively) qudit indexes at
     * the current engine state
     *
     * \return Vector of already measured qudit (destructively) indexes at
     * the current engine state
     */
    std::vector<idx> get_measured_destructively() const {
        std::vector<idx> result;
        for (idx i = 0; i < this->qc_ptr_->get_nq(); ++i) {
            if (was_measured_destructively(i)) {
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
    std::vector<idx> get_non_measured_destructively() const {
        std::vector<idx> result;
        for (idx i = 0; i < this->qc_ptr_->get_nq(); ++i) {
            if (!was_measured_destructively(i)) {
                result.emplace_back(i);
            }
        }

        return result;
    }

    /**
     * \brief Measurement statistics for multiple runs
     *
     * \return Hash table with collected measurement statistics for multiple
     * runs, with hash key being the string representation of the vector of
     * measurement results and value being the number of occurrences (of the
     * vector of measurement results), with the most significant bit located
     * at index 0 (i.e., top/left).
     */
    internal::QEngineStatistics get_stats() const { return stats_; }

    /**
     * \brief Returns true if post-selection was successful (or absent), false
     * otherwise
     */
    bool post_select_ok() const { return qeng_st_.post_select_ok_; }

    /**
     * \brief Returns true if post-selection is ensured to succeed, false
     * otherwise
     */
    bool ensure_post_selection() const {
        return qeng_st_.ensure_post_selection_;
    }
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
    QEngineT& set_dit(idx i, idx value) {
        // EXCEPTION CHECKS
        if (i >= this->qc_ptr_->get_nc()) {
            throw exception::OutOfRange("qpp::QEngineT::set_dit()", "i");
        }
        // END EXCEPTION CHECKS
        qeng_st_.dits_[i] = value;

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
    QEngineT& set_dits(std::vector<idx> dits) {
        // EXCEPTION CHECKS
        if (dits.size() != qeng_st_.dits_.size()) {
            throw exception::SizeMismatch("qpp::QEngineT::set_dits()", "dits");
        }
        // END EXCEPTION CHECKS
        qeng_st_.dits_ = std::move(dits);

        return *this;
    }

    /**
     * \brief Sets the underlying quantum state to \a state
     *
     * \note The order is lexicographical with respect to the remaining
     * non-measured qudits
     *
     * \param[out] state State vector
     * \return Reference to the current instance
     */
    QEngineT& set_state(const T& state) override {
        // EXCEPTION CHECKS
        idx n = get_non_measured_destructively().size();
        idx D = internal::safe_pow(this->qc_ptr_->get_d(), n);
        if constexpr (std::is_same_v<T, ket>) {
            if (static_cast<idx>(state.rows()) != D) {
                throw exception::DimsNotEqual("qpp::QEngineT::set_state()",
                                              "state");
            }
        } else {
            if (!internal::check_square_mat(state)) {
                throw exception::MatrixNotSquare("qpp::QEngineT::set_state()");
            }
            if (static_cast<idx>(state.rows()) != D) {
                throw exception::DimsNotEqual("qpp::QEngineT::set_state()",
                                              "state");
            }
        }
        // END EXCEPTION CHECKS

        qeng_st_.qstate_ = state;

        return *this;
    }
    // end setters

    /**
     * \brief Resets the collected measurement statistics hash table
     *
     * \return Reference to the current instance
     */
    QEngineT& reset_stats() {
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
     * \param ensure_post_selection Optional (false by default). If true,
     * executes a measurement step repeatedly until the post-selection result(s)
     * agree
     * \return Reference to the current instance
     */
    virtual QEngineT& reset(bool reset_stats = true,
                            bool ensure_post_selection = false) {
        qeng_st_.reset(qeng_st_.qstate_);
        qeng_st_.ensure_post_selection_ = ensure_post_selection;

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
    QEngineT& execute(
        const typename QCircuitTraits<QCircuit>::value_type& elem) override {
        // EXCEPTION CHECKS
        // iterator must point to the same quantum circuit description
        if (elem.get_qc_ptr() != this->qc_ptr_) {
            throw exception::InvalidIterator(
                "qpp::QEngineT::execute()",
                "Iterator does not point to the same circuit description");
        }
        // the rest of exceptions are caught by the iterator::operator*()
        // END EXCEPTION CHECKS

        auto gate_step_visitor =
            [&](const internal::QCircuitGateStep& gate_step) {
                return this->execute_gate_step_(gate_step);
            };

        auto measurement_step_visitor =
            [&](const internal::QCircuitMeasurementStep& measurement_step) {
                return this->execute_measurement_step_(measurement_step);
            };
        auto nop_step_visitor = [&](const internal::QCircuitNOPStep& nop_step) {
            return this->execute_nop_step_(nop_step);
        };

        std::visit(
            overloaded{
                gate_step_visitor,
                measurement_step_visitor,
                nop_step_visitor,
            },
            elem.get_step());

        return *this;
    }

    /**
     * \brief Executes the entire quantum circuit description
     *
     * \param reps Number of repetitions
     * \return Reference to the current instance
     */
    QEngineT& execute(idx reps = 1) override {
        // EXCEPTION CHECKS
        if (reps == 0) {
            throw exception::OutOfRange("qpp::QEngine::execute()", "reps");
        }
        // END EXCEPTION CHECKS

        this->reset(false, qeng_st_.ensure_post_selection_);
        auto steps = reps > 1 ? internal::canonical_form(*this->qc_ptr_)
                              : internal::circuit_as_iterators(*this->qc_ptr_);
        if (steps.empty()) {
            return *this;
        }

        LOG << "Circuit in canonical form" << std::endl;
        for (auto elem : steps) {
            LOG << "---> " << *elem << std::endl;
        }

        auto [can_sample, first_measurement_discard_reset_pos] =
            can_sample_from_(steps);
        qeng_st_.can_sample_ = reps > 1 && can_sample;

        // execute everything ONCE in the interval
        // [0, first_measurement_discard_reset_pos)
        for (idx i = 0; i < first_measurement_discard_reset_pos; ++i) {
            execute(steps[i]);
        }

        // TODO: remove this after debugging
        // qeng_st_.can_sample_ = false;

        // execute repeatedly everything in the remaining interval
        // can sample: every step from now on is a projective measurement
        if (qeng_st_.can_sample_) {
            execute_sample_(qeng_st_, steps,
                            first_measurement_discard_reset_pos, reps);
        }
        // cannot sample
        else {
            execute_no_sample_(qeng_st_, steps,
                               first_measurement_discard_reset_pos, reps);
        }

        return *this;
    }

    /**
     * \brief qpp::IJSON::to_JSON() override
     *
     * Displays the state of the engine in JSON format
     *
     * \param enclosed_in_curly_brackets If true, encloses the result in
     * curly brackets
     * \return String containing the JSON representation of the state of the
     * engine
     */
    std::string to_JSON(bool enclosed_in_curly_brackets = true) const override {
        std::string result;

        if (enclosed_in_curly_brackets) {
            result += "{";
        }

        result += QBaseEngine<T, QCircuit>::to_JSON(false);

        std::ostringstream ss;
        ss << ", ";

        ss << "\"sampling\": " << (qeng_st_.can_sample_ ? "true" : "false")
           << ", ";
        ss << "\"measured/discarded (destructively)\": ";
        ss << disp(get_measured_destructively(),
                   IOManipContainerOpts{}.set_sep(", "));
        result += ss.str() + ", ";

        ss.str("");
        ss.clear();
        ss << disp(get_non_measured_destructively(),
                   IOManipContainerOpts{}.set_sep(", "));
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
     * Writes to the output stream a textual representation of the state of
     * the engine
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

        os << '[' << traits_get_name() << ']';
        if (qeng_st_.can_sample_) {
            os << " (Sampling)";
        }
        os << '\n';

        os << "<QCircuit nq: " << this->get_circuit().get_nq()
           << ", nc: " << this->get_circuit().get_nc()
           << ", d: " << this->get_circuit().get_d();

        if (this->get_circuit().get_name().has_value()) {
            os << ", name: ";
            os << "\"" << this->get_circuit().get_name().value() << "\"";
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
}; /* class QEngineT */

/**
 * \class qpp::QEngine
 * \brief Pure state quantum engine
 * \note Kept for backwards compatibility, use qpp::QKetEngine
 * \see qpp::QEngineT
 */
struct QEngine : public QEngineT<ket> {
    using QEngineT<ket>::QEngineT;
    // traits
    /**
     * \brief qpp::IQEngineTraits::traits_get_name() override
     */
    std::string traits_get_name() const override { return "QEngine"; }
    // end traits
};

/**
 * \class qpp::QKetEngine
 * \brief Pure state quantum engine
 * \see qpp::QEngineT
 */
struct QKetEngine : public QEngineT<ket> {
    using QEngineT<ket>::QEngineT;
    // traits
    /**
     * \brief qpp::IQEngineTraits::traits_get_name() override
     */
    std::string traits_get_name() const override { return "QKetEngine"; }
    // end traits
};

/**
 * \class qpp::QDensityEngine
 * \brief Mixed state quantum engine
 * \see qpp::QEngineT
 */
struct QDensityEngine : public QEngineT<cmat> {
    using QEngineT<cmat>::QEngineT;
    // traits
    /**
     * \brief qpp::IQEngineTraits::traits_get_name() override
     */
    std::string traits_get_name() const override { return "QDensityEngine"; }
    // end traits
};

} /* namespace qpp */

#endif /* QPP_CLASSES_QENGINE_HPP_ */
