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

#include <algorithm>
#include <cassert>
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
     * \brief Compute the position of the first
     * measurement/discard/reset/conditional step
     * \see qpp::internal::canonical_form()
     *
     * \param steps Vector of qpp::QCircuit::iterator
     * \return Position of the first measurement/discard/reset/conditional step
     */
    idx compute_optimize_up_to_pos_(
        const std::vector<QCircuit::iterator>& steps) const {
        // find the position of the first
        // measurement/reset/discard/conditional step
        auto optimize_up_to_pos_it =
            std::find_if(steps.begin(), steps.end(), [](auto&& elem) {
                return internal::is_measurement(elem) ||
                       internal::is_discard(elem) || internal::is_reset(elem) ||
                       internal::is_conditional(elem);
            });
        idx optimize_up_to_pos =
            std::distance(steps.begin(), optimize_up_to_pos_it);

        return optimize_up_to_pos;
    }

    /**
     * \brief Returns pair of (bool, idx), first true if the canonical form of
     * the circuit can be sampled from, second denoting the position of the
     * first measurement/reset/discard/conditional step
     * \see qpp::internal::canonical_form()
     *
     * \note A circuit can be sampled from if and only if it its canonical form
     * have only projective measurements (including post-selection) after the
     * first measurement/reset/discard/conditional step. In other words, it's of
     * the form [Gate step(s)] * [Projective measurements step(s)].
     *
     * \param steps Vector of qpp::QCircuit::iterator
     * \return Pair of (bool, idx), first true if the canonical form of the
     * circuit can be sampled from, second denoting the position of the first
     * measurement/reset/discard/conditional step
     */
    std::pair<bool, idx>
    can_sample_from_(const std::vector<QCircuit::iterator>& steps) const {

        // in the following, we will partition the circuit as
        // [0 ... optimize_up_to_pos ... end)

        // find the position of the first
        // measurement/reset/discard/conditional step
        idx optimize_up_to_pos = compute_optimize_up_to_pos_(steps);

        // decide if we can sample (every step after optimize_up_to_pos
        // must be a projective measurement)
        for (idx i = optimize_up_to_pos; i < steps.size();
             ++i) {
            if (!(internal::is_projective_measurement(steps[i])) ||
                internal::is_discard(steps[i])) {
                return {false, optimize_up_to_pos};
            }
        }

        return {true, optimize_up_to_pos};
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
                            return {measured, false};
                        }
                    }
                }
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
            }
        }

        return {measured, true};
    }

  protected:
    internal::QEngineState<T> qeng_st_; ///< current state of the engine
    internal::QEngineStatistics
        stats_{}; ///< measurement statistics for multiple runs

    // FIXME: doc
    /**
     * \brief Executes a contiguous series of projective measurement steps
     *
     * \param steps Vector of qpp::QCircuit::iterator
     * \param pos Index from where the execution starts
     * \param ensure_post_selection When \a ensure_post_selection is true, the
     * step is executed repeatedly until the post-selection succeeds, or until
     * the maximum number of post-selection reps is reached,
     * see qpp::QEngineT::set_max_post_selection_reps(), in which case the
     * post-selection is not guaranteed to succeed; check the state of the
     * engine, see qpp::QEngineT::post_select_ok()
     */
    void execute_circuit_steps_once_(std::vector<QCircuit::iterator>& steps,
                                     idx pos, bool ensure_post_selection) {
        // save the engine state
        auto engine_state_copy = qeng_st_;

        // lambda that restores qeng_st_.ensure_post_selection_ flag
        auto restore_ensure_post_selection_flag = [&] {
            if (ensure_post_selection) {
                qeng_st_.ensure_post_selection_ =
                    engine_state_copy.ensure_post_selection_;
            }
        };

        // sets the state of the engine to the entry state
        if (ensure_post_selection) {
            qeng_st_.ensure_post_selection_ = true;
        }

        bool measured = false;
        for (auto it = steps[pos]; it.get_ip() < steps.size(); ++it) {
            if (internal::is_measurement(it)) {
                measured = true;
            }
            this->execute(it);
            // post-selection failed, stop executing
            if (!qeng_st_.post_select_ok_) {
                restore_ensure_post_selection_flag();
                return;
            }
        }

        // at least one qudit was measured
        if (measured) {
            ++stats_.data()[get_dits()];
        }

        restore_ensure_post_selection_flag();
    }

    // FIXME: doc
    /**
     * \brief Executes a contiguous series of projective measurement steps
     *
     * \param steps Vector of qpp::QCircuit::iterator
     * \param pos Index from where the execution starts
     * \param reps Number of repetitions
     */
    void execute_no_sample_(std::vector<QCircuit::iterator>& steps, idx pos,
                            idx reps) {
        // save the engine state
        internal::QEngineState<T> engine_state_copy = qeng_st_;

        // this state will store the state of the engine for the first
        // successful post-selection
        internal::QEngineState<T> engine_state_ps_ok = qeng_st_;

        bool found_successful_rep = false;
        for (idx rep = 0; rep < reps; ++rep) {
            execute_circuit_steps_once_(steps, pos,
                                        qeng_st_.ensure_post_selection_);
            bool post_select_ok = this->post_select_ok();
            // save the first successful post-selection engine state
            if (post_select_ok && !found_successful_rep) {
                found_successful_rep = true;
                engine_state_ps_ok = qeng_st_;
            }
            // if not last rep, restore the engine state for the next run
            if (rep + 1 < reps) {
                qeng_st_ = engine_state_copy;
            }
        }

        // restore the engine state to the first successful post-selection
        // engine state
        if (found_successful_rep) {
            qeng_st_ = engine_state_ps_ok;
        }
    }

    /**
     * \brief Executes a contiguous series of projective measurement steps via
     * sampling
     *
     * \param steps Vector of qpp::QCircuit::iterator
     * \param pos Index from where the execution starts; from this index
     * further, all steps are assumed to be a projective measurement (including
     * post-selection)
     * \param reps Number of repetitions
     */
    void execute_sample_(std::vector<QCircuit::iterator>& steps, idx pos,
                         idx reps) {
        std::map<idx, std::pair<idx, std::optional<idx>>>
            c_q; // records the c <- (q, [OPTIONAL ps_val]) map
        std::map<idx, idx> q_ps_val; // records the q <- ps_val map
        std::set<idx> q_m;           // qudits that were measured projectively,
                                     // including post-selected

        auto [measured, post_select_compatible] =
            build_sampling_maps_(steps, pos, c_q, q_ps_val, q_m);

        if (!post_select_compatible) {
            qeng_st_.post_select_ok_ = false;
            return;
        }

        if (!measured) {
            return;
        }

        // at least one qudit was measured, add it to stats_

        // build the vector of measured qudits that we must sample
        // from
        std::set<idx> sample_from_set;
        for (auto [dit, qudit_ps_val] : c_q) {
            sample_from_set.emplace(qudit_ps_val.first);
        }

        std::vector<idx> sample_from(sample_from_set.begin(),
                                     sample_from_set.end());

        // save the engine state
        internal::QEngineState<T> engine_state_copy = qeng_st_;

        // this state will store the state of the engine for the first
        // successful post-selection
        internal::QEngineState<T> engine_state_ps_ok = qeng_st_;

        bool found_successful_rep = false;
        for (idx rep = 0; rep < reps; ++rep) {
            idx max_post_select_idx = 0;
        REPEAT_POST_SELECT_SAMPLE:
            ++max_post_select_idx;
            // sample from the quantum state
            std::vector<idx> sample_result_restricted_support = sample(
                engine_state_copy.qstate_, sample_from, this->qc_ptr_->get_d());

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
                qeng_st_.post_select_ok_ = false;
                if (qeng_st_.ensure_post_selection_) {
                    if (max_post_select_idx <
                        qeng_st_.max_post_selection_reps_) {
                        goto REPEAT_POST_SELECT_SAMPLE;
                    } else {
                        continue;
                    }
                } else {
                    continue;
                }
            } else {
                qeng_st_.post_select_ok_ = true;
            }

            // at this point we know for sure that the current rep succeeded
            assert(qeng_st_.post_select_ok_ == true);

            // save the first successful post-selection engine state and
            // repeatedly re-execute the rep until post-selection succeeds, so
            // we can compute the final quantum state
            if (!found_successful_rep) {
                found_successful_rep = true;
                qeng_st_.max_post_selection_reps_ =
                    std::numeric_limits<idx>::max();
                execute_circuit_steps_once_(steps, pos, true);
                qeng_st_.max_post_selection_reps_ =
                    engine_state_copy.max_post_selection_reps_;
                engine_state_ps_ok = qeng_st_;
            } else {
                // extend sample_result to full support
                std::vector<idx> sample_result = this->get_dits();
                idx i = 0;
                for (auto [dit, qudit_ps_val] : c_q) {
                    sample_result[dit] = sample_result_restricted_support[i++];
                }
                ++stats_.data()[sample_result];
            }
        }

        // set the engine state to the first successful post-selection
        // engine state
        if (found_successful_rep) {
            qeng_st_ = engine_state_ps_ok;
        }
    }

    /**
     * \brief Marks qudit \a i as measured then re-label accordingly the
     * remaining non-measured qudits
     *
     * \param i Qudit index
     */
    void set_measured_destructively_(idx i) {
        // EXCEPTION CHECKS
        if (was_measured_d(i)) {
            throw exception::QuditAlreadyMeasured(
                "qpp::QEngineT::set_measured_()", "i");
        }
        // END EXCEPTION CHECKS
        qeng_st_.subsys_[i] = std::numeric_limits<idx>::max(); // set qudit i to
                                                               // measured state
        for (idx m = i; m < this->qc_ptr_->get_nq(); ++m) {
            if (!was_measured_d(m)) {
                --qeng_st_.subsys_[m];
            }
        }
    }

    // giving a vector of non-measured qudits, get their relative position
    // w.r.t. the measured qudits
    /**
     * \brief Giving a vector \a v of non-measured qudits, gets their
     * relative position with respect to the measured qudits
     *
     * \param v Vector of non-measured qudit indexes
     * \return Vector of qudit indexes
     */
    std::vector<idx> get_relative_pos_(std::vector<idx> v) {
        idx vsize = v.size();
        for (idx i = 0; i < vsize; ++i) {
            // EXCEPTION CHECKS
            if (was_measured_d(v[i])) {
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
            // TODO: check if this can happen (st_.dits_.empty())
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

        idx max_post_select_idx = 0;
        switch (mstep_type_) {
            case MType::NONE:
                break;

            case MType::MEASURE:
            case MType::MEASURE_ND:
            case MType::POST_SELECT:
            case MType::POST_SELECT_ND:
                max_post_select_idx = 0;
            REPEAT_POST_SELECT:
                std::tie(results, probs, qstate) = measure_seq(
                    qeng_st_.qstate_, target_rel_pos, d, destructive);

                // POST_SELECT(_ND)
                if (mstep_type_ == MType::POST_SELECT ||
                    mstep_type_ == MType::POST_SELECT_ND) {
                    // post-selection failed
                    if (results[0] != measurement_step.ps_vals_.value()[0]) {
                        ++max_post_select_idx;
                        if (!qeng_st_.ensure_post_selection_) {
                            qeng_st_.post_select_ok_ = false;
                        } else if (max_post_select_idx <
                                   qeng_st_.max_post_selection_reps_) {
                            goto REPEAT_POST_SELECT;
                        }
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
                max_post_select_idx = 0;
            REPEAT_POST_SELECT_MANY:
                std::tie(results, probs, qstate) = measure_seq(
                    qeng_st_.qstate_, target_rel_pos, d, destructive);

                // POST_SELECT_MANY(_ND)
                if (mstep_type_ == MType::POST_SELECT_MANY ||
                    mstep_type_ == MType::POST_SELECT_MANY_ND) {
                    // post-selection failed
                    if (results != measurement_step.ps_vals_.value()) {
                        ++max_post_select_idx;
                        if (!qeng_st_.ensure_post_selection_) {
                            qeng_st_.post_select_ok_ = false;
                        } else if (max_post_select_idx <
                                   qeng_st_.max_post_selection_reps_) {
                            goto REPEAT_POST_SELECT_MANY;
                        } else {
                            qeng_st_.post_select_ok_ = false;
                            return;
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
                max_post_select_idx = 0;
            REPEAT_POST_SELECT_V:
                std::tie(mres, probs, states) = measure(
                    qeng_st_.qstate_, h_tbl[measurement_step.mats_hash_[0]],
                    target_rel_pos, d, destructive);

                // POST_SELECT_V(_ND)
                if (mstep_type_ == MType::POST_SELECT_V ||
                    mstep_type_ == MType::POST_SELECT_V_ND) {
                    // post-selection failed
                    if (mres != measurement_step.ps_vals_.value()[0]) {
                        ++max_post_select_idx;
                        if (!qeng_st_.ensure_post_selection_) {
                            qeng_st_.post_select_ok_ = false;
                        } else if (max_post_select_idx <
                                   qeng_st_.max_post_selection_reps_) {
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
                max_post_select_idx = 0;
            REPEAT_POST_SELECT_V_JOINT:
                std::tie(mres, probs, states) = measure(
                    qeng_st_.qstate_, h_tbl[measurement_step.mats_hash_[0]],
                    target_rel_pos, d, destructive);

                // POST_SELECT_V(_ND)
                if (mstep_type_ == MType::POST_SELECT_V_JOINT ||
                    mstep_type_ == MType::POST_SELECT_V_JOINT_ND) {
                    // post-selection failed
                    if (mres != measurement_step.ps_vals_.value()[0]) {
                        ++max_post_select_idx;
                        if (!qeng_st_.ensure_post_selection_) {
                            qeng_st_.post_select_ok_ = false;
                        } else if (max_post_select_idx <
                                   qeng_st_.max_post_selection_reps_) {
                            goto REPEAT_POST_SELECT_V_JOINT;
                        }
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
     * \brief Executes qpp::internal::QCircuitConditionalStep
     *
     * \param conditional_step Instance of
     * qpp::internal::QCircuitConditionalStep
     * \param it Iterator pointing to the step to be executed
     *
     */
    virtual void execute_conditional_step_(
        const internal::QCircuitConditionalStep& conditional_step,
        QCircuitTraits<QCircuit>::iterator_type it) {
        using Type = internal::QCircuitConditionalStep::Type;
        auto if_expr = conditional_step.ctx_.if_expr;
        auto else_expr = conditional_step.ctx_.else_expr;
        auto endif_expr = conditional_step.ctx_.endif_expr;
        auto dit_shift = conditional_step.ctx_.dit_shift;
        bool is_true = true;
        switch (conditional_step.condition_type_) {
            case Type::IF:
                if (if_expr.has_value()) {
                    auto lambda = if_expr.value().second;
                    auto it_begin =
                        std::next(qeng_st_.dits_.begin(), dit_shift);
                    auto it_end = qeng_st_.dits_.end();
                    is_true = lambda(std::vector<idx>(it_begin, it_end));
                }
                LOG << "Executing IF statement -> " << std::boolalpha << is_true
                    << "\n";
                // jump on false
                if (!is_true) {
                    idx adv = else_expr.has_value() ? else_expr.value()
                                                    : endif_expr.value();
                    adv -= if_expr.value().first;

                    it.advance(adv - 1);
                    LOG << "\t\texecute_conditional_step_(): " << it
                        << std::endl;
                } else if (else_expr.has_value()) {
                    idx adv = endif_expr.value();
                    adv -= else_expr.value();
                    it.advance(adv - 1);
                }
                break;
            case Type::ELSE: {
                LOG << "Executing ELSE statement\n";
                break;
            }
            case Type::ENDIF:
                LOG << "Executing ENDIF statement\n";
                break;
            case Type::NONE:
                break;
        }

        return;
    }

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
     * \brief Constructs a quantum engine out of a quantum circuit
     * description
     *
     * \note The quantum circuit description must be an lvalue
     * \see qpp::QEngine(QCircuit&&)
     *
     * \note The initial underlying quantum state is set to
     * \f$|0\rangle^{\otimes n}\f$
     *
     * \param qc Quantum circuit description
     */
    explicit QEngineT(const QCircuit& qc)
        : QBaseEngine<T, QCircuit>{qc}, qeng_st_{this->qc_ptr_}, stats_{} {}

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
     * \see qpp::QEngineT::set_dits()
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
     * \see qpp::QEngineT::set_dit()
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
    bool was_measured_d(idx i) const {
        return qeng_st_.subsys_[i] == std::numeric_limits<idx>::max();
    }

    /**
     * \brief Vector of already measured (destructively) qudit indexes at
     * the current engine state
     *
     * \return Vector of already measured qudit (destructively) indexes at
     * the current engine state
     */
    std::vector<idx> get_measured_d() const {
        std::vector<idx> result;
        for (idx i = 0; i < this->qc_ptr_->get_nq(); ++i) {
            if (was_measured_d(i)) {
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
        for (idx i = 0; i < this->qc_ptr_->get_nq(); ++i) {
            if (!was_measured_d(i)) {
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
     * \brief Returns true if post-selection was successful (or absent),
     * false otherwise
     */
    bool post_select_ok() const { return qeng_st_.post_select_ok_; }

    /**
     * \brief Returns true if post-selection is enforced (must succeed), false
     * otherwise
     */
    bool get_ensure_post_selection() const {
        return qeng_st_.ensure_post_selection_;
    }

    /**
     * \brief Maximum number of executions of a circuit post-selection step
     * until success
     *
     * \return Maximum number of executions of a post-selection circuit step
     * until success
     */
    idx get_max_post_selection_reps() const {
        return qeng_st_.max_post_selection_reps_;
    }
    // end getters

    // setters
    /**
     * \brief Sets the classical dit at position \a i
     * \see qpp::QEngineT::get_dit()
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
     * \see qpp::QEngineT::get_dits()
     *
     * \note When interfacing with OpenQASM, the classical dits/registers
     * are evaluated in little-endian order, with  the least significant bit
     * being stored first. For example, [1,0,0] is interpreted as 1 (and not
     * 4).
     *
     * \param dits Vector of classical dits, must have the same size as the
     * internal vector of classical dits returned by
     * qpp::QEngineT::get_dits()
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
     * \brief Enforces post-selection (must succeed) when \a val is true
     *
     * \param val If true, repeatedly executes post-selection steps until the
     * post-selection result(s) agree, or until the maximum number of
     * post-selection repetitions is reached, see
     * qpp::QEngineT::set_max_post_selection_reps(), in which case the
     * post-selection is not guaranteed to succeed; check the state of the
     * engine, see qpp::QEngineT::post_select_ok().
     * \return Reference to the current instance
     */
    QEngineT& set_ensure_post_selection(bool val) {
        qeng_st_.ensure_post_selection_ = val;
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
        idx n = get_non_measured_d().size();
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

    /**
     * \brief Sets the maximum number of executions of a circuit
     * post-selection step until success
     *
     * \param max_post_selection_reps Maximum number of executions of a
     * post-selection step until success
     * \return Reference to the current instance
     */
    QEngineT& set_max_post_selection_reps(idx max_post_selection_reps) {
        qeng_st_.max_post_selection_reps_ = max_post_selection_reps;
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
     */
    virtual QEngineT& reset(bool reset_stats = true) {
        qeng_st_.reset(qeng_st_.qstate_);

        if (reset_stats) {
            this->reset_stats();
        }

        return *this;
    }

    /**
     * \brief Executes one step in the quantum circuit description
     *
     * \note Override only this member function in every derived class to
     * achieve the desired behaviour
     *
     * \param it Iterator pointing to the step to be executed
     * \return Reference to the current instance
     */
    QEngineT& execute(QCircuitTraits<QCircuit>::iterator_type it) override {
        // EXCEPTION CHECKS
        // iterator must point to the same quantum circuit description
        if (it->get_qc_ptr() != this->qc_ptr_) {
            throw exception::InvalidIterator(
                "qpp::QEngineT::execute()",
                "Iterator does not point to the same circuit description");
        }
        if (!this->qc_ptr_->validate_conditionals()) {
            throw exception::InvalidConditional("qpp::QEngineT::execute())",
                                                "Missing ENDIF");
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

        auto conditional_step_visitor =
            [&](const internal::QCircuitConditionalStep& conditional_step) {
                return this->execute_conditional_step_(conditional_step, it);
            };

        std::visit(
            overloaded{
                conditional_step_visitor,
                gate_step_visitor,
                measurement_step_visitor,
                nop_step_visitor,
            },
            it->get_step());

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
            throw exception::OutOfRange("qpp::QEngineT::execute()", "reps");
        }
        if (!this->qc_ptr_->validate_conditionals()) {
            throw exception::InvalidConditional("qpp::QEngineT::execute()",
                                                "Missing ENDIF");
        }
        // END EXCEPTION CHECKS

        // save the engine state
        auto engine_state_copy = qeng_st_;

        this->reset(false); // keep the statistics
        this->set_ensure_post_selection(
            engine_state_copy.ensure_post_selection_);
        this->set_max_post_selection_reps(
            engine_state_copy.max_post_selection_reps_);

        // TODO: change back first to canonical_form
        auto steps_as_iterators =
            reps > 1 ? internal::circuit_as_iterators(*this->qc_ptr_)
                     : internal::circuit_as_iterators(*this->qc_ptr_);
        if (steps_as_iterators.empty()) {
            return *this;
        }

        auto [can_sample, optimize_up_to_pos] =
            can_sample_from_(steps_as_iterators);
        qeng_st_.can_sample_ = reps > 1 && can_sample;

        // execute everything ONCE in the interval
        // [0, first_measurement_discard_reset_pos)
        for (auto it = steps_as_iterators[0];
             it.get_ip() < optimize_up_to_pos; ++it) {
            execute(it);
        }

        // TODO: comment the line below in production
        qeng_st_.can_sample_ = false;

        // execute repeatedly everything in the remaining interval
        // can sample: every step from now on is a projective measurement
        if (qeng_st_.can_sample_) {
            execute_sample_(steps_as_iterators,
                            optimize_up_to_pos, reps);
        }
        // cannot sample
        else {
            execute_no_sample_(steps_as_iterators,
                               optimize_up_to_pos, reps);
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
        ss << disp(get_measured_d(), IOManipContainerOpts{}.set_sep(", "));
        result += ss.str() + ", ";

        ss.str("");
        ss.clear();
        ss << disp(get_non_measured_d(), IOManipContainerOpts{}.set_sep(", "));
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
           << disp(get_dits(), IOManipContainerOpts{}.set_sep(", "));

        // compute the statistics
        if (!stats_.data().empty()) {
            os << '\n' << stats_;
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
