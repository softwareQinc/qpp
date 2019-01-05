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
template <class T>
class QCirc {
    idx d_;                      // dimension, qubit circuit (d = 2) by default
    idx nq_, nc_;                // quantum bits, classical bits/dits
    ket psi_;                    // current state of the circuit
    std::vector<bool> measured_; // set elements to one if the corresponding
                                 // qubit had been measured
    std::vector<idx> results_;   // get_results of the measurements
    std::vector<idx> bits_;      // store the values of the classical bits/dits

  protected:
    std::vector<idx> update_subsys_(const std::vector<idx>& subsys) {
        std::vector<idx> result = subsys;
        idx subsys_size = subsys.size();
        for (idx i = 0; i < subsys_size; ++i) {
            for (idx m = 0; m < subsys[i]; ++m) {
                if (measured_[m]) { // if the qubit m was measured
                    --result[i];
                }
            }
        }
        std::sort(std::begin(result), std::end(result), std::less<idx>{});
        return result;
    }

  public:
    QCirc(idx nq, idx nc, idx d = 2)
        : d_{d}, nq_{nq}, nc_{nc}, psi_{st.zero(nq_, d_)},
          measured_(nq_, false), results_(nq_, -1), bits_(nc_, 0) {}

    // destructive measurement
    void measure(std::vector<idx> subsys) {
        // sort subsys in decreasing order
        idx subsys_size = subsys.size();
        for (idx i = 0; i < subsys_size; ++i) {
            if (measured_[subsys[i]])
                throw exception::CustomException(
                    "qpp::QCirc::measure()", "Subsystem was measured before");
        }
        // update subsystem labels
        std::vector<idx> subsys_updated = update_subsys_(subsys);
        for (idx i = 0; i < subsys_size; ++i) {
            measured_[subsys[i]] = true;
        }
        auto m = measure_seq(psi_, subsys_updated, d_);
        auto result = std::get<0>(m); // measurement result
        for (idx i = 0; i < subsys_size; ++i)
            results_[subsys[i]] = result[i];
        // update psi_
        psi_ = std::get<2>(m);
        // std::cout << disp(psi_) << '\n';
    }

    // destructive measurement of all remaining qubits
    void measure_all() {
        std::vector<idx> subsys;
        for (idx i = 0; i < nq_; ++i) {
            if (!measured_[i])
                subsys.push_back(i);
        }
        // update subsystem labels
        std::vector<idx> subsys_updated = update_subsys_(subsys);
        if (subsys.size() != 0)
            this->measure(subsys);
        else
            throw exception::CustomException("qpp::QCirc::measure_all()",
                                             "All qubits were measured before");
    }

    void apply(const cmat& gate, const std::vector<idx>& subsys) {
        psi_ = qpp::apply(psi_, gate, update_subsys_(subsys), d_);
    }

    void applyCTRL(const cmat& gate, const std::vector<idx>& ctrl,
                   const std::vector<idx>& target) {
        psi_ = qpp::applyCTRL(psi_, gate, ctrl, target, d_);
    }

    void apply_all(const cmat& gate) {
        for (idx i = 0; i < nq_; ++i)
            if (!measured_[i]) {
                psi_ = qpp::apply(psi_, gate, update_subsys_({i}), d_);
            }
    }

    // resets the circuit
    void reset() {
        psi_ = st.zero(nq_);
        measured_ = std::vector<bool>(nq_, false);
        results_ = std::vector<idx>(nq_, -1);
        bits_ = std::vector<idx>(nc_, 0);
    }

    /* getters */
    // local dimension
    idx d() const noexcept { return d_; }

    // total number of qubits, regardless of being measured or not
    idx get_nq() const noexcept { return nq_; }

    // total number of classical bits/dits
    idx get_nc() const noexcept { return nc_; }

    // total number of qubits AND classical bits/dits
    idx get_size() const noexcept { return nq_ + nc_; }

    // total number of measured qubits
    idx get_num_measured_qubits() const noexcept {
        return std::count(std::begin(measured_), std::end(measured_), true);
    }

    // total number of non-measured qubits
    idx get_num_active_qubits() const noexcept {
        return this->get_nq() - this->get_num_measured_qubits();
    }

    // returns the up-to-date quantum state
    ket get_psi() const { return psi_; }

    // returns the results of the measurements
    std::vector<idx> get_results() const { return results_; }

    // returns the current result of the measurement as an integer,
    // ignores non-measured qubits
    idx get_results_as_N() const {
        std::vector<idx> tmp;
        for (idx i = 0; i < nq_; ++i) {
            if (measured_[i])
                tmp.push_back(results_[i]);
        }

        return multiidx2n(tmp, std::vector<idx>(tmp.size(), d_));
    }

    // returns the classical bits/dits array
    std::vector<idx>& bits() noexcept { return bits_; }
    /* end getters */
};

class LogicalCircuit : public IDisplay {
    using idx_vec = std::vector<idx>;
    // gate, gate name, control, targer
    using elem_type = std::tuple<cmat, std::string, idx_vec, idx_vec>;
    std::vector<elem_type> gates_;
    idx gate_count = 0;

  public:
    void add(const cmat& gate, const std::string& gate_name,
             const idx_vec& ctrl, const idx_vec& target) {
        auto elem = std::make_tuple(gate, gate_name, ctrl, target);
        gates_.push_back(elem);
        ++gate_count;
    }

    idx get_gate_count() const { return gate_count; }

    std::ostream& display(std::ostream& os) const override {
        os << "[";
        bool first = true;
        for (auto&& elem : gates_) {
            if (first) {
                os << "(";
                first = false;
            } else {
                os << ", (";
            }
            os << std::get<1>(elem) << ", ";
            os << disp(std::get<2>(elem), ", ");
            os << ", ";
            os << disp(std::get<3>(elem), ", ");
            os << ")";
        }
        os << "]";

        return os;
    }
};

class Test {
    idx nq_;                     ///< number of qubits/qudits
    idx nc_;                     ///< number of classical "dits"
    idx d_;                      ///< dimension
    ket psi_;                    ///< state vector
    std::vector<bool> measured_; ///< keeps track of measured qubits/qudits
    std::vector<idx> results_;   ///< keeps track of the measurement results
    std::vector<idx> dits_;      ///< classical bits/dits

    /**
     * \brief Type of operation being executed at one step
     */
    enum class Operation {
        GATE,       ///< unitary gate
        CTRL_GATE,  ///< controlled unitary gate
        MEASURE_Z,  ///< Z measurement
        MEASURE_V,  ///< measurement in the orthonormal basis or rank-1 POVM
                    ///< specified by the columns of matrix \a V
        MEASURE_KS, ///< generalized measurement with Kraus operators Ks
    };

    /**
     * \brief One step in the circuit
     */
    struct StepType {
        Operation op_;            ///< operation
        std::vector<cmat> mats_;  ///< matrix/matrices being applied/measured
        std::vector<idx> ctrl_;   ///< control (empty for measurements)
        std::vector<idx> target_; ///< target where the gate/measurement is
                                  ///< being applied
        std::string name_;        ///< custom name of the step
        StepType(Operation op, const std::vector<cmat>& mats,
                 const std::vector<idx>& ctrl, const std::vector<idx>& target,
                 const std::string& name)
            : op_{op}, mats_{mats}, ctrl_{ctrl}, target_{target}, name_{name} {}
    };

    std::vector<StepType> circuit_; ///< quantum circuit representation
  protected:
    std::vector<idx> update_subsys_(const std::vector<idx>& subsys) {
        std::vector<idx> result = subsys;
        idx subsys_size = subsys.size();
        for (idx i = 0; i < subsys_size; ++i) {
            for (idx m = 0; m < subsys[i]; ++m) {
                if (measured_[m]) { // if the qubit m was measured
                    --result[i];
                }
            }
        }
        std::sort(std::begin(result), std::end(result), std::less<idx>{});
        return result;
    }

  public:
    Test(idx nq, idx nc = 0, idx d = 2)
        : nq_{nq}, nc_{nc}, d_{d}, psi_{st.zero(nq_, d_)},
          measured_(nq_, false), results_(nq_, -1), dits_(nc_, 0) {}

    /**
     * \brief Apply quantum gate
     *
     * \param A Eigen expression (quantum gate)
     * \param target Subsystem indexes where the gate \a A is applied
     * \param name Optional name of the operation
     */
    template <typename Derived>
    void apply(const Eigen::MatrixBase<Derived>& A,
               const std::vector<idx>& target, const std::string& name = "") {
        circuit_.emplace_back(Operation::GATE, std::vector<cmat>{A},
                              std::vector<idx>{}, target, name);
    }
    /**
     * \brief Apply controlled quantum gate

     * \param A Eigen expression (quantum gate)
     * \param ctrl Control subsystem indexes
     * \param target Subsystem indexes where the gate \a A is applied
     * \param name Optional name of the operation
     */
    template <typename Derived>
    void applyCTRL(const Eigen::MatrixBase<Derived>& A,
                   const std::vector<idx>& ctrl, const std::vector<idx>& target,
                   const std::string& name = "") {
        circuit_.emplace_back(Operation::CTRL_GATE, std::vector<cmat>{A}, ctrl,
                              target, name);
    }

    /**
     * \brief Measures the part \a target of state vector in the Z basis
     *
     * \param target Subsystem indexes that are measured
     * \param name Optional name of the operation
     */
    void measureZ(const std::vector<idx>& target,
                  const std::string& name = "") {
        circuit_.emplace_back(Operation::MEASURE_Z, std::vector<cmat>{},
                              std::vector<idx>{}, target, name);
    }

    /**
     * \brief Measures the part \a target of state vector in the orthonormal
     * basis or rank-1 POVM specified by the matrix \a V
     *
     * \param V Matrix whose columns represent the measurement basis vectors or
     *        the bra parts of the rank-1 POVM
     * \param target Subsystem indexes that are measured
     * \param name Optional name of the operation
     */
    void measureV(const cmat& V, const std::vector<idx>& target,
                  const std::string& name = "") {
        circuit_.emplace_back(Operation::MEASURE_V, std::vector<cmat>{V},
                              std::vector<idx>{}, target, name);
    }

    /**
     * \brief Measures the part \a target of state vector using the set of
     * Kraus operators \a Ks
     *
     * \param Ks Set of Kraus operators
     * \param target Subsystem indexes that are measured
     */
    void measureKs(const std::vector<cmat>& Ks, const std::vector<idx>& target,
                   const std::string& name = "") {
        circuit_.emplace_back(Operation::MEASURE_KS, Ks, std::vector<idx>{},
                              target, name);
    }
};

} /* namespace experimental */
} /* namespace qpp */

#endif /* EXPERIMENTAL_EXPERIMENTAL_H_ */
