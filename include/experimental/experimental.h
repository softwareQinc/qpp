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
struct QCircuit : public IDisplay {
    idx nq_;                     ///< number of qudits
    idx nc_;                     ///< number of classical "dits"
    idx d_;                      ///< dimension
    ket psi_;                    ///< state vector
    std::vector<bool> measured_; ///< keeps track of measured qudits
    std::vector<idx> dits_;      ///< classical dits

    /**
     * \brief Type of operation being executed at one step
     */
    enum class GateType {
        SINGLE, ///< unitary gate on a single qudit

        TWO, ///< unitary gate on 2 qudits

        THREE, ///< unitary gate on 3 qudits

        FAN, ///< same unitary gate on multiple qudits

        CUSTOM, ///< custom gate on multiple qudits

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

        NONE, ///< signals no gate
    };
    enum class MeasureType {
        MEASURE_Z,  ///< Z measurement of single qubit/qudit
        MEASURE_V,  ///< measurement of single qubit/qudit in the orthonormal
                    ///< basis or rank-1 POVM specified by the columns of matrix
                    ///< \a V
        MEASURE_KS, ///< generalized measurement of single qubit/qudit with
                    ///< Kraus operators Ks
        NONE,       ///< signals no measurement
    };

    /**
     * \brief One step consisting of unitary gates in the circuit
     */
    struct GateStep {
        GateType gate_type_ = GateType::NONE; ///< gate type
        cmat gate_;                           ///< gate
        std::vector<idx> ctrl_;               ///< control
        std::vector<idx> target_; ///< target where the gate is being applied
        std::string name_;        ///< custom name of the step
        GateStep() = default;
        GateStep(GateType gate_type, const cmat& gate,
                 const std::vector<idx>& ctrl, const std::vector<idx>& target,
                 const std::string& name = "")
            : gate_type_{gate_type}, gate_{gate}, ctrl_{ctrl}, target_{target},
              name_{name} {}
    };

    struct MeasureStep {
        MeasureType measurement_type_ = MeasureType::NONE; ///< measurement type
        std::vector<cmat> mats_;  ///< matrix/matrices being measured
        std::vector<idx> target_; ///< target where the measurement is applied
        std::string name_;        ///< custom name of the step
        MeasureStep() = default;
        MeasureStep(MeasureType measurement_type, const std::vector<cmat>& mats,
                    const std::vector<idx>& target,
                    const std::string& name = "")
            : measurement_type_{measurement_type}, mats_{mats}, target_{target},
              name_{name} {}
    };

    friend std::ostream& operator<<(std::ostream& os,
                                    const GateType& gate_type) {
        switch (gate_type) {
        case GateType::SINGLE:
            return os << "SINGLE";
        case GateType::TWO:
            return os << "TWO";
        case GateType::THREE:
            return os << "THREE";
        case GateType::FAN:
            return os << "FAN";
        case GateType::CUSTOM:
            return os << "CUSTOM";
        case GateType::SINGLE_CTRL_SINGLE_TARGET:
            return os << "SINGLE_CTRL_SINGLE_TARGET";
        case GateType::SINGLE_CTRL_MULTIPLE_TARGET:
            return os << "SINGLE_CTRL_MULTIPLE_TARGET";
        case GateType::MULTIPLE_CTRL_SINGLE_TARGET:
            return os << "MULTIPLE_CTRL_SINGLE_TARGET";
        case GateType::MULTIPLE_CTRL_MULTIPLE_TARGET:
            return os << "MULTIPLE_CTRL_MULTIPLE_TARGET";
        case GateType::CUSTOM_CTRL:
            return os << "CUSTOM_CTRL";
        case GateType::SINGLE_cCTRL_SINGLE_TARGET:
            return os << "SINGLE_cCTRL_SINGLE_TARGET";
        case GateType::SINGLE_cCTRL_MULTIPLE_TARGET:
            return os << "SINGLE_cCTRL_MULTIPLE_TARGET";
        case GateType::MULTIPLE_cCTRL_SINGLE_TARGET:
            return os << "MULTIPLE_cCTRL_SINGLE_TARGET";
        case GateType::MULTIPLE_cCTRL_MULTIPLE_TARGET:
            return os << "MULTIPLE_cCTRL_MULTIPLE_TARGET";
        case GateType::CUSTOM_cCTRL:
            return os << "CUSTOM_cCTRL";
        case GateType::NONE:
            return os << "NONE";
        }
    }

    friend std::ostream& operator<<(std::ostream& os,
                                    const MeasureType& measure_type) {
        switch (measure_type) {
        case MeasureType::MEASURE_Z:
            return os << "MEASURE_Z";
        case MeasureType::MEASURE_V:
            return os << "MEASURE_V";
        case MeasureType::MEASURE_KS:
            return os << "MEASURE_KS";
        case MeasureType::NONE:
            return os << "NONE";
        }
    }

    ///< quantum circuit representation
    std::vector<GateStep> gates_;
    std::vector<MeasureStep> measurements_;

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
    QCircuit(idx nq, idx nc = 0, idx d = 2)
        : nq_{nq}, nc_{nc}, d_{d}, psi_{st.zero(nq_, d_)},
          measured_(nq_, false), dits_(nc_, 0) {}

    // single gate single qubit/qudit
    void apply(const cmat& U, idx i, const std::string& name = "") {
        gates_.emplace_back(GateType::SINGLE, U, std::vector<idx>{},
                            std::vector<idx>{i}, name);
    }
    // single gate 2 qubits/qudits
    void apply(const cmat& U, idx i, idx j, const std::string& name = "") {
        gates_.emplace_back(GateType::TWO, U, std::vector<idx>{},
                            std::vector<idx>{i, j}, name);
    }
    // single gate 3 qubits/qudits
    void apply(const cmat& U, idx i, idx j, idx k,
               const std::string& name = "") {
        gates_.emplace_back(GateType::THREE, U, std::vector<idx>{},
                            std::vector<idx>{i, j, k}, name);
    }

    // multiple qubits/qudits same gate
    void apply(const cmat& U, const std::vector<idx>& target,
               const std::string& name = "") {
        gates_.emplace_back(GateType::FAN, U, std::vector<idx>{}, target, name);
    }

    // custom gate
    void apply_custom(const cmat& U, const std::vector<idx>& target,
                      const std::string& name = "") {
        gates_.emplace_back(GateType::CUSTOM, U, std::vector<idx>{}, target,
                            name);
    }

    // single ctrl single target
    void CTRL(const cmat& U, idx ctrl, idx target, const std::string& name) {
        gates_.emplace_back(GateType::SINGLE_CTRL_SINGLE_TARGET, U,
                            std::vector<idx>{ctrl}, std::vector<idx>{target},
                            name);
    }

    // single ctrl multiple target
    void CTRL(const cmat& U, idx ctrl, const std::vector<idx>& target,
              const std::string& name) {
        gates_.emplace_back(GateType::SINGLE_CTRL_MULTIPLE_TARGET, U,
                            std::vector<idx>{ctrl}, target, name);
    }

    // multiple ctrl single target
    void CTRL(const cmat& U, const std::vector<idx>& ctrl, idx target,
              const std::string& name) {
        gates_.emplace_back(GateType::MULTIPLE_CTRL_SINGLE_TARGET, U, ctrl,
                            std::vector<idx>{target}, name);
    }

    // multiple ctrl multiple target
    void CTRL(const cmat& U, const std::vector<idx>& ctrl,
              const std::vector<idx>& target, const std::string& name) {
        gates_.emplace_back(GateType::MULTIPLE_CTRL_MULTIPLE_TARGET, U, ctrl,
                            std::vector<idx>{target}, name);
    }

    //  custom controlled gate with multiple controls and multiple targets
    void CTRL_custom(const cmat& U, const std::vector<idx>& ctrl,
                     const std::vector<idx>& target, const std::string& name) {
        gates_.emplace_back(GateType::CUSTOM_CTRL, U, ctrl, target, name);
    }

    // single ctrl single target
    void cCTRL(const cmat& U, idx ctrl, idx target, const std::string& name) {
        gates_.emplace_back(GateType::SINGLE_cCTRL_SINGLE_TARGET, U,
                            std::vector<idx>{ctrl}, std::vector<idx>{target},
                            name);
    }

    // single ctrl multiple target
    void cCTRL(const cmat& U, idx ctrl, const std::vector<idx>& target,
               const std::string& name) {
        gates_.emplace_back(GateType::SINGLE_cCTRL_MULTIPLE_TARGET, U,
                            std::vector<idx>{ctrl}, target, name);
    }

    // multiple ctrl single target
    void cCTRL(const cmat& U, const std::vector<idx>& ctrl, idx target,
               const std::string& name) {
        gates_.emplace_back(GateType::MULTIPLE_cCTRL_SINGLE_TARGET, U, ctrl,
                            std::vector<idx>{target}, name);
    }

    // multiple ctrl multiple target
    void cCTRL(const cmat& U, const std::vector<idx>& ctrl,
               const std::vector<idx>& target, const std::string& name) {
        gates_.emplace_back(GateType::MULTIPLE_cCTRL_MULTIPLE_TARGET, U, ctrl,
                            std::vector<idx>{target}, name);
    }

    //  custom controlled gate with multiple controls and multiple targets
    void cCTRL_custom(const cmat& U, const std::vector<idx>& ctrl,
                      const std::vector<idx>& target, const std::string& name) {
        gates_.emplace_back(GateType::CUSTOM_cCTRL, U, ctrl, target, name);
    }

    // Z measurement on single qudit
    void measureZ(idx i, const std::string& name) {
        measurements_.emplace_back(MeasureType::MEASURE_Z, std::vector<cmat>{},
                                   std::vector<idx>{i}, name);
    }

    std::ostream& display(std::ostream& os) const override {
        for (auto&& elem : gates_) {
            os << elem.gate_type_ << " "
               << "'" << elem.name_ << "'"
               << " ";
            os << disp(elem.ctrl_, ",") << " " << disp(elem.target_, ",");
            os << '\n';
        }
        for (auto&& elem : measurements_) {
            os << elem.measurement_type_ << " "
               << "'" << elem.name_ << "'"
               << " ";
            os << disp(elem.target_, ",");
            os << '\n';
        }

        return os;
    }
}; // namespace experimental

} // namespace experimental
} /* namespace qpp */

#endif /* EXPERIMENTAL_EXPERIMENTAL_H_ */
