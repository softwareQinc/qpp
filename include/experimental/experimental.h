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

// TODO: display in JSon format

/**
 * \class qpp::QCircuitDescription
 * \brief Quantum circuit description class
 * \see qpp::QCircuit
 */
class QCircuitDescription : public IDisplay {
    const idx nq_;                         ///< number of qudits
    const idx nc_;                         ///< number of classical "dits"
    const idx d_;                          ///< dimension
    std::vector<idx> measurement_steps_{}; ///< keeps track of where the
    ///< measurements take place
    std::string name_;           ///< optional circuit name
    std::vector<bool> measured_; ///< keeps track of the measured qudits
    idx step_cnt_;               ///< step counter

  public:
    /**
     * \brief Type of gate being executed at one step
     */
    enum class GateType {
        NONE, ///< represents no gate

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
    /**
     * \brief Type of measurement being executed at one step
     */
    enum class MeasureType {
        NONE, ///< represents no measurement

        MEASURE_Z, ///< Z measurement of single qudit

        MEASURE_V, ///< measurement of single qudit in the orthonormal basis
        ///< or rank-1 projectors specified by the columns of matrix \a V

        MEASURE_V_MANY, ///< measurement of multiple qudits in the orthonormal
        ///< basis or rank-1 projectors specified by the columns of matrix
        ///< \a V
    };

  private:
    /**
     * \brief One step consisting only of gates/operators in the circuit
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
                 idx step_no, std::string name = "")
            : gate_type_{gate_type}, gate_{gate}, ctrl_{ctrl}, target_{target},
              step_no_{step_no}, name_{name} {}
    };

    /**
     * \brief One step consisting only of measurements in the circuit
     */
    struct MeasureStep {
        MeasureType measurement_type_ = MeasureType::NONE; ///< measurement type
        std::vector<cmat> mats_;  ///< matrix/matrices that specify the
                                  /// measurement
        std::vector<idx> target_; ///< target where the measurement is applied
        idx c_reg_{}; ///< index of the classical register where the measurement
        ///< result is being stored
        idx step_no_;      ///< step number
        std::string name_; ///< custom name of the step
        MeasureStep() = default;
        MeasureStep(MeasureType measurement_type, const std::vector<cmat>& mats,
                    const std::vector<idx>& target, idx c_reg, idx step_no,
                    std::string name = "")
            : measurement_type_{measurement_type}, mats_{mats}, target_{target},
              c_reg_{c_reg}, step_no_{step_no}, name_{name} {}
    };

    /**
     * \brief Extraction operator overload for
     * qpp::QCircuitDescription::GateType enum class
     *
     * \param os Output stream
     * \param gate_type qpp::QCircuitDescription::GateType enum class
     * \return Output stream
     */
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

    /**
     * \brief Extraction operator overload for
     * qpp::QCircuitDescription::MeasureType enum class
     *
     * \param os Output stream
     * \param gate_type qpp::QCircuitDescription::MeasureType enum class
     * \return Output stream
     */
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

    /**
     * \brief Extraction operator overload for
     * qpp::QCircuitDescription::GateStep class
     *
     * \param os Output stream
     * \param gate_type qpp::QCircuitDescription::GateStep class
     * \return Output stream
     */
    friend std::ostream& operator<<(std::ostream& os,
                                    const GateStep& gate_step) {
        os << gate_step.gate_type_ << ", ";
        if (gate_step.gate_type_ >= GateType::SINGLE_CTRL_SINGLE_TARGET)
            os << "ctrl = " << disp(gate_step.ctrl_, ",") << ", ";
        os << "target = " << disp(gate_step.target_, ",") << ", ";
        os << "name = " << '\"' << gate_step.name_ << '\"';

        return os;
    }

    /**
     * \brief Extraction operator overload for
     * qpp::QCircuitDescription::MeasureStep class
     *
     * \param os Output stream
     * \param gate_type qpp::QCircuitDescription::MeasureStep enum class
     * \return Output stream
     */
    friend std::ostream& operator<<(std::ostream& os,
                                    const MeasureStep& measure_step) {
        os << measure_step.measurement_type_ << ", ";
        os << "target = " << disp(measure_step.target_, ",") << ", ";
        os << "c_reg = " << measure_step.c_reg_ << ", ";
        os << "name = " << '\"' << measure_step.name_ << '\"';
        os << " ";

        return os;
    }

    std::vector<GateStep> gates_{};           ///< gates
    std::vector<MeasureStep> measurements_{}; ///< measurements

  public:
    /**
     * \brief Constructs a quantum circuit description
     *
     * \note The measurement results can only be stored in the classical dits
     * of which number is specified by \a nc
     *
     * \param nq Number of qbits
     * \param nc Number of classical dits
     * \param d Subsystem dimensions (optional, default is qubit, i.e. \a d = 2)
     * \param name Circuit description name (optional)
     */
    QCircuitDescription(idx nq, idx nc = 0, idx d = 2, std::string name = "")
        : nq_{nq}, nc_{nc}, d_{d}, name_{name}, measured_(nq, false),
          step_cnt_{0} {
        // EXCEPTION CHECKS

        if (nq == 0)
            throw exception::ZeroSize(
                "qpp::QCircuitDescription::QCircuitDescription()");
        if (d < 2)
            throw exception::OutOfRange(
                "qpp::QCircuitDescription::QCircuitDescription()");
        // END EXCEPTION CHECKS
    }

    // getters
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
     * \brief Local dimension of the comprising qudits
     *
     * \return Local dimension
     */
    idx get_d() const noexcept { return d_; }

    /**
     * \brief Vector of measurement positions in the circuit, i.e. the indexes
     * where the measurements take place
     *
     * \note If there are more consecutive measurements after step S, then their
     * indexes will all be S, i.e. it is always assumed that the measurements
     * taking place immediately after a gate step have the same index as the
     * preceding gate step.
     *
     * \return Vector of measurement positions
     */
    std::vector<idx> get_measurement_steps() const {
        return measurement_steps_;
    }

    /**
     * \brief Vector of qpp::QCircuitDescription::MeasureStep
     *
     * \return Vector of qpp::QCircuitDescription::MeasureStep
     */
    const std::vector<MeasureStep>& get_measurements() const noexcept {
        return measurements_;
    }

    /**
     * \brief Vector of qpp::QCircuitDescription::GateStep
     *
     * \return Vector of qpp::QCircuitDescription::GateStep
     */
    const std::vector<GateStep>& get_gates() const noexcept { return gates_; }

    /**
     * \brief Quantum circuit name
     *
     * \return Quantum circuit name
     */
    std::string get_name() const { return name_; }

    /**
     * \brief Check whether qudit \a i was already measured
     * \param i Qudit index
     * \return True if qudit \a i was already measured, false othwewise
     */
    idx get_measured(idx i) const {
        // EXCEPTION CHECKS

        if (i > nq_)
            throw exception::OutOfRange(
                "qpp::QCircuitDescription::get_measured()");
        // END EXCEPTION CHECKS

        return measured_[i];
    }

    /**
     * \brief Vector of already measured qudit indexes
     *
     * \return Vector of already measured qudit indexes
     */
    std::vector<idx> get_measured() const {
        std::vector<idx> result;
        for (idx i = 0; i < nq_; ++i)
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
        for (idx i = 0; i < nq_; ++i)
            if (!get_measured(i))
                result.emplace_back(i);

        return result;
    }

    /**
     * \brief Quantum circuit total gate count
     *
     * \return Total gate count
     */
    idx get_gate_count() const noexcept { return gates_.size(); }

    /**
     * \brief Quantum circuit total measurement count
     *
     * \return Total measurement count
     */
    idx get_measurement_count() const noexcept { return measurements_.size(); }

    /**
     * \brief Quantum circuit total count, i.e. the sum of gate count and
     * measurement count
     *
     * \return Total (gates + measurements) count
     */
    idx get_total_count() const noexcept { return step_cnt_; }
    // end getters

    /**
     * \brief Applies the single qudit gate \a U on single qudit \a i
     *
     * \param U Single qudit quantum gate
     * \param i Qudit index
     * \param name Optional gate name
     * \return Reference to the current instance
     */
    QCircuitDescription& gate(const cmat& U, idx i, std::string name = "") {
        // EXCEPTION CHECKS

        try {
            // check valid target
            if (i > nq_)
                throw exception::OutOfRange("qpp::QCircuitDescription::gate()");
            // check not measured before
            if (get_measured(i))
                throw exception::QuditAlreadyMeasured(
                    "qpp::QCircuitDescription::gate()");

            // check square matrix for the gate
            if (!internal::check_square_mat(U))
                throw exception::MatrixNotSquare(
                    "qpp::QCircuitDescription::gate()");
            // check correct dimension
            if (static_cast<idx>(U.rows()) != d_)
                throw exception::DimsMismatchMatrix(
                    "qpp::QCircuitDescription::gate()");
        } catch (qpp::exception::Exception& e) {
            std::cerr << "At STEP " << step_cnt_ << "\n";
            throw;
        }
        // END EXCEPTION CHECKS

        if (name == "")
            name = qpp::Gates::get_instance().get_name(U);
        gates_.emplace_back(GateType::SINGLE, U, std::vector<idx>{},
                            std::vector<idx>{i}, step_cnt_++, name);

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
    QCircuitDescription& gate(const cmat& U, idx i, idx j,
                              std::string name = "") {
        // EXCEPTION CHECKS

        try {
            // check valid target
            if (i > nq_ || j > nq_ || i == j)
                throw exception::OutOfRange("qpp::QCircuitDescription::gate()");
            if (get_measured(i) || get_measured(j))
                throw exception::QuditAlreadyMeasured(
                    "qpp::QCircuitDescription::gate()");

            // check square matrix for the gate
            if (!internal::check_square_mat(U))
                throw exception::MatrixNotSquare(
                    "qpp::QCircuitDescription::gate()");
            // check correct dimension
            if (static_cast<idx>(U.rows()) != d_ * d_)
                throw exception::DimsMismatchMatrix(
                    "qpp::QCircuitDescription::gate()");
        } catch (qpp::exception::Exception& e) {
            std::cerr << "At STEP " << step_cnt_ << "\n";
            throw;
        }
        // END EXCEPTION CHECKS

        if (name == "")
            name = qpp::Gates::get_instance().get_name(U);
        gates_.emplace_back(GateType::TWO, U, std::vector<idx>{},
                            std::vector<idx>{i, j}, step_cnt_++, name);

        return *this;
    }

    /**
     * \brief Applies the three qudit gate \a U on qudits \a i, \a j and \a k
     *
     * \param U Three qudit quantum gate
     * \param i Qudit index
     * \param j Qudit index
     * \param k Qudit index
     * \param name Optional gate name
     * \return Reference to the current instance
     */
    QCircuitDescription& gate(const cmat& U, idx i, idx j, idx k,
                              std::string name = "") {
        // EXCEPTION CHECKS

        try {
            // check valid target
            if (i > nq_ || j > nq_ || k > nq_ || (i == j) || (i == k) ||
                (j == k))
                throw exception::OutOfRange("qpp::QCircuitDescription::gate()");
            if (get_measured(i) || get_measured(j) || get_measured(k))
                throw exception::QuditAlreadyMeasured(
                    "qpp::QCircuitDescription::gate()");

            // check square matrix for the gate
            if (!internal::check_square_mat(U))
                throw exception::MatrixNotSquare(
                    "qpp::QCircuitDescription::gate()");
            // check correct dimension
            if (static_cast<idx>(U.rows()) != d_ * d_ * d_)
                throw exception::DimsMismatchMatrix(
                    "qpp::QCircuitDescription::gate()");
        } catch (qpp::exception::Exception& e) {
            std::cerr << "At STEP " << step_cnt_ << "\n";
            throw;
        }
        // END EXCEPTION CHECKS

        if (name == "")
            name = qpp::Gates::get_instance().get_name(U);
        gates_.emplace_back(GateType::THREE, U, std::vector<idx>{},
                            std::vector<idx>{i, j, k}, step_cnt_++, name);

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
    QCircuitDescription& gate_fan(const cmat& U, const std::vector<idx>& target,
                                  std::string name = "") {
        // EXCEPTION CHECKS

        idx n = static_cast<idx>(target.size());
        try {
            // check valid target
            if (n == 0)
                throw exception::ZeroSize(
                    "qpp::QCircuitDescription::gate_fan()");
            for (idx i = 0; i < n; ++i)
                if (target[i] > nq_)
                    throw exception::OutOfRange(
                        "qpp::QCircuitDescription::gate_fan()");
            // check no duplicates
            if (!internal::check_no_duplicates(target))
                throw exception::Duplicates(
                    "qpp::QCircuitDescription::gate_fan()");
            // check target was not measured before
            for (idx i = 0; i < n; ++i)
                if (get_measured(target[i]))
                    throw exception::QuditAlreadyMeasured(
                        "qpp::QCircuitDescription::gate_fan()");

            // check square matrix for the gate
            if (!internal::check_square_mat(U))
                throw exception::MatrixNotSquare(
                    "qpp::QCircuitDescription::gate_fan()");
            // check correct dimension
            if (static_cast<idx>(U.rows()) != d_)
                throw exception::DimsMismatchMatrix(
                    "qpp::QCircuitDescription::gate_fan()");
        } catch (qpp::exception::Exception& e) {
            std::cerr << "At STEP " << step_cnt_ << "\n";
            throw;
        }
        // END EXCEPTION CHECKS

        if (name == "")
            name = qpp::Gates::get_instance().get_name(U);
        gates_.emplace_back(GateType::FAN, U, std::vector<idx>{}, target,
                            step_cnt_++, name);

        return *this;
    }

    /**
     * \brief Applies the single qudit gate \a U on every remaining non-measured
     * qudit
     *
     * \param U Single qudit quantum gate
     * \param name Optional gate name
     * \return Reference to the current instance
     */
    QCircuitDescription& gate_fan(const cmat& U, std::string name = "") {
        // EXCEPTION CHECKS

        try {
            // check square matrix for the gate
            if (!internal::check_square_mat(U))
                throw exception::MatrixNotSquare(
                    "qpp::QCircuitDescription::gate_fan()");
            // check correct dimension
            if (static_cast<idx>(U.rows()) != d_)
                throw exception::DimsMismatchMatrix(
                    "qpp::QCircuitDescription::gate_fan()");
        } catch (qpp::exception::Exception& e) {
            std::cerr << "At STEP " << step_cnt_ << "\n";
            throw;
        }
        // END EXCEPTION CHECKS

        if (name == "")
            name = qpp::Gates::get_instance().get_name(U);
        gates_.emplace_back(GateType::FAN, U, std::vector<idx>{},
                            get_non_measured(), step_cnt_++, name);

        return *this;
    }

    /**
     * \brief Jointly applies the custom multiple qudit gate \a U on the qudit
     * indexes specified by \a target
     *
     * \param U Multiple qudit quantum gate
     * \param target Subsystem indexes where the gate \a U is applied
     * \param name Optional gate name
     *
     * \return Reference to the current instance
     */
    QCircuitDescription& gate_custom(const cmat& U,
                                     const std::vector<idx>& target,
                                     std::string name = "") {
        // EXCEPTION CHECKS

        idx n = static_cast<idx>(target.size());
        idx D = static_cast<idx>(std::llround(std::pow(d_, n)));

        try {
            // check valid target
            if (n == 0)
                throw exception::ZeroSize(
                    "qpp::QCircuitDescription::gate_custom()");
            for (idx i = 0; i < n; ++i)
                if (target[i] > nq_)
                    throw exception::OutOfRange(
                        "qpp::QCircuitDescription::gate()");
            // check no duplicates
            if (!internal::check_no_duplicates(target))
                throw exception::Duplicates(
                    "qpp::QCircuitDescription::gate_custom()");
            // check target was not measured before
            for (idx i = 0; i < n; ++i)
                if (get_measured(target[i]))
                    throw exception::QuditAlreadyMeasured(
                        "qpp::QCircuitDescription::gate_custom()");

            // check square matrix for the gate
            if (!internal::check_square_mat(U))
                throw exception::MatrixNotSquare(
                    "qpp::QCircuitDescription::gate_custom()");
            // check correct dimension
            if (static_cast<idx>(U.rows()) != D)
                throw exception::DimsMismatchMatrix(
                    "qpp::QCircuitDescription::gate_custom()");
        } catch (qpp::exception::Exception& e) {
            std::cerr << "At STEP " << step_cnt_ << "\n";
            throw;
        }
        // END EXCEPTION CHECKS

        if (name == "")
            name = qpp::Gates::get_instance().get_name(U);
        gates_.emplace_back(GateType::CUSTOM, U, std::vector<idx>{}, target,
                            step_cnt_++, name);

        return *this;
    }

    // QFT
    /**
     * \brief Applies the quantum Fourier transform (as a series of gates) on
     * the qudit indexes specified by \a target
     *
     * \param target Subsystem indexes where the quantum Fourier transform is
     * applied
     * \param swap Swaps the qubits at the end (true by default)
     * \return Reference to the current instance
     */
    QCircuitDescription& QFT(const std::vector<idx>& target,
                             bool swap QPP_UNUSED_ = true) {
        // EXCEPTION CHECKS

        try {
            throw exception::NotImplemented("qpp::QCircuitDescription::QFT()");
        } catch (qpp::exception::Exception& e) {
            std::cerr << "At STEP " << step_cnt_ << "\n";
            throw;
        }
        // END EXCEPTION CHECKS

        gates_.emplace_back(GateType::QFT, cmat{}, std::vector<idx>{}, target,
                            step_cnt_++, "QFT");

        return *this;
    }

    // TFQ
    /**
     * \brief Applies the inverse quantum Fourier transform (as a series of
     * gates) on the qudit indexes specified by \a target
     *
     * \param target Subsystem indexes where the inverse quantum Fourier
     * transform is applied
     * \param swap Swaps the qubits at the end (true by default)
     * \return Reference to the current instance
     */
    QCircuitDescription& TFQ(const std::vector<idx>& target,
                             bool swap QPP_UNUSED_ = true) {
        // EXCEPTION CHECKS

        try {
            throw exception::NotImplemented("qpp::QCircuitDescription::TFQ()");
        } catch (qpp::exception::Exception& e) {
            std::cerr << "At STEP " << step_cnt_ << "\n";
            throw;
        }
        // END EXCEPTION CHECKS
        gates_.emplace_back(GateType::TFQ, cmat{}, std::vector<idx>{}, target,
                            step_cnt_++, "TFQ");

        return *this;
    }

    // single ctrl single target
    /**
     * \brief Applies the single qudit controlled gate \a U with control qudit
     * \a ctrl and target qudit \a target
     *
     * \param U Single qudit quantum gate
     * \param ctrl Control qudit index
     * \param target Target qudit index
     * \param name Optional gate name
     * \return Reference to the current instance
     */
    QCircuitDescription& CTRL(const cmat& U, idx ctrl, idx target,
                              std::string name = "") {
        // EXCEPTION CHECKS

        try {
            // check valid ctrl and target
            if (ctrl > nq_ || target > nq_ || ctrl == target)
                throw exception::OutOfRange("qpp::QCircuitDescription::CTRL()");
            if (get_measured(ctrl) || get_measured(target))
                throw exception::QuditAlreadyMeasured(
                    "qpp::QCircuitDescription::CTRL()");

            // check square matrix for the gate
            if (!internal::check_square_mat(U))
                throw exception::MatrixNotSquare(
                    "qpp::QCircuitDescription::CTRL()");
            // check correct dimension
            if (static_cast<idx>(U.rows()) != d_)
                throw exception::DimsMismatchMatrix(
                    "qpp::QCircuitDescription::CTRL()");
        } catch (qpp::exception::Exception& e) {
            std::cerr << "At STEP " << step_cnt_ << "\n";
            throw;
        }
        // END EXCEPTION CHECKS

        if (name == "")
            name = qpp::Gates::get_instance().get_name(U);
        gates_.emplace_back(GateType::SINGLE_CTRL_SINGLE_TARGET, U,
                            std::vector<idx>{ctrl}, std::vector<idx>{target},
                            step_cnt_++, name);

        return *this;
    }

    // single ctrl multiple target
    /**
     * \brief Applies the single qudit controlled gate \a U with control qudit
     * \a ctrl on every qudit listed in \a target
     *
     * \param U Single qudit quantum gate
     * \param ctrl Control qudit index
     * \param target Target qudit indexes; the gate \a U is applied on every one
     * of them depending on the values of the control qudits
     * \param name Optional gate name
     * \return Reference to the current instance
     */
    QCircuitDescription& CTRL(const cmat& U, idx ctrl,
                              const std::vector<idx>& target,
                              std::string name = "") {
        // EXCEPTION CHECKS

        idx n = static_cast<idx>(target.size());

        try {
            // check valid ctrl
            if (ctrl > nq_)
                throw exception::OutOfRange("qpp::QCircuitDescription::CTRL()");
            if (get_measured(ctrl))
                throw exception::QuditAlreadyMeasured(
                    "qpp::QCircuitDescription::CTRL()");

            // check valid target
            if (n == 0)
                throw exception::ZeroSize("qpp::QCircuitDescription::CTRL()");
            for (idx i = 0; i < n; ++i)
                if (target[i] > nq_)
                    throw exception::OutOfRange(
                        "qpp::QCircuitDescription::CTRL()");
            // check no duplicates target
            if (!internal::check_no_duplicates(target))
                throw exception::Duplicates("qpp::QCircuitDescription::CTRL()");
            // check target was not measured before and that ctrl is not part of
            // target
            for (idx i = 0; i < n; ++i) {
                if (get_measured(target[i]))
                    throw exception::QuditAlreadyMeasured(
                        "qpp::QCircuitDescription::CTRL()");
                if (ctrl == target[i])
                    throw exception::OutOfRange(
                        "qpp::QCircuitDescription::CTRL()");
            }

            // check square matrix for the gate
            if (!internal::check_square_mat(U))
                throw exception::MatrixNotSquare(
                    "qpp::QCircuitDescription::CTRL()");
            // check correct dimension
            if (static_cast<idx>(U.rows()) != d_)
                throw exception::DimsMismatchMatrix(
                    "qpp::QCircuitDescription::CTRL()");
        } catch (qpp::exception::Exception& e) {
            std::cerr << "At STEP " << step_cnt_ << "\n";
            throw;
        }
        // END EXCEPTION CHECKS

        if (name == "")
            name = qpp::Gates::get_instance().get_name(U);
        gates_.emplace_back(GateType::SINGLE_CTRL_MULTIPLE_TARGET, U,
                            std::vector<idx>{ctrl}, target, step_cnt_++, name);

        return *this;
    }

    // multiple ctrl single target
    /**
     * \brief Applies the single qudit controlled gate \a U with multiple
     * control qudits listed in \a ctrl on the target qudit \a target
     *
     * \param U Single qudit quantum gate
     * \param ctrl Control qudits indexes
     * \param target Target qudit index
     * \param name Optional gate name
     * \return Reference to the current instance
     */
    QCircuitDescription& CTRL(const cmat& U, const std::vector<idx>& ctrl,
                              idx target, std::string name = "") {
        // EXCEPTION CHECKS

        try {
            // check valid ctrl
            for (idx i = 0; i < static_cast<idx>(ctrl.size()); ++i)
                if (ctrl[i] > nq_)
                    throw exception::OutOfRange(
                        "qpp::QCircuitDescription::CTRL()");
            // check no duplicates ctrl
            if (!internal::check_no_duplicates(ctrl))
                throw exception::Duplicates("qpp::QCircuitDescription::CTRL()");
            // check ctrl was not measured before and that ctrl is not part of
            // target
            for (idx i = 0; i < static_cast<idx>(ctrl.size()); ++i) {
                if (get_measured(ctrl[i]))
                    throw exception::QuditAlreadyMeasured(
                        "qpp::QCircuitDescription::CTRL()");
                if (target == ctrl[i])
                    throw exception::OutOfRange(
                        "qpp::QCircuitDescription::CTRL()");
            }

            // check valid target
            if (target > nq_)
                throw exception::OutOfRange("qpp::QCircuitDescription::CTRL()");
            if (get_measured(target))
                throw exception::QuditAlreadyMeasured(
                    "qpp::QCircuitDescription::CTRL()");

            // check square matrix for the gate
            if (!internal::check_square_mat(U))
                throw exception::MatrixNotSquare(
                    "qpp::QCircuitDescription::CTRL()");
            // check correct dimension
            if (static_cast<idx>(U.rows()) != d_)
                throw exception::DimsMismatchMatrix(
                    "qpp::QCircuitDescription::CTRL()");
        } catch (qpp::exception::Exception& e) {
            std::cerr << "At STEP " << step_cnt_ << "\n";
            throw;
        }
        // END EXCEPTION CHECKS

        if (name == "")
            name = qpp::Gates::get_instance().get_name(U);
        gates_.emplace_back(GateType::MULTIPLE_CTRL_SINGLE_TARGET, U, ctrl,
                            std::vector<idx>{target}, step_cnt_++, name);

        return *this;
    }

    // multiple ctrl multiple target
    // FIXME
    /**
     * \brief Applies the single qudit controlled gate \a U with multiple
     * control qudits listed in \a ctrl on every qudit listed in \a target
     *
     * \param U Single qudit quantum gate
     * \param ctrl Control qudits indexes
     * \param target Target qudit indexes; the gate \a U is applied on every one
     * of them depending on the values of the control qudits
     * \param name Optional gate name
     * \return Reference to the current instance
     */
    QCircuitDescription& CTRL(const cmat& U, const std::vector<idx>& ctrl,
                              const std::vector<idx>& target,
                              std::string name = "") {
        // EXCEPTION CHECKS

        idx n = static_cast<idx>(target.size());

        try {
            // check valid ctrl
            for (idx i = 0; i < static_cast<idx>(ctrl.size()); ++i)
                if (ctrl[i] > nq_)
                    throw exception::OutOfRange(
                        "qpp::QCircuitDescription::CTRL()");
            // check no duplicates ctrl
            if (!internal::check_no_duplicates(ctrl))
                throw exception::Duplicates("qpp::QCircuitDescription::CTRL()");
            // check ctrl was not measured before
            for (idx i = 0; i < static_cast<idx>(ctrl.size()); ++i) {
                if (get_measured(ctrl[i]))
                    throw exception::QuditAlreadyMeasured(
                        "qpp::QCircuitDescription::CTRL()");
            }

            // check valid target
            if (n == 0)
                throw exception::ZeroSize("qpp::QCircuitDescription::CTRL()");
            for (idx i = 0; i < n; ++i)
                if (target[i] > nq_)
                    throw exception::OutOfRange(
                        "qpp::QCircuitDescription::CTRL()");
            // check no duplicates target
            if (!internal::check_no_duplicates(target))
                throw exception::Duplicates("qpp::QCircuitDescription::CTRL()");
            // check target was not measured before
            for (idx i = 0; i < n; ++i) {
                if (get_measured(target[i]))
                    throw exception::QuditAlreadyMeasured(
                        "qpp::QCircuitDescription::CTRL()");
            }

            // check ctrl and target don't share common elements
            for (idx i = 0; i < static_cast<idx>(ctrl.size()); ++i)
                for (idx j = 0; j < n; ++j)
                    if (ctrl[i] == target[j])
                        throw exception::OutOfRange(
                            "qpp::QCircuitDescription::CTRL()");

            // check square matrix for the gate
            if (!internal::check_square_mat(U))
                throw exception::MatrixNotSquare(
                    "qpp::QCircuitDescription::CTRL()");
            // check correct dimension
            if (static_cast<idx>(U.rows()) != d_)
                throw exception::DimsMismatchMatrix(
                    "qpp::QCircuitDescription::CTRL()");
        } catch (qpp::exception::Exception& e) {
            std::cerr << "At STEP " << step_cnt_ << "\n";
            throw;
        }
        // END EXCEPTION CHECKS

        if (name == "")
            name = qpp::Gates::get_instance().get_name(U);
        gates_.emplace_back(GateType::MULTIPLE_CTRL_MULTIPLE_TARGET, U, ctrl,
                            std::vector<idx>{target}, step_cnt_++, name);

        return *this;
    }

    // custom multiple control composed target
    /**
     * \brief Jointly applies the custom multiple-qudit controlled gate \a U
     * with multiple control qudits listed in \a ctrl on the qudit indexes
     * specified by \a target
     *
     * \param U Multiple-qudit quantum gate
     * \param ctrl Control qudits indexes
     * \param target Target qudit indexes where the gate \a U is applied
     * depending on the values of the control qudits
     * \param name Optional gate name
     * \return Reference to the current instance
     */
    QCircuitDescription& CTRL_custom(const cmat& U,
                                     const std::vector<idx>& ctrl,
                                     const std::vector<idx>& target,
                                     std::string name = "") {
        // EXCEPTION CHECKS

        idx n = static_cast<idx>(target.size());
        idx D = static_cast<idx>(std::llround(std::pow(d_, n)));

        try {
            // check valid ctrl
            for (idx i = 0; i < static_cast<idx>(ctrl.size()); ++i)
                if (ctrl[i] > nq_)
                    throw exception::OutOfRange(
                        "qpp::QCircuitDescription::CTRL_custom()");
            // check no duplicates ctrl
            if (!internal::check_no_duplicates(ctrl))
                throw exception::Duplicates(
                    "qpp::QCircuitDescription::CTRL_custom()");
            // check ctrl was not measured before
            for (idx i = 0; i < static_cast<idx>(ctrl.size()); ++i) {
                if (get_measured(ctrl[i]))
                    throw exception::QuditAlreadyMeasured(
                        "qpp::QCircuitDescription::CTRL_custom()");
            }

            // check valid target
            if (n == 0)
                throw exception::ZeroSize(
                    "qpp::QCircuitDescription::CTRL_custom()");
            for (idx i = 0; i < n; ++i)
                if (target[i] > nq_)
                    throw exception::OutOfRange(
                        "qpp::QCircuitDescription::CTRL_custom()");
            // check no duplicates target
            if (!internal::check_no_duplicates(target))
                throw exception::Duplicates(
                    "qpp::QCircuitDescription::CTRL_custom()");
            // check target was not measured before
            for (idx i = 0; i < n; ++i) {
                if (get_measured(target[i]))
                    throw exception::QuditAlreadyMeasured(
                        "qpp::QCircuitDescription::CTRL_custom()");
            }

            // check ctrl and target don't share common elements
            for (idx i = 0; i < static_cast<idx>(ctrl.size()); ++i)
                for (idx j = 0; j < static_cast<idx>(target.size()); ++j)
                    if (ctrl[i] == target[j])
                        throw exception::OutOfRange(
                            "qpp::QCircuitDescription::CTRL_custom()");

            // check square matrix for the gate
            if (!internal::check_square_mat(U))
                throw exception::MatrixNotSquare(
                    "qpp::QCircuitDescription::CTRL_custom()");
            // check correct dimension
            if (static_cast<idx>(U.rows()) != D)
                throw exception::DimsMismatchMatrix(
                    "qpp::QCircuitDescription::CTRL_custom()");
        } catch (qpp::exception::Exception& e) {
            std::cerr << "At STEP " << step_cnt_ << "\n";
            throw;
        }
        // END EXCEPTION CHECKS

        if (name == "")
            name = qpp::Gates::get_instance().get_name(U);
        gates_.emplace_back(GateType::CUSTOM_CTRL, U, ctrl, target, step_cnt_++,
                            name);

        return *this;
    }

    // single ctrl single target
    // FIXME, use the corresponding dits
    /**
     * \brief Applies the single qubit controlled gate \a U with classical
     * control dit \a ctrl and target qudit \a target
     *
     * \param U Single qudit quantum gate
     * \param ctrl_dit Classical control dit index
     * \param target Target qudit index
     * \param name Optional gate name
     * \return Reference to the current instance
     */
    QCircuitDescription& cCTRL(const cmat& U, idx ctrl_dit, idx target,
                               std::string name = "") {
        // EXCEPTION CHECKS

        try {
            // check valid ctrl_dit and target
            if (ctrl_dit > nc_ || target > nq_)
                throw exception::OutOfRange(
                    "qpp::QCircuitDescription::cCTRL()");
            if (get_measured(target))
                throw exception::QuditAlreadyMeasured(
                    "qpp::QCircuitDescription::cCTRL()");

            // check square matrix for the gate
            if (!internal::check_square_mat(U))
                throw exception::MatrixNotSquare(
                    "qpp::QCircuitDescription::cCTRL()");
            // check correct dimension
            if (static_cast<idx>(U.rows()) != d_)
                throw exception::DimsMismatchMatrix(
                    "qpp::QCircuitDescription::cCTRL()");
        } catch (qpp::exception::Exception& e) {
            std::cerr << "At STEP " << step_cnt_ << "\n";
            throw;
        }
        // END EXCEPTION CHECKS

        if (name == "")
            name = qpp::Gates::get_instance().get_name(U);
        gates_.emplace_back(GateType::SINGLE_cCTRL_SINGLE_TARGET, U,
                            std::vector<idx>{ctrl_dit},
                            std::vector<idx>{target}, step_cnt_++, name);

        return *this;
    }

    // single ctrl multiple targets
    /**
     * \brief Applies the single qudit controlled gate \a U with classical
     * control dit \a ctrl on every qudit listed in \a target
     *
     * \param U Single qudit quantum gate
     * \param ctrl_dit Classical control dit index
     * \param target Target qudit indexes; the gate \a U is applied on every one
     * of them depending on the values of the classical control dits
     * \param name Optional gate name
     * \return Reference to the current instance
     */
    QCircuitDescription& cCTRL(const cmat& U, idx ctrl_dit,
                               const std::vector<idx>& target,
                               std::string name = "") {
        // EXCEPTION CHECKS

        idx n = static_cast<idx>(target.size());

        try {
            // check valid ctrl_dit
            if (ctrl_dit > nc_)
                throw exception::OutOfRange(
                    "qpp::QCircuitDescription::cCTRL()");

            // check valid target
            if (n == 0)
                throw exception::ZeroSize("qpp::QCircuitDescription::cCTRL()");
            for (idx i = 0; i < n; ++i)
                if (target[i] > nq_)
                    throw exception::OutOfRange(
                        "qpp::QCircuitDescription::cCTRL()");
            // check no duplicates target
            if (!internal::check_no_duplicates(target))
                throw exception::Duplicates(
                    "qpp::QCircuitDescription::cCTRL()");
            // check target was not measured before
            for (idx i = 0; i < n; ++i) {
                if (get_measured(target[i]))
                    throw exception::QuditAlreadyMeasured(
                        "qpp::QCircuitDescription::cCTRL()");
            }

            // check square matrix for the gate
            if (!internal::check_square_mat(U))
                throw exception::MatrixNotSquare(
                    "qpp::QCircuitDescription::cCTRL()");
            // check correct dimension
            if (static_cast<idx>(U.rows()) != d_)
                throw exception::DimsMismatchMatrix(
                    "qpp::QCircuitDescription::cCTRL()");
        } catch (qpp::exception::Exception& e) {
            std::cerr << "At STEP " << step_cnt_ << "\n";
            throw;
        }
        // END EXCEPTION CHECKS

        if (name == "")
            name = qpp::Gates::get_instance().get_name(U);
        gates_.emplace_back(GateType::SINGLE_cCTRL_MULTIPLE_TARGET, U,
                            std::vector<idx>{ctrl_dit}, target, step_cnt_++,
                            name);

        return *this;
    }

    // multiple ctrl single target
    /**
     * \brief Applies the single qudit controlled gate \a U with multiple
     * classical control dits listed in \a ctrl on the target qudit \a target
     *
     * \param U Single qudit quantum gate
     * \param ctrl_dits Classical control dits indexes
     * \param target Target qudit index
     * \param name Optional gate name
     * \return Reference to the current instance
     */
    QCircuitDescription& cCTRL(const cmat& U, const std::vector<idx>& ctrl_dits,
                               idx target, std::string name = "") {
        // EXCEPTION CHECKS

        try {
            // check valid ctrl_dits
            for (idx i = 0; i < ctrl_dits.size(); ++i)
                if (ctrl_dits[i] > nc_)
                    throw exception::OutOfRange(
                        "qpp::QCircuitDescription::cCTRL()");
            // check no duplicates ctrl_dits
            if (!internal::check_no_duplicates(ctrl_dits))
                throw exception::Duplicates(
                    "qpp::QCircuitDescription::cCTRL()");

            // check valid target
            if (target > nq_)
                throw exception::OutOfRange(
                    "qpp::QCircuitDescription::cCTRL()");
            // check target was not measured before
            if (get_measured(target))
                throw exception::QuditAlreadyMeasured(
                    "qpp::QCircuitDescription::cCTRL()");

            // check square matrix for the gate
            if (!internal::check_square_mat(U))
                throw exception::MatrixNotSquare(
                    "qpp::QCircuitDescription::cCTRL()");
            // check correct dimension
            if (static_cast<idx>(U.rows()) != d_)
                throw exception::DimsMismatchMatrix(
                    "qpp::QCircuitDescription::cCTRL()");
        } catch (qpp::exception::Exception& e) {
            std::cerr << "At STEP " << step_cnt_ << "\n";
            throw;
        }
        // END EXCEPTION CHECKS

        if (name == "")
            name = qpp::Gates::get_instance().get_name(U);
        gates_.emplace_back(GateType::MULTIPLE_cCTRL_SINGLE_TARGET, U,
                            ctrl_dits, std::vector<idx>{target}, step_cnt_++,
                            name);

        return *this;
    }

    // multiple ctrl multiple targets
    /**
     * \brief Applies the single qudit controlled gate \a U with multiple
     * classical control dits listed in \a ctrl on every qudit listed in
     * \a target
     *
     * \param U Single qudit quantum gate
     * \param ctrl_dits Classical control dits indexes
     * \param target Target qudit indexes; the gate \a U is applied on every one
     * of them depending on the values of the classical control dits
     * \param name Optional gate name
     * \return Reference to the current instance
     */
    QCircuitDescription& cCTRL(const cmat& U, const std::vector<idx>& ctrl_dits,
                               const std::vector<idx>& target,
                               std::string name = "") {
        // EXCEPTION CHECKS

        idx n = static_cast<idx>(target.size());

        try {
            // check valid ctrl_dits
            for (idx i = 0; i < ctrl_dits.size(); ++i)
                if (ctrl_dits[i] > nc_)
                    throw exception::OutOfRange(
                        "qpp::QCircuitDescription::cCTRL()");
            // check no duplicates ctrl_dits
            if (!internal::check_no_duplicates(ctrl_dits))
                throw exception::Duplicates(
                    "qpp::QCircuitDescription::cCTRL()");

            // check valid target
            if (n == 0)
                throw exception::ZeroSize("qpp::QCircuitDescription::cCTRL()");
            for (idx i = 0; i < n; ++i)
                if (target[i] > nq_)
                    throw exception::OutOfRange(
                        "qpp::QCircuitDescription::cCTRL()");
            // check no duplicates target
            if (!internal::check_no_duplicates(target))
                throw exception::Duplicates(
                    "qpp::QCircuitDescription::cCTRL()");
            // check target was not measured before
            for (idx i = 0; i < n; ++i) {
                if (get_measured(target[i]))
                    throw exception::QuditAlreadyMeasured(
                        "qpp::QCircuitDescription::cCTRL()");
            }

            // check square matrix for the gate
            if (!internal::check_square_mat(U))
                throw exception::MatrixNotSquare(
                    "qpp::QCircuitDescription::cCTRL()");
            // check correct dimension
            if (static_cast<idx>(U.rows()) != d_)
                throw exception::DimsMismatchMatrix(
                    "qpp::QCircuitDescription::cCTRL()");
        } catch (qpp::exception::Exception& e) {
            std::cerr << "At STEP " << step_cnt_ << "\n";
            throw;
        }
        // END EXCEPTION CHECKS

        if (name == "")
            name = qpp::Gates::get_instance().get_name(U);
        gates_.emplace_back(GateType::MULTIPLE_cCTRL_MULTIPLE_TARGET, U,
                            ctrl_dits, std::vector<idx>{target}, step_cnt_++,
                            name);

        return *this;
    }

    //  custom controlled gate with multiple controls and multiple targets
    /**
     * \brief Jointly applies the custom multiple-qudit controlled gate \a U
     * with multiple classical control dits listed in \a ctrl on the qudit
     * indexes specified by \a target
     *
     * \param U Multiple-qudit quantum gate
     * \param ctrl_dits Classical control dits indexes
     * \param target Target qudit indexes where the gate \a U is applied
     * depending on the values of the classical control dits
     * \param name Optional gate name
     * \return Reference to the current instance
     */
    QCircuitDescription& cCTRL_custom(const cmat& U,
                                      const std::vector<idx>& ctrl_dits,
                                      const std::vector<idx>& target,
                                      std::string name = "") {
        // EXCEPTION CHECKS

        idx n = static_cast<idx>(target.size());
        idx D = static_cast<idx>(std::llround(std::pow(d_, n)));

        try {
            // check valid ctrl_dits
            for (idx i = 0; i < ctrl_dits.size(); ++i)
                if (ctrl_dits[i] > nc_)
                    throw exception::OutOfRange(
                        "qpp::QCircuitDescription::cCTRL_custom()");
            // check no duplicates ctrl_dits
            if (!internal::check_no_duplicates(ctrl_dits))
                throw exception::Duplicates(
                    "qpp::QCircuitDescription::cCTRL_custom()");

            // check valid target
            if (n == 0)
                throw exception::ZeroSize(
                    "qpp::QCircuitDescription::cCTRL_custom()");
            for (idx i = 0; i < n; ++i)
                if (target[i] > nq_)
                    throw exception::OutOfRange(
                        "qpp::QCircuitDescription::cCTRL_custom()");
            // check no duplicates target
            if (!internal::check_no_duplicates(target))
                throw exception::Duplicates(
                    "qpp::QCircuitDescription::cCTRL_custom()");
            // check target was not measured before
            for (idx i = 0; i < n; ++i) {
                if (get_measured(target[i]))
                    throw exception::QuditAlreadyMeasured(
                        "qpp::QCircuitDescription::cCTRL_custom()");
            }

            // check square matrix for the gate
            if (!internal::check_square_mat(U))
                throw exception::MatrixNotSquare(
                    "qpp::QCircuitDescription::cCTRL_custom()");
            // check correct dimension
            if (static_cast<idx>(U.rows()) != D)
                throw exception::DimsMismatchMatrix(
                    "qpp::QCircuitDescription::cCTRL_custom()");
        } catch (qpp::exception::Exception& e) {
            std::cerr << "At STEP " << step_cnt_ << "\n";
            throw;
        }
        // END EXCEPTION CHECKS

        if (name == "")
            name = qpp::Gates::get_instance().get_name(U);
        gates_.emplace_back(GateType::CUSTOM_cCTRL, U, ctrl_dits, target,
                            step_cnt_++, name);

        return *this;
    }

    // Z measurement of single qudit
    /**
     * \brief Measurement of single qudit in the computational basis (Z-basis)
     *
     * \param i Qudit index
     * \param c_reg Classical register where the value of the measurement is
     * being stored
     * \param name Optional measurement name, default is "Measure Z"
     * \return Reference to the current instance
     */
    QCircuitDescription& measureZ(idx i, idx c_reg, std::string name = "") {
        // EXCEPTION CHECKS

        try {
            // measuring non-existing qudit
            if (i >= nq_)
                throw exception::OutOfRange(
                    "qpp::QCircuitDescription::measureZ()");
            // trying to put the result into an non-existing classical slot
            if (c_reg >= nc_)
                throw exception::OutOfRange(
                    "qpp::QCircuitDescription::measureZ()");
            // qudit was measured before
            if (get_measured(i))
                throw exception::QuditAlreadyMeasured(
                    "qpp:QCircuitDescription::measureZ");
        } catch (qpp::exception::Exception& e) {
            std::cerr << "At STEP " << step_cnt_ << "\n";
            throw;
        }
        // END EXCEPTION CHECKS

        if (name == "")
            name = "Measure Z";
        measured_[i] = true;
        measurements_.emplace_back(MeasureType::MEASURE_Z, std::vector<cmat>{},
                                   std::vector<idx>{i}, c_reg, step_cnt_++,
                                   name);
        measurement_steps_.emplace_back(gates_.size());

        return *this;
    }

    // measurement of single qudit in the orthonormal basis or rank-1 projectors
    // specified by the columns of matrix V
    /**
     * \brief Measurement of single qudit in the orthonormal basis or rank-1
     * projectors specified by the columns of matrix \a V
     *
     * \param V Orthonormal basis or rank-1 projectors specified by the
     * columns of matrix V
     * \param i Qudit index
     * \param c_reg Classical register where the value of the measurement is
     * stored
     * \param name Optional measurement name
     * \return Reference to the current instance
     */
    QCircuitDescription& measureV(const cmat& V, idx i, idx c_reg,
                                  std::string name = "") {
        // EXCEPTION CHECKS

        try {
            // measuring non-existing qudit
            if (i >= nq_)
                throw exception::OutOfRange(
                    "qpp::QCircuitDescription::measureV()");
            // trying to put the result into an non-existing classical slot
            if (c_reg >= nc_)
                throw exception::OutOfRange(
                    "qpp::QCircuitDescription::measureV()");
            // qudit was measured before
            if (get_measured(i))
                throw exception::QuditAlreadyMeasured(
                    "qpp:QCircuitDescription::measureV");
        } catch (qpp::exception::Exception& e) {
            std::cerr << "At STEP " << step_cnt_ << "\n";
            throw;
        }
        // END EXCEPTION CHECKS

        if (name == "")
            name = qpp::Gates::get_instance().get_name(V);
        measured_[i] = true;
        measurements_.emplace_back(MeasureType::MEASURE_V, std::vector<cmat>{V},
                                   std::vector<idx>{i}, c_reg, step_cnt_++,
                                   name);
        measurement_steps_.emplace_back(gates_.size());

        return *this;
    }

    // measurement of multiple qudits in the orthonormal basis or rank-1
    // projectors specified by the columns of matrix V
    /**
     * \brief Joint measurement of multiple qudits in the orthonormal basis or
     * rank-1 projectors specified by the columns of matrix \a V
     *
     * \param V Orthonormal basis or rank-1 projectors specified by the
     * columns of matrix V
     * \param target Target qudit indexes that are jointly measured
     * \param c_reg Classical register where the value of the measurement is
     * stored
     * \param name Optional measurement name
     * \return Reference to the current instance
     */
    QCircuitDescription& measureV(const cmat& V, const std::vector<idx>& target,
                                  idx c_reg, std::string name = "") {
        // EXCEPTION CHECKS

        try {
            // measuring non-existing qudit
            for (auto&& i : target)
                if (i >= nq_)
                    throw exception::OutOfRange(
                        "qpp::QCircuitDescription::measureV()");
            // trying to put the result into an non-existing classical slot
            if (c_reg >= nc_)
                throw exception::OutOfRange(
                    "qpp::QCircuitDescription::measureV()");
            // qudit was measured before
            for (auto&& i : target)
                if (get_measured(i))
                    throw exception::QuditAlreadyMeasured(
                        "qpp::QCircuitDescription::measureV");
        } catch (qpp::exception::Exception& e) {
            std::cerr << "At STEP " << step_cnt_ << "\n";
            throw;
        }
        // END EXCEPTION CHECKS

        if (name == "")
            name = qpp::Gates::get_instance().get_name(V);
        for (auto&& i : target)
            measured_[i] = true;
        measurements_.emplace_back(MeasureType::MEASURE_V_MANY,
                                   std::vector<cmat>{V}, target, c_reg,
                                   step_cnt_++, name);
        measurement_steps_.emplace_back(gates_.size());

        return *this;
    }

    /**
     * \brief qpp::IDisplay::display() override
     *
     * Writes to the output stream a textual representation of the
     * quantum circuit description
     *
     * \param os Output stream
     * \return Output stream
     */
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
    void set_measured_(idx i) {
        if (was_measured(i))
            throw exception::QuditAlreadyMeasured(
                "qpp::QCircuit::set_measured_()");
        subsys_[i] = idx_infty; // set qudit i to measured state
        for (idx m = i; m < qcd_.get_nq(); ++m) {
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
        : qcd_{qcd}, psi_{st.zero(qcd.get_nq(), qcd.get_d())},
          dits_(qcd.get_nc(), 0), probs_(qcd.get_nc(), 0),
          subsys_(qcd.get_nq(), 0), m_ip_{0}, q_ip_{0}, ip_{0} {
        std::iota(std::begin(subsys_), std::end(subsys_), 0);
    }

    // getters

    ket get_psi() const { return psi_; }

    std::vector<idx> get_dits() const { return dits_; }

    idx get_dit(idx i) const {
        if (i > qcd_.get_nc())
            throw exception::OutOfRange("qpp::QCircuit::get_dit()");

        return dits_[i];
    }

    std::vector<double> get_probs() const { return probs_; }

    idx was_measured(idx i) const { return subsys_[i] == idx_infty; }

    std::vector<idx> get_measured() const {
        std::vector<idx> result;
        for (idx i = 0; i < qcd_.get_nq(); ++i)
            if (was_measured(i))
                result.emplace_back(i);

        return result;
    }

    // qudits that were not (yet) measured
    std::vector<idx> get_not_measured() const {
        std::vector<idx> result;
        for (idx i = 0; i < qcd_.get_nq(); ++i)
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

    // setters
    QCircuit& set_dit(idx i, idx val) {
        if (i > qcd_.get_nc())
            throw exception::OutOfRange("qpp::QCircuit::set_dit()");
        dits_[i] = val;

        return *this;
    }
    // end setters

    void reset() {
        psi_ = st.zero(qcd_.get_nq(), qcd_.get_d());
        dits_ = std::vector<idx>(qcd_.get_nc(), 0);
        probs_ = std::vector<double>(qcd_.get_nc(), 0);
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

        for (; q_ip_ <= qcd_.get_gates().size(); ++q_ip_) {
            // check for measurements
            if (m_ip_ < qcd_.get_measurements().size()) {
                idx m_step = qcd_.get_measurement_steps()[m_ip_];
                // we have a measurement at step i
                if (m_step == q_ip_) {
                    while (qcd_.get_measurement_steps()[m_ip_] == m_step) {
                        if (ip_ - saved_ip_ == no_steps)
                            return;

                        if (verbose) {
                            std::cout << std::left;
                            std::cout
                                << std::setw(8)
                                << qcd_.get_measurements()[m_ip_].step_no_;
                            std::cout << std::right;
                            std::cout << "|> " << qcd_.get_measurements()[m_ip_]
                                      << '\n';
                        }

                        std::vector<idx> target_rel_pos = get_relative_pos_(
                            qcd_.get_measurements()[m_ip_].target_);

                        std::vector<idx> resZ;
                        double probZ;

                        idx mres = 0;
                        std::vector<double> probs;
                        std::vector<cmat> states;

                        switch (
                            qcd_.get_measurements()[m_ip_].measurement_type_) {
                        case QCircuitDescription::MeasureType::NONE:
                            break;
                        case QCircuitDescription::MeasureType::MEASURE_Z:
                            std::tie(resZ, probZ, psi_) =
                                measure_seq(psi_, target_rel_pos, qcd_.get_d());
                            dits_[qcd_.get_measurements()[m_ip_].c_reg_] =
                                resZ[0];
                            probs_[qcd_.get_measurements()[m_ip_].c_reg_] =
                                probZ;
                            set_measured_(
                                qcd_.get_measurements()[m_ip_].target_[0]);
                            break;
                        case QCircuitDescription::MeasureType::MEASURE_V:
                            std::tie(mres, probs, states) = measure(
                                psi_, qcd_.get_measurements()[m_ip_].mats_[0],
                                target_rel_pos, qcd_.get_d());
                            psi_ = states[mres];
                            dits_[qcd_.get_measurements()[m_ip_].c_reg_] = mres;
                            probs_[qcd_.get_measurements()[m_ip_].c_reg_] =
                                probs[mres];
                            set_measured_(
                                qcd_.get_measurements()[m_ip_].target_[0]);
                            break;
                        case QCircuitDescription::MeasureType::MEASURE_V_MANY:
                            std::tie(mres, probs, states) = measure(
                                psi_, qcd_.get_measurements()[m_ip_].mats_[0],
                                target_rel_pos, qcd_.get_d());
                            psi_ = states[mres];
                            dits_[qcd_.get_measurements()[m_ip_].c_reg_] = mres;
                            probs_[qcd_.get_measurements()[m_ip_].c_reg_] =
                                probs[mres];
                            for (auto&& i :
                                 qcd_.get_measurements()[m_ip_].target_)
                                set_measured_(i);
                            break;
                        }
                        ++ip_;
                        if (++m_ip_ == qcd_.get_measurements().size())
                            break;
                    }
                }
            }

            // check for gates
            if (q_ip_ < qcd_.get_gates().size()) {
                if (ip_ - saved_ip_ == no_steps)
                    return;
                if (verbose) {
                    std::cout << std::left;
                    std::cout << std::setw(8)
                              << qcd_.get_gates()[q_ip_].step_no_;
                    std::cout << std::right;
                    std::cout << qcd_.get_gates()[q_ip_] << '\n';
                }

                std::vector<idx> ctrl_rel_pos;
                std::vector<idx> target_rel_pos =
                    get_relative_pos_(qcd_.get_gates()[q_ip_].target_);

                switch (qcd_.get_gates()[q_ip_].gate_type_) {
                case QCircuitDescription::GateType::NONE:
                    break;
                case QCircuitDescription::GateType::SINGLE:
                case QCircuitDescription::GateType::TWO:
                case QCircuitDescription::GateType::THREE:
                case QCircuitDescription::GateType::CUSTOM:
                    psi_ = apply(psi_, qcd_.get_gates()[q_ip_].gate_,
                                 target_rel_pos, qcd_.get_d());
                    break;
                case QCircuitDescription::GateType::FAN:
                    for (idx m = 0; m < qcd_.get_gates()[q_ip_].target_.size();
                         ++m)
                        psi_ = apply(psi_, qcd_.get_gates()[q_ip_].gate_,
                                     {target_rel_pos[m]}, qcd_.get_d());
                    break;
                case QCircuitDescription::GateType::QFT:
                case QCircuitDescription::GateType::TFQ:
                case QCircuitDescription::GateType::SINGLE_CTRL_SINGLE_TARGET:
                case QCircuitDescription::GateType::SINGLE_CTRL_MULTIPLE_TARGET:
                case QCircuitDescription::GateType::MULTIPLE_CTRL_SINGLE_TARGET:
                case QCircuitDescription::GateType::
                    MULTIPLE_CTRL_MULTIPLE_TARGET:
                case QCircuitDescription::GateType::CUSTOM_CTRL:
                    ctrl_rel_pos =
                        get_relative_pos_(qcd_.get_gates()[q_ip_].ctrl_);
                    psi_ =
                        applyCTRL(psi_, qcd_.get_gates()[q_ip_].gate_,
                                  ctrl_rel_pos, target_rel_pos, qcd_.get_d());
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
                        psi_ = apply(psi_, qcd_.get_gates()[q_ip_].gate_,
                                     target_rel_pos, qcd_.get_d());
                    } else {
                        bool should_apply = true;
                        idx first_dit =
                            dits_[(qcd_.get_gates()[q_ip_].ctrl_)[0]];
                        for (idx m = 0;
                             m < qcd_.get_gates()[q_ip_].ctrl_.size(); ++m) {
                            if (dits_[(qcd_.get_gates()[q_ip_].ctrl_)[m]] !=
                                first_dit) {
                                should_apply = false;
                                break;
                            }
                        }
                        if (should_apply) {
                            psi_ =
                                apply(psi_, powm(qcd_.get_gates()[q_ip_].gate_,
                                                 first_dit),
                                      target_rel_pos, qcd_.get_d());
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
