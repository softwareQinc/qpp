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
 * \file classes/circuits.h
 * \brief Support for qudit quantum circuits
 */

#ifndef CLASSES_CIRCUITS_H_
#define CLASSES_CIRCUITS_H_

namespace qpp {
/**
 * \class qpp::QCircuitDescription
 * \brief Quantum circuit description class
 * \see qpp::QCircuit
 */
class QCircuitDescription : public IDisplay {
    const idx nq_;                         ///< number of qudits
    const idx nc_;                         ///< number of classical "dits"
    const idx d_;                          ///< qudit dimension
    std::vector<idx> measurement_steps_{}; ///< keeps track of where the
    ///< measurements take place
    std::string name_;           ///< optional circuit name
    std::vector<bool> measured_; ///< keeps track of the measured qudits
    idx steps_cnt_;              ///< step counter
  public:
    /**
     * \brief Type of gate being executed in a gate step
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
        ///< multiple controls and multiple targets

        CUSTOM_CTRL, ///< custom controlled gate with multiple controls
        ///< and multiple targets

        SINGLE_cCTRL_SINGLE_TARGET, ///< controlled 1 qudit unitary gate with
        ///< one classical control and one target

        SINGLE_cCTRL_MULTIPLE_TARGET, ///< controlled 1 qudit unitary gate with
        ///< one classical control and multiple targets

        MULTIPLE_cCTRL_SINGLE_TARGET, ///< controlled 1 qudit unitary gate with
        ///< multiple classical controls and single target

        MULTIPLE_cCTRL_MULTIPLE_TARGET, ///< controlled 1 qudit unitary gate
        ///< with multiple classical controls and multiple targets

        CUSTOM_cCTRL, ///< custom controlled gate with multiple controls and
        ///< multiple targets
    };
    /**
     * \brief Type of measurement being executed in a measurement step
     */
    enum class MeasureType {
        NONE, ///< represents no measurement

        MEASURE_Z, ///< Z measurement of single qudit

        MEASURE_V, ///< measurement of single qudit in the orthonormal basis
        ///< or rank-1 projectors specified by the columns of matrix \a V

        MEASURE_V_MANY, ///< measurement of multiple qudits in the orthonormal
        ///< basis or rank-1 projectors specified by the columns of matrix \a V
    };

  public:
    /**
     * \class qpp::QCircuitDescription::iterator
     * \brief Quantum circuit description bound-checking (safe) iterator
     * \note The iterator is a const_iterator by default
     */
    class iterator {
        friend QCircuitDescription;
        friend class IQCircuit;

        ///< non-owning pointer to const circuit description
        const QCircuitDescription* qcd_{nullptr};

        ///< iterator value type
        struct value_type_ : public IDisplay {
            bool is_measurement_{false}; ///< current step is a measurement
            idx m_ip_{idx_infty};        ///< measurements instruction pointer
            idx q_ip_{idx_infty};        ///< gates instruction pointer
            idx ip_{idx_infty}; ///< total (measurements + gates) instruction
            ///< pointer

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

            ///< non-owning pointer to the parent iterator
            const QCircuitDescription* value_type_qcd_;
            explicit value_type_(const QCircuitDescription* value_type_qcd)
                : value_type_qcd_{value_type_qcd} {}

            /**
             * \brief qpp::IDisplay::display() override
             *
             * Writes to the output stream the textual representation of the
             * iterator de-referenced element
             *
             * \param os Output stream passed by reference
             * \return Reference to the output stream
             */
            std::ostream& display(std::ostream& os) const override {
                // os << "m_ip_: " << m_ip_;
                // os << ", q_ip_: " << q_ip_;
                // os << ", ip_: " << ip_;

                // field spacing for the step number
                idx text_width =
                    std::to_string(value_type_qcd_->get_steps_count()).size() +
                    1;

                if (is_measurement_) {
                    os << std::left;
                    os << std::setw(text_width) << ip_;
                    os << std::right;
                    os << "|> " << value_type_qcd_->get_measurements()[m_ip_];
                } else {
                    os << std::left;
                    os << std::setw(text_width) << ip_;
                    os << std::right;
                    os << value_type_qcd_->get_gates()[q_ip_];
                }

                return os;
            }
        };

        value_type_ elem_{nullptr}; ///< de-referenced iterator element

        /**
         * \brief Sets the internal quantum circuit description pointer
         *
         * \param qcd Constant pointer to a quantum circuit description
         */
        void set_(const QCircuitDescription* qcd) {
            qcd_ = qcd;
            elem_ = value_type_{qcd};

            // if the circuit is empty, then all m_ip_, q_ip_ and ip_ are
            // equal to idx_infty (i.e. idx(-1))

            // if the circuit has measurements
            if (qcd->get_measurement_count() != 0) {
                elem_.m_ip_ = 0;
                elem_.ip_ = 0;
                if (qcd->get_measurement_steps()[0] == 0) {
                    elem_.is_measurement_ = true;
                }
            }

            // if the circuit has gates
            if (qcd->get_gate_count() != 0) {
                elem_.q_ip_ = 0;
                elem_.ip_ = 0;
            }
        }

      public:
        /**
         * \brief Default constructor
         */
        iterator() = default;

        // silence -Weffc++ class has pointer data members
        /**
         * \brief Default copy constructor
         */
        iterator(const iterator&) = default;

        // silence -Weffc++ class has pointer data members
        /**
         * \brief Default copy assignment operator
         *
         * \return Reference to the current instance
         */
        iterator& operator=(const iterator&) = default;

        /**
         * \brief Prefix increment operator
         *
         * \return Reference to the current instance
         */
        iterator& operator++() {
            // EXCEPTION CHECKS

            // protects against incrementing invalid iterators
            if (qcd_ == nullptr) {
                throw qpp::exception::InvalidIterator(
                    "qpp::QCircuitDescription::iterator::operator++()");
            }

            // protects against incrementing past the end
            if (elem_.ip_ == qcd_->get_steps_count()) {
                throw qpp::exception::InvalidIterator(
                    "qpp::QCircuitDescription::iterator::operator++()");
            }

            // protect against incrementing an empty circuit iterator
            if (elem_.q_ip_ == idx_infty && elem_.m_ip_ == idx_infty) {
                throw qpp::exception::InvalidIterator(
                    "qpp::QCircuitDescription::iterator::operator++()");
            }
            // END EXCEPTION CHECKS

            // increment the instruction pointer
            ++elem_.ip_;

            // only measurements, no gates
            if (elem_.q_ip_ == idx_infty) {
                elem_.is_measurement_ = true;
                ++elem_.m_ip_;
                return *this;
            }

            // only gates, no measurements
            if (elem_.m_ip_ == idx_infty) {
                elem_.is_measurement_ = false;
                ++elem_.q_ip_;
                return *this;
            }

            // current step is a measurement
            if (elem_.m_ip_ < qcd_->get_measurement_count() &&
                qcd_->measurement_steps_[elem_.m_ip_] == elem_.q_ip_) {
                // next step is a measurement
                if (elem_.m_ip_ + 1 < qcd_->get_measurement_count() &&
                    qcd_->measurement_steps_[elem_.m_ip_ + 1] == elem_.q_ip_) {
                    elem_.is_measurement_ = true;
                    ++elem_.m_ip_;
                } else
                // next step is a gate
                {
                    elem_.is_measurement_ = false;
                    ++elem_.m_ip_;
                }
            } else
            // current step is a gate
            {
                // next step is a measurement
                if (elem_.m_ip_ < qcd_->get_measurement_count() &&
                    qcd_->measurement_steps_[elem_.m_ip_] == elem_.q_ip_ + 1) {
                    elem_.is_measurement_ = true;
                    ++elem_.q_ip_;
                } else
                // next step is a gate
                {
                    elem_.is_measurement_ = false;
                    ++elem_.q_ip_;
                }
            }

            return *this;
        }

        /**
         * \brief Postfix increment operator
         *
         * \return Copy of the current instance before the increment
         */
        iterator operator++(int) {
            iterator retval = *this;
            ++(*this);
            return retval;
        }

        /**
         * \brief Equality operator
         *
         * \param rhs Iterator against which the equality is being tested
         * \return True if the iterators are equal, false otherwise
         */
        bool operator==(const iterator& rhs) const {
            return elem_.ip_ == rhs.elem_.ip_;
        }

        /**
         * \brief Inequality operator
         *
         * \param rhs Iterator against which the inequality is being tested
         * \return True if the iterators are not equal (bit by bit), false
         * otherwise
         */
        bool operator!=(iterator rhs) const { return !(*this == rhs); }

        /**
         * \brief Safe de-referencing operator
         *
         * \return Constant reference to the iterator element
         */
        const value_type_& operator*() const {
            // EXCEPTION CHECKS

            // protect against de-referencing past the last element or against
            // de-referencing invalid iterators
            if (qcd_ == nullptr || elem_.ip_ == qcd_->get_steps_count())
                throw exception::InvalidIterator(
                    "qpp::QCircuitDescription::iterator::operator*()");
            // END EXCEPTION CHECKS

            return elem_;
        }

        // iterator traits
        using difference_type = long long;                   ///< iterator trait
        using value_type = value_type_;                      ///< iterator trait
        using pointer = const value_type*;                   ///< iterator trait
        using reference = const value_type&;                 ///< iterator trait
        using iterator_category = std::forward_iterator_tag; ///< iterator trait
    };

    using const_iterator = iterator; ///< both iterators are const_iterators

    /**
     * \brief Iterator to the first element
     *
     * \return Iterator to the first element
     */
    iterator begin() {
        iterator it;
        it.set_(this);
        return it;
    }

    /**
     * \brief Constant iterator to the first element
     *
     * \return Constant iterator to the first element
     */
    const_iterator begin() const noexcept {
        iterator it;
        it.set_(this);
        return it;
    }

    /**
     * \brief Constant iterator to the first element
     *
     * \return Constant iterator to the first element
     */
    const_iterator cbegin() const noexcept {
        iterator it;
        it.set_(this);
        return it;
    }

    /**
     * \brief Iterator to the next to the last element
     *
     * \return Iterator to the next to the last element
     */
    iterator end() {
        iterator it;
        it.set_(this);
        if (steps_cnt_ == 0) // empty circuit
            return it;
        else
            it.elem_.ip_ = steps_cnt_;
        return it;
    }

    /**
     * \brief Constant iterator to the next to the last element
     *
     * \return Constant iterator to the next to the last element
     */
    const_iterator end() const noexcept {
        iterator it;
        it.set_(this);
        if (steps_cnt_ == 0) // empty circuit
            return it;
        else
            it.elem_.ip_ = steps_cnt_;
        return it;
    }

    /**
     * \brief Constant iterator to the next to the last element
     *
     * \return Constant iterator to the next to the last element
     */
    const_iterator cend() const noexcept {
        iterator it;
        it.set_(this);
        if (steps_cnt_ == 0) // empty circuit
            return it;
        else
            it.elem_.ip_ = steps_cnt_;
        return it;
    }

  private:
    /**
     * \brief One step consisting only of gates/operators in the circuit
     */
    struct GateStep {
        GateType gate_type_ = GateType::NONE; ///< gate type
        cmat gate_;                           ///< gate
        std::vector<idx> ctrl_;               ///< control
        std::vector<idx> target_; ///< target where the gate is applied
        idx step_no_;             ///< step number
        std::string name_;        ///< custom name of the step
        /**
         * \brief Default constructor
         */
        GateStep() = default;
        /**
         * \brief Constructs a gate step instance
         *
         * \param gate_type Gate type
         * \param gate Quantum gate
         * \param ctrl Control qudit indexes
         * \param target Target qudit indexes
         * \param step_no Circuit step number
         * \param name Optional gate name
         */
        explicit GateStep(GateType gate_type, const cmat& gate,
                          const std::vector<idx>& ctrl,
                          const std::vector<idx>& target, idx step_no,
                          std::string name = "")
            : gate_type_{gate_type}, gate_{gate}, ctrl_{ctrl}, target_{target},
              step_no_{step_no}, name_{name} {}
    };

    /**
     * \brief One step consisting only of measurements in the circuit
     */
    struct MeasureStep {
        MeasureType measurement_type_ = MeasureType::NONE; ///< measurement type
        std::vector<cmat> mats_; ///< matrix/matrices that specify the
        /// measurement
        std::vector<idx> target_; ///< target where the measurement is applied
        idx c_reg_{}; ///< index of the classical register where the measurement
        ///< result is being stored
        idx step_no_;      ///< step number
        std::string name_; ///< custom name of the step
        /**
         * \brief Default constructor
         */
        MeasureStep() = default;
        /**
         * \brief Constructs a measurement step instance
         *
         * \param measurement_type Measurement type
         * \param mats Vector of measurement matrices (can be only one or many
         * for Kraus measurements)
         * \param target Target qudit indexes
         * \param c_reg Classical register where the value of the measurement is
         * stored
         * \param step_no Circuit step number
         * \param name Optional gate name
         */
        explicit MeasureStep(MeasureType measurement_type,
                             const std::vector<cmat>& mats,
                             const std::vector<idx>& target, idx c_reg,
                             idx step_no, std::string name = "")
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
            os << "ctrl = " << disp(gate_step.ctrl_, ", ") << ", ";
        os << "target = " << disp(gate_step.target_, ", ") << ", ";
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
        os << "target = " << disp(measure_step.target_, ", ") << ", ";
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
    explicit QCircuitDescription(idx nq, idx nc = 0, idx d = 2,
                                 std::string name = "")
        : nq_{nq}, nc_{nc}, d_{d}, name_{name},
          measured_(nq, false), steps_cnt_{0} {
        // EXCEPTION CHECKS

        if (nq == 0)
            throw exception::ZeroSize(
                "qpp::QCircuitDescription::QCircuitDescription()");
        if (d < 2)
            throw exception::OutOfRange(
                "qpp::QCircuitDescription::QCircuitDescription()");
        // END EXCEPTION CHECKS
    }

    /**
     * \brief Default virtual destructor
     */
    virtual ~QCircuitDescription() = default;

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
     * \brief Quantum circuit total steps count, i.e. the sum of gate count and
     * measurement count
     *
     * \return Total (gates + measurements) count
     */
    idx get_steps_count() const noexcept { return steps_cnt_; }
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
        } catch (qpp::exception::Exception&) {
            std::cerr << "At STEP " << steps_cnt_ << "\n";
            throw;
        }
        // END EXCEPTION CHECKS

        if (name == "")
            name = qpp::Gates::get_instance().get_name(U);
        gates_.emplace_back(GateType::SINGLE, U, std::vector<idx>{},
                            std::vector<idx>{i}, steps_cnt_++, name);

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
        } catch (qpp::exception::Exception&) {
            std::cerr << "At STEP " << steps_cnt_ << "\n";
            throw;
        }
        // END EXCEPTION CHECKS

        if (name == "")
            name = qpp::Gates::get_instance().get_name(U);
        gates_.emplace_back(GateType::TWO, U, std::vector<idx>{},
                            std::vector<idx>{i, j}, steps_cnt_++, name);

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
        } catch (qpp::exception::Exception&) {
            std::cerr << "At STEP " << steps_cnt_ << "\n";
            throw;
        }
        // END EXCEPTION CHECKS

        if (name == "")
            name = qpp::Gates::get_instance().get_name(U);
        gates_.emplace_back(GateType::THREE, U, std::vector<idx>{},
                            std::vector<idx>{i, j, k}, steps_cnt_++, name);

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
        } catch (qpp::exception::Exception&) {
            std::cerr << "At STEP " << steps_cnt_ << "\n";
            throw;
        }
        // END EXCEPTION CHECKS

        if (name == "")
            name = qpp::Gates::get_instance().get_name(U);
        gates_.emplace_back(GateType::FAN, U, std::vector<idx>{}, target,
                            steps_cnt_++, name);

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
        } catch (qpp::exception::Exception&) {
            std::cerr << "At STEP " << steps_cnt_ << "\n";
            throw;
        }
        // END EXCEPTION CHECKS

        if (name == "")
            name = qpp::Gates::get_instance().get_name(U);
        gates_.emplace_back(GateType::FAN, U, std::vector<idx>{},
                            get_non_measured(), steps_cnt_++, name);

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
        } catch (qpp::exception::Exception&) {
            std::cerr << "At STEP " << steps_cnt_ << "\n";
            throw;
        }
        // END EXCEPTION CHECKS

        if (name == "")
            name = qpp::Gates::get_instance().get_name(U);
        gates_.emplace_back(GateType::CUSTOM, U, std::vector<idx>{}, target,
                            steps_cnt_++, name);

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
        } catch (qpp::exception::Exception&) {
            std::cerr << "At STEP " << steps_cnt_ << "\n";
            throw;
        }
        // END EXCEPTION CHECKS

        gates_.emplace_back(GateType::QFT, cmat{}, std::vector<idx>{}, target,
                            steps_cnt_++, "QFT");

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
        } catch (qpp::exception::Exception&) {
            std::cerr << "At STEP " << steps_cnt_ << "\n";
            throw;
        }
        // END EXCEPTION CHECKS
        gates_.emplace_back(GateType::TFQ, cmat{}, std::vector<idx>{}, target,
                            steps_cnt_++, "TFQ");

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
        } catch (qpp::exception::Exception&) {
            std::cerr << "At STEP " << steps_cnt_ << "\n";
            throw;
        }
        // END EXCEPTION CHECKS

        if (name == "")
            name = qpp::Gates::get_instance().get_name(U);
        gates_.emplace_back(GateType::SINGLE_CTRL_SINGLE_TARGET, U,
                            std::vector<idx>{ctrl}, std::vector<idx>{target},
                            steps_cnt_++, name);

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
        } catch (qpp::exception::Exception&) {
            std::cerr << "At STEP " << steps_cnt_ << "\n";
            throw;
        }
        // END EXCEPTION CHECKS

        if (name == "")
            name = qpp::Gates::get_instance().get_name(U);
        gates_.emplace_back(GateType::SINGLE_CTRL_MULTIPLE_TARGET, U,
                            std::vector<idx>{ctrl}, target, steps_cnt_++, name);

        return *this;
    }

    // multiple ctrl single target
    /**
     * \brief Applies the single qudit controlled gate \a U with multiple
     * control qudits listed in \a ctrl on the target qudit \a target
     *
     * \param U Single qudit quantum gate
     * \param ctrl Control qudit indexes
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
        } catch (qpp::exception::Exception&) {
            std::cerr << "At STEP " << steps_cnt_ << "\n";
            throw;
        }
        // END EXCEPTION CHECKS

        if (name == "")
            name = qpp::Gates::get_instance().get_name(U);
        gates_.emplace_back(GateType::MULTIPLE_CTRL_SINGLE_TARGET, U, ctrl,
                            std::vector<idx>{target}, steps_cnt_++, name);

        return *this;
    }

    // multiple ctrl multiple target
    // FIXME
    /**
     * \brief Applies the single qudit controlled gate \a U with multiple
     * control qudits listed in \a ctrl on every qudit listed in \a target
     *
     * \param U Single qudit quantum gate
     * \param ctrl Control qudit indexes
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
        } catch (qpp::exception::Exception&) {
            std::cerr << "At STEP " << steps_cnt_ << "\n";
            throw;
        }
        // END EXCEPTION CHECKS

        if (name == "")
            name = qpp::Gates::get_instance().get_name(U);
        gates_.emplace_back(GateType::MULTIPLE_CTRL_MULTIPLE_TARGET, U, ctrl,
                            std::vector<idx>{target}, steps_cnt_++, name);

        return *this;
    }

    // custom multiple control composed target
    /**
     * \brief Jointly applies the custom multiple-qudit controlled gate \a U
     * with multiple control qudits listed in \a ctrl on the qudit indexes
     * specified by \a target
     *
     * \param U Multiple-qudit quantum gate
     * \param ctrl Control qudit indexes
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
        } catch (qpp::exception::Exception&) {
            std::cerr << "At STEP " << steps_cnt_ << "\n";
            throw;
        }
        // END EXCEPTION CHECKS

        if (name == "")
            name = qpp::Gates::get_instance().get_name(U);
        gates_.emplace_back(GateType::CUSTOM_CTRL, U, ctrl, target,
                            steps_cnt_++, name);

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
        } catch (qpp::exception::Exception&) {
            std::cerr << "At STEP " << steps_cnt_ << "\n";
            throw;
        }
        // END EXCEPTION CHECKS

        if (name == "")
            name = qpp::Gates::get_instance().get_name(U);
        gates_.emplace_back(GateType::SINGLE_cCTRL_SINGLE_TARGET, U,
                            std::vector<idx>{ctrl_dit},
                            std::vector<idx>{target}, steps_cnt_++, name);

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
        } catch (qpp::exception::Exception&) {
            std::cerr << "At STEP " << steps_cnt_ << "\n";
            throw;
        }
        // END EXCEPTION CHECKS

        if (name == "")
            name = qpp::Gates::get_instance().get_name(U);
        gates_.emplace_back(GateType::SINGLE_cCTRL_MULTIPLE_TARGET, U,
                            std::vector<idx>{ctrl_dit}, target, steps_cnt_++,
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
        } catch (qpp::exception::Exception&) {
            std::cerr << "At STEP " << steps_cnt_ << "\n";
            throw;
        }
        // END EXCEPTION CHECKS

        if (name == "")
            name = qpp::Gates::get_instance().get_name(U);
        gates_.emplace_back(GateType::MULTIPLE_cCTRL_SINGLE_TARGET, U,
                            ctrl_dits, std::vector<idx>{target}, steps_cnt_++,
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
        } catch (qpp::exception::Exception&) {
            std::cerr << "At STEP " << steps_cnt_ << "\n";
            throw;
        }
        // END EXCEPTION CHECKS

        if (name == "")
            name = qpp::Gates::get_instance().get_name(U);
        gates_.emplace_back(GateType::MULTIPLE_cCTRL_MULTIPLE_TARGET, U,
                            ctrl_dits, std::vector<idx>{target}, steps_cnt_++,
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
        } catch (qpp::exception::Exception&) {
            std::cerr << "At STEP " << steps_cnt_ << "\n";
            throw;
        }
        // END EXCEPTION CHECKS

        if (name == "")
            name = qpp::Gates::get_instance().get_name(U);
        gates_.emplace_back(GateType::CUSTOM_cCTRL, U, ctrl_dits, target,
                            steps_cnt_++, name);

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
        } catch (qpp::exception::Exception&) {
            std::cerr << "At STEP " << steps_cnt_ << "\n";
            throw;
        }
        // END EXCEPTION CHECKS

        if (name == "")
            name = "Measure Z";
        measured_[i] = true;
        measurements_.emplace_back(MeasureType::MEASURE_Z, std::vector<cmat>{},
                                   std::vector<idx>{i}, c_reg, steps_cnt_++,
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
        } catch (qpp::exception::Exception&) {
            std::cerr << "At STEP " << steps_cnt_ << "\n";
            throw;
        }
        // END EXCEPTION CHECKS

        if (name == "")
            name = qpp::Gates::get_instance().get_name(V);
        measured_[i] = true;
        measurements_.emplace_back(MeasureType::MEASURE_V, std::vector<cmat>{V},
                                   std::vector<idx>{i}, c_reg, steps_cnt_++,
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
        } catch (qpp::exception::Exception&) {
            std::cerr << "At STEP " << steps_cnt_ << "\n";
            throw;
        }
        // END EXCEPTION CHECKS

        if (name == "")
            name = qpp::Gates::get_instance().get_name(V);
        for (auto&& i : target)
            measured_[i] = true;
        measurements_.emplace_back(MeasureType::MEASURE_V_MANY,
                                   std::vector<cmat>{V}, target, c_reg,
                                   steps_cnt_++, name);
        measurement_steps_.emplace_back(gates_.size());

        return *this;
    }

    /**
     * \brief qpp::IDisplay::display() override
     *
     * Writes to the output stream a textual representation of the quantum
     * circuit description
     *
     * \param os Output stream passed by reference
     * \return Reference to the output stream
     */
    std::ostream& display(std::ostream& os) const override {
        os << "nq = " << nq_ << ", nc = " << nc_ << ", d = " << d_;

        if (name_ != "") // if the circuit is named
            os << ", name = \"" << name_ << "\"\n";
        else
            os << ", name = \"\"\n";

        for (auto&& elem : *this) {
            os << elem << '\n';
        }

        os << "measurement steps: " << disp(get_measurement_steps(), ", ")
           << '\n';
        os << "measured positions: " << disp(get_measured(), ", ") << '\n';
        os << "non-measured positions: " << disp(get_non_measured(), ", ");

        return os;
    }

    /**
     * \brief Displays the circuit description in JSON format
     *
     * \param enclosed_in_curly_brackets Encloses the result in curly brackets
     * if true
     * \return String containing the JSON representation of the circuit
     * description
     */
    std::string to_JSON(bool enclosed_in_curly_brackets = true) const {
        std::string result;

        if (enclosed_in_curly_brackets)
            result += "{";

        result += "\"nq\" : " + std::to_string(nq_);
        result += ", \"nc\" : " + std::to_string(nc_);
        result += ", \"d\" : " + std::to_string(d_);
        result += ", \"name\" : \"" + name_ + "\"";

        bool is_first = true;
        std::ostringstream ss;
        result += ", \"steps\" : [";
        for (auto&& elem : *this) {
            if (is_first) {
                is_first = false;
            } else {
                result += ", ";
            }
            result += "{\"step\" : " + std::to_string(elem.ip_) + ", ";
            result += "\"type\" : ";
            if (elem.is_measurement_) {
                ss.str("");
                ss.clear();
                ss << measurements_[elem.m_ip_].measurement_type_;
                result += "\"" + ss.str() + "\", ";
                ss.str("");
                ss.clear();
                ss << disp(measurements_[elem.m_ip_].target_, ", ");
                result += "\"target\" : " + ss.str() + ", ";
                result += "\"c_reg\" : " +
                          std::to_string(measurements_[elem.m_ip_].c_reg_) +
                          ", ";
                result += "\"name\" : ";
                result += "\"" + measurements_[elem.m_ip_].name_ + "\"" + "}";
            } else {
                ss.str("");
                ss.clear();
                ss << gates_[elem.q_ip_].gate_type_;
                result += "\"" + ss.str() + "\", ";
                if (gates_[elem.q_ip_].ctrl_.size() != 0) {
                    ss.str("");
                    ss.clear();
                    ss << disp(gates_[elem.q_ip_].ctrl_, ", ");
                    result += "\"ctrl\" : " + ss.str() + ", ";
                }
                ss.str("");
                ss.clear();
                ss << disp(gates_[elem.q_ip_].target_, ", ");
                result += "\"target\" : " + ss.str() + ", ";
                result += "\"name\" : ";
                result += "\"" + gates_[elem.q_ip_].name_ + "\"" + "}";
            }
        }              // end for
        result += "]"; // end steps

        ss.str("");
        ss.clear();
        ss << disp(get_measurement_steps(), ", ");
        result += ", \"measurement steps\" : " + ss.str() + ", ";

        ss.str("");
        ss.clear();
        ss << disp(get_measured(), ", ");
        result += "\"measured positions\" : " + ss.str() + ", ";

        ss.str("");
        ss.clear();
        ss << disp(get_non_measured(), ", ");
        result += "\"non-measured positions\" : " + ss.str();

        if (enclosed_in_curly_brackets)
            result += "}";

        return result;
    }
}; /* class QCircuitDescription */

/**
 * \class qpp::IQCircuit
 * \brief Quantum circuit simulator abstract class
 * \see qpp::QCircuitDescription
 * \note Every further derived class has to override the run() member
 * function
 */
class IQCircuit : public IDisplay {
  protected:
    const QCircuitDescription& qcd_; ///< quantum circuit description
    ket psi_;                        ///< state vector
    std::vector<idx> dits_;          ///< classical dits
    std::vector<double> probs_;      ///< measurement probabilities
    std::vector<idx> subsys_; ///< keeps track of the measured subsystems,
    ///< relabel them after measurements

    QCircuitDescription::const_iterator it_; ///< iterator to current step

    /**
     * \brief Marks qudit \a i as measured then re-label accordingly the
     * remaining non-measured qudits
     * \param i Qudit index
     */
    void set_measured_(idx i) {
        if (get_measured(i))
            throw exception::QuditAlreadyMeasured(
                "qpp::QCircuit::set_measured_()");
        subsys_[i] = idx_infty; // set qudit i to measured state
        for (idx m = i; m < qcd_.get_nq(); ++m) {
            if (!get_measured(m)) {
                --subsys_[m];
            }
        }
    }

    // giving a vector of non-measured qudits, get their relative
    // position wrt the measured qudits
    /**
     * \brief Giving a vector \a V of non-measured qudits, get their
     * relative position with respect to the measured qudits \param v
     * Qudit index
     */
    std::vector<idx> get_relative_pos_(std::vector<idx> v) {
        idx vsize = v.size();
        for (idx i = 0; i < vsize; ++i) {
            if (get_measured(v[i]))
                throw exception::QuditAlreadyMeasured(
                    "qpp::QCircuit::get_relative_pos_()");
            v[i] = subsys_[v[i]];
        }
        return v;
    }

  public:
    /**
     * \brief Constructs a quantum circuit out of a quantum circuit
     * description
     *
     * \note The quantum circuit description must be an lvalue
     * \see qpp::QCircuit(QCircuitDescription&&)
     *
     * \note The initial underlying quantum state is set to
     * \f$|0\rangle^{\otimes n}\f$
     *
     * \param qcd Quantum circuit description
     */
    explicit IQCircuit(const QCircuitDescription& qcd)
        : qcd_{qcd}, psi_{States::get_instance().zero(qcd.get_nq(),
                                                      qcd.get_d())},
          dits_(qcd.get_nc(), 0), probs_(qcd.get_nc(), 0),
          subsys_(qcd.get_nq(), 0), it_{qcd_.begin()} {
        std::iota(std::begin(subsys_), std::end(subsys_), 0);
    }

    /**
     * \brief Disables rvalue QCircuitDescription
     */
    IQCircuit(QCircuitDescription&&) = delete;

    /**
     * \brief Default virtual destructor
     */
    virtual ~IQCircuit() = default;

    // getters
    /**
     * \brief Underlying quantum state
     *
     * \return Underlying quantum state
     */
    ket get_psi() const { return psi_; }

    /**
     * \brief Reference to the underlying quantum state
     *
     * \return Reference to the underlying quantum state
     */
    ket& get_ref_psi() { return psi_; }

    /**
     * \brief Vector with the values of the underlying classical dits
     *
     * \return Vector of underlying classical dits
     */
    std::vector<idx> get_dits() const { return dits_; }

    /**
     * \brief Value of the classical dit at position \a i
     *
     * \param i Classical dit index
     *
     * \return Value of the classical dit at position \a i
     */
    idx get_dit(idx i) const {
        if (i > qcd_.get_nc())
            throw exception::OutOfRange("qpp::QCircuit::get_dit()");

        return dits_[i];
    }

    /**
     * \brief Vector of underlying measurement outcome probabilities
     *
     * \note The probability vector has the same length as the vector of
     * classical dits. If the measurement result is stored at the index
     * \a c_reg, then the outcome probability is automatically stored at
     * the same index \a c_reg in the probability vector.
     *
     * \return Vector of underlying measurement outcome probabilities
     */
    std::vector<double> get_probs() const { return probs_; }

    /**
     * \brief Check whether qudit \a i was already measured
     *
     * \param i Qudit index
     * \return True if qudit \a i was already measured, false othwewise
     */
    idx get_measured(idx i) const { return subsys_[i] == idx_infty; }

    /**
     * \brief Vector of already measured qudit indexes
     *
     * \return Vector of already measured qudit indexes
     */
    std::vector<idx> get_measured() const {
        std::vector<idx> result;
        for (idx i = 0; i < qcd_.get_nq(); ++i)
            if (get_measured(i))
                result.emplace_back(i);

        return result;
    }

    /**
     * \brief Vector of non-measured qudit indexes
     *
     * \return Vector of non-measured qudit indexes
     */
    std::vector<idx> get_not_measured() const {
        std::vector<idx> result;
        for (idx i = 0; i < qcd_.get_nq(); ++i)
            if (!get_measured(i))
                result.emplace_back(i);

        return result;
    }

    /**
     * \brief Checks whether the current step in the circuit is a
     * measurement step
     *
     * \return True if measurement step, false otherwise
     */
    bool is_measurement_step() const { return it_.elem_.is_measurement_; }

    /**
     * \brief Measurement instruction pointer
     *
     * Points to the index of the next measurement to be executed from
     * the std::vector<MeasureStep> of measurements in the circuit
     * description
     *
     * \return Measurement instruction pointer
     */
    idx get_m_ip() const { return it_.elem_.m_ip_; }

    /**
     * \brief Quantum instruction pointer
     *
     * Points to the index of the next quantum gate to be executed from
     * the std::vector<GateStep> of quantum gates in the circuit
     * description
     *
     * \return Quantum instruction pointer
     */
    idx get_q_ip() const { return it_.elem_.q_ip_; }

    /**
     * \brief Total instruction pointer
     *
     * \return The sum of measurement instruction pointer and quantum
     * instruction pointer
     */
    idx get_ip() const { return it_.elem_.ip_; }

    /**
     * \brief Iterator to current step
     *
     * \return Iterator to current step in the circuit
     */
    QCircuitDescription::const_iterator get_iter() const { return it_; }

    /**
     * \brief Quantum circuit description
     *
     * \return Quantum circuit description
     */
    const QCircuitDescription& get_circuit_description() const { return qcd_; }
    // end getters

    // setters
    /**
     * \brief Sets the classical dit at position \a i
     *
     * \param i Classical dit index
     * \param value Classical dit value
     * \return Reference to the current instance
     */
    IQCircuit& set_dit(idx i, idx value) {
        if (i > qcd_.get_nc())
            throw exception::OutOfRange("qpp::QCircuit::set_dit()");
        dits_[i] = value;

        return *this;
    }
    // end setters

    /**
     * \brief Resets the quantum circuit
     *
     * Re-initializes everything to zero and sets the initial state to
     * \f$|0\rangle^{\otimes n}\f$
     */
    void reset() {
        psi_ = States::get_instance().zero(qcd_.get_nq(), qcd_.get_d());
        dits_ = std::vector<idx>(qcd_.get_nc(), 0);
        probs_ = std::vector<double>(qcd_.get_nc(), 0);
        std::iota(std::begin(subsys_), std::end(subsys_), 0);
        it_ = qcd_.begin();
    }

    /**
     * \brief qpp::IDisplay::display() override
     *
     * Writes to the output stream a textual representation of the
     * quantum circuit
     *
     * \param os Output stream passed by reference
     * \return Reference to the output stream
     */
    std::ostream& display(std::ostream& os) const override {
        os << "measured: " << disp(get_measured(), ", ") << '\n';
        os << "dits: " << disp(get_dits(), ", ") << '\n';
        os << "probs: " << disp(get_probs(), ", ");

        return os;
    }

    /**
     * \brief Executes the quantum circuit, pure virtual member function
     * that must be overridden by all derived classes
     *
     * \param verbose If true, displays at console every executed step
     * \param step How many steps to execute, by default executes until
     * the end
     */
    virtual void run(bool verbose = false, idx step = idx_infty) = 0;
}; /* class IQCircuit */

/**
 * \class qpp::QCircuit
 * \brief Quantum circuit simulator class
 * \see qpp::QCircuitDescription
 */
class QCircuit : public IQCircuit {
  public:
    using IQCircuit::IQCircuit; ///< Uses the base IQCircuit constructor

    /**
     * \brief Executes the quantum circuit, qpp::IQCircuit::run() override
     *
     * \param verbose If true, displays at console every executed step
     * \param step How many steps to execute, by default executes until the end
     */
    void run(bool verbose = false, idx step = idx_infty) override {
        // EXCEPTION CHECKS

        // trying to run an empty circuit
        if (qcd_.get_steps_count() == 0) // same as if(get_ip() == idx_infty)
            throw exception::ZeroSize("qpp::QCircuit::run()");

        idx no_steps;
        if (step == idx_infty)
            no_steps = qcd_.get_steps_count() - get_ip();
        else
            no_steps = step;

        if (get_ip() + no_steps > qcd_.get_steps_count())
            throw exception::OutOfRange("qpp::QCircuit::run()");
        // END EXCEPTION CHECKS

        // no steps to run
        if (step == 0)
            return;

        // main loop
        for (idx i = 0; i < no_steps; ++i) {

            if (verbose) {
                std::cout << *it_ << '\n';
            }

            // unpack the iterator
            bool is_measurement_ = (*it_).is_measurement_;
            idx m_ip_ = (*it_).m_ip_;
            idx q_ip_ = (*it_++).q_ip_; // post-increment the iterator

            // we have a measurement step
            if (is_measurement_) {
                std::vector<idx> target_rel_pos =
                    get_relative_pos_(qcd_.get_measurements()[m_ip_].target_);

                std::vector<idx> resZ;
                double probZ;

                idx mres = 0;
                std::vector<double> probs;
                std::vector<cmat> states;

                switch (qcd_.get_measurements()[m_ip_].measurement_type_) {
                case QCircuitDescription::MeasureType::NONE:
                    break;
                case QCircuitDescription::MeasureType::MEASURE_Z:
                    std::tie(resZ, probZ, psi_) =
                        measure_seq(psi_, target_rel_pos, qcd_.get_d());
                    dits_[qcd_.get_measurements()[m_ip_].c_reg_] = resZ[0];
                    probs_[qcd_.get_measurements()[m_ip_].c_reg_] = probZ;
                    set_measured_(qcd_.get_measurements()[m_ip_].target_[0]);
                    break;
                case QCircuitDescription::MeasureType::MEASURE_V:
                    std::tie(mres, probs, states) =
                        measure(psi_, qcd_.get_measurements()[m_ip_].mats_[0],
                                target_rel_pos, qcd_.get_d());
                    psi_ = states[mres];
                    dits_[qcd_.get_measurements()[m_ip_].c_reg_] = mres;
                    probs_[qcd_.get_measurements()[m_ip_].c_reg_] = probs[mres];
                    set_measured_(qcd_.get_measurements()[m_ip_].target_[0]);
                    break;
                case QCircuitDescription::MeasureType::MEASURE_V_MANY:
                    std::tie(mres, probs, states) =
                        measure(psi_, qcd_.get_measurements()[m_ip_].mats_[0],
                                target_rel_pos, qcd_.get_d());
                    psi_ = states[mres];
                    dits_[qcd_.get_measurements()[m_ip_].c_reg_] = mres;
                    probs_[qcd_.get_measurements()[m_ip_].c_reg_] = probs[mres];
                    for (auto&& i : qcd_.get_measurements()[m_ip_].target_)
                        set_measured_(i);
                    break;
                }  // end switch on measurement type
            }      // end if measurement step
            else { // we have a gate step
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
                            psi_ = apply(
                                psi_,
                                powm(qcd_.get_gates()[q_ip_].gate_, first_dit),
                                target_rel_pos, qcd_.get_d());
                        }
                    }
                    break;
                } // end switch on gate type
            }     // end else gate step
        }         // end main for loop
    }             // end QCircuit::run(idx, bool)
};                /* class QCircuit */

} /* namespace qpp */

#endif /* CLASSES_CIRCUITS_H_ */
