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
 * \file classes/circuits/circuits.h
 * \brief Qudit quantum circuits
 */

#ifndef CLASSES_CIRCUITS_CIRCUITS_H_
#define CLASSES_CIRCUITS_CIRCUITS_H_

namespace qpp {
/**
 * \class qpp::QCircuit
 * \brief Quantum circuit description
 * \see qpp::QEngine
 */
class QCircuit : public IDisplay, public IJSON {
    friend class QEngine;

    idx nq_;                     ///< number of qudits
    idx nc_;                     ///< number of classical "dits"
    idx d_;                      ///< qudit dimension
    std::string name_;           ///< optional circuit name
    std::vector<bool> measured_; ///< keeps track of the measured qudits

    std::unordered_map<std::size_t, cmat>
        cmat_hash_tbl_{}; ///< hash table with the matrices used in the circuit,
                          ///< with [Key = std::size_t, Value = cmat]
    std::unordered_map<std::string, idx> count_{}; ///< gate counts
    std::unordered_map<std::string, idx>
        measurement_count_{}; ///< measurement counts

    /**
     * \brief Adds matrix to the hash table
     *
     * \note Throws if a hash collision is detected., i.e., if two different
     * matrices have the same hash
     *
     * \param U Complex matrix
     * \param hashU Hash value of U
     */
    void add_hash_(const cmat& U, std::size_t hashU) {
        // EXCEPTION CHECKS

        auto search = cmat_hash_tbl_.find(hashU);
        static internal::EqualEigen equal_eigen;
        if (search != cmat_hash_tbl_.end()) // found the hash in the table
        {
            // have a hash collision
            if (!equal_eigen(search->second, U))
                throw exception::CustomException("qpp::QCircuit::add_hash_()",
                                                 "Matrix hash collision");
        }
        // END EXCEPTION CHECKS
        cmat_hash_tbl_.insert({hashU, U});
    }

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

        SINGLE_CTRL_SINGLE_TARGET, ///< controlled 1 qudit unitary gate with
                                   ///< one control and one target

        SINGLE_CTRL_MULTIPLE_TARGET, ///< controlled 1 qudit unitary gate with
                                     ///< one control and multiple targets

        MULTIPLE_CTRL_SINGLE_TARGET, ///< controlled 1 qudit unitary gate with
                                     ///< multiple controls and single target

        MULTIPLE_CTRL_MULTIPLE_TARGET, ///< controlled 1 qudit unitary gate with
                                       ///< multiple controls and multiple
                                       ///< targets

        CUSTOM_CTRL, ///< custom controlled gate with multiple controls
                     ///< and multiple targets

        SINGLE_cCTRL_SINGLE_TARGET, ///< controlled 1 qudit unitary gate with
                                    ///< one classical control and one target

        SINGLE_cCTRL_MULTIPLE_TARGET, ///< controlled 1 qudit unitary gate with
                                      ///< one classical control and multiple
                                      ///< targets

        MULTIPLE_cCTRL_SINGLE_TARGET, ///< controlled 1 qudit unitary gate with
                                      ///< multiple classical controls and
                                      ///< single target

        MULTIPLE_cCTRL_MULTIPLE_TARGET, ///< controlled 1 qudit unitary gate
                                        ///< with multiple classical controls
                                        ///< and multiple targets

        CUSTOM_cCTRL, ///< custom controlled gate with multiple classical
                      ///< controls and multiple targets
    };

    /**
     * \brief Extraction operator overload for qpp::QCircuit::GateType enum
     * class
     *
     * \param os Output stream
     * \param gate_type qpp::QCircuit::GateType enum class
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
     * \brief One step consisting only of gates/operators in the circuit
     */
    struct GateStep {
        GateType gate_type_ = GateType::NONE; ///< gate type
        std::size_t gate_hash_{};             ///< gate hash
        std::vector<idx> ctrl_{};             ///< control
        std::vector<idx> target_{}; ///< target where the gate is applied
        std::vector<idx> shift_{};  ///< shifts in CTRL gates
        std::string name_{};        ///< custom name of the gate(s)

        /**
         * \brief Default constructor
         */
        GateStep() = default;

        /**
         * \brief Constructs a gate step instance
         *
         * \param gate_type Gate type
         * \param gate_hash Hash of the quantum gate
         * \param ctrl Control qudit indexes
         * \param target Target qudit indexes
         * \param name Optional gate name
         */
        explicit GateStep(GateType gate_type, std::size_t gate_hash,
                          const std::vector<idx>& ctrl,
                          const std::vector<idx>& target,
                          const std::vector<idx>& shift = {},
                          std::string name = {})
            : gate_type_{gate_type}, gate_hash_{gate_hash}, ctrl_{ctrl},
              target_{target}, shift_{shift}, name_{name} {}
    };

    /**
     * \brief Extraction operator overload for qpp::QCircuit::GateStep class
     *
     * \param os Output stream
     * \param gate_step qpp::QCircuit::GateStep class
     * \return Output stream
     */
    friend std::ostream& operator<<(std::ostream& os,
                                    const GateStep& gate_step) {
        os << gate_step.gate_type_ << ", ";
        if (gate_step.gate_type_ >= GateType::SINGLE_cCTRL_SINGLE_TARGET)
            os << "c_ctrl = " << disp(gate_step.ctrl_, ", ") << ", ";
        else if (gate_step.gate_type_ >= GateType::SINGLE_CTRL_SINGLE_TARGET)
            os << "ctrl = " << disp(gate_step.ctrl_, ", ") << ", ";
        os << "target = " << disp(gate_step.target_, ", ") << ", ";
        if (!gate_step.shift_.empty())
            os << "shift = " << disp(gate_step.shift_, ", ") << ", ";
        os << "name = " << '\"' << gate_step.name_ << '\"';

        return os;
    }

    /**
     * \brief Type of measurement being executed in a measurement step
     */
    enum class MeasureType {
        NONE, ///< represents no measurement

        MEASURE_Z, ///< Z measurement of single qudit

        MEASURE_Z_MANY, ///< Z measurement of multiple qudit

        MEASURE_V, ///< measurement of single qudit in the orthonormal basis
                   ///< or rank-1 projectors specified by the columns of matrix
                   ///< \a V

        MEASURE_V_MANY, ///< joint measurement of multiple qudits in the
                        ///< orthonormal basis or rank-1 projectors specified
                        ///< by the columns of the matrix \a V

        MEASURE_Z_ND, ///< Z measurement of single qudit, non-destructive

        MEASURE_Z_MANY_ND, ///< Z measurement of multiple qudit, non-destructive

        MEASURE_V_ND, ///< measurement of single qudit in the orthonormal basis
                      ///< or rank-1 projectors specified by the columns of
                      ///< matrix \a V, non-destructive

        MEASURE_V_MANY_ND, ///< joint measurement of multiple qudits in the
                           ///< orthonormal basis or rank-1 projectors specified
                           ///< by the columns of the matrix \a V,
                           ///< non-destructive

        RESET, ///< resets single qudit

        RESET_MANY, ///< resets multiple qudits

        DISCARD, ///< discards single qudit

        DISCARD_MANY, ///< discards multiple qudits
    };

    /**
     * \brief Extraction operator overload for qpp::QCircuit::MeasureType enum
     * class
     *
     * \param os Output stream
     * \param measure_type qpp::QCircuit::MeasureType enum class
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
            case MeasureType::MEASURE_Z_MANY:
                os << "MEASURE_Z_MANY";
                break;
            case MeasureType::MEASURE_V:
                os << "MEASURE_V";
                break;
            case MeasureType::MEASURE_V_MANY:
                os << "MEASURE_V_MANY";
                break;
            case MeasureType::MEASURE_Z_ND:
                os << "MEASURE_Z_ND";
                break;
            case MeasureType::MEASURE_Z_MANY_ND:
                os << "MEASURE_Z_MANY_ND";
                break;
            case MeasureType::MEASURE_V_ND:
                os << "MEASURE_V_ND";
                break;
            case MeasureType::MEASURE_V_MANY_ND:
                os << "MEASURE_V_MANY_ND";
                break;
            case MeasureType::RESET:
                os << "RESET";
                break;
            case MeasureType::RESET_MANY:
                os << "RESET_MANY";
                break;
            case MeasureType::DISCARD:
                os << "DISCARD";
                break;
            case MeasureType::DISCARD_MANY:
                os << "DISCARD_MANY";
                break;
        }

        return os;
    }

    /**
     * \brief One step consisting only of measurements in the circuit
     */
    struct MeasureStep {
        MeasureType measurement_type_ = MeasureType::NONE; ///< measurement type
        std::vector<std::size_t> mats_hash_{}; ///< hashes of measurement
                                               ///< matrix/matrices
        std::vector<idx> target_{}; ///< target where the measurement is applied
        idx c_reg_{}; ///< index of the classical register where the measurement
                      ///< result is being stored
        std::string name_{}; ///< custom name of the measurement(s)

        /**
         * \brief Default constructor
         */
        MeasureStep() = default;

        /**
         * \brief Constructs a measurement step instance
         *
         * \param measurement_type Measurement type
         * \param mats_hash Vector of hashes of the measurement matrix/matrices
         * \param target Target qudit indexes
         * \param c_reg Classical register where the value of the measurement is
         * stored
         * \param name Optional gate name
         */
        explicit MeasureStep(MeasureType measurement_type,
                             const std::vector<std::size_t>& mats_hash,
                             const std::vector<idx>& target, idx c_reg,
                             std::string name = {})
            : measurement_type_{measurement_type}, mats_hash_{mats_hash},
              target_{target}, c_reg_{c_reg}, name_{name} {}
    };

    /**
     * \brief Extraction operator overload for qpp::QCircuit::MeasureStep class
     *
     * \param os Output stream
     * \param measure_step qpp::QCircuit::MeasureStep enum class
     * \return Output stream
     */
    friend std::ostream& operator<<(std::ostream& os,
                                    const MeasureStep& measure_step) {
        os << measure_step.measurement_type_ << ", ";
        os << "target = " << disp(measure_step.target_, ", ") << ", ";
        if (measure_step.measurement_type_ != MeasureType::RESET &&
            measure_step.measurement_type_ != MeasureType::RESET_MANY &&
            measure_step.measurement_type_ != MeasureType::DISCARD &&
            measure_step.measurement_type_ != MeasureType::DISCARD_MANY)
            os << "c_reg = " << measure_step.c_reg_ << ", ";
        os << "name = " << '\"' << measure_step.name_ << '\"';

        return os;
    }

    /**
     * \brief Types of each step in the quantum circuit description
     */
    enum class StepType {
        NONE,        ///< represents no step
        GATE,        ///< quantum gate(s)
        MEASUREMENT, ///< measurement
        NOP,         ///< no-op
    };

  private:
    std::vector<GateStep> gates_{};           ///< gates
    std::vector<MeasureStep> measurements_{}; ///< measurements
    std::vector<StepType> step_types_{};      ///< type of each step

    /**
     * \brief Vector of qpp::QCircuit::MeasureStep
     *
     * \return Vector of qpp::QCircuit::MeasureStep
     */
    const std::vector<MeasureStep>& get_measurements_() const noexcept {
        return measurements_;
    }

    /**
     * \brief Vector of qpp::QCircuit::GateStep
     *
     * \return Vector of qpp::QCircuit::GateStep
     */
    const std::vector<GateStep>& get_gates_() const noexcept { return gates_; }

    /**
     * \brief Hash table with the matrices used in the circuit
     *
     * \return Hash table with the matrices used in the circuit
     */
    const std::unordered_map<std::size_t, cmat>& get_cmat_hash_tbl_() const
        noexcept {
        return cmat_hash_tbl_;
    }

  public:
    /**
     * \class qpp::QCircuit::iterator
     * \brief Quantum circuit description bound-checking (safe) iterator
     *
     * \note The iterator is a const_iterator by default
     */
    class iterator {
        ///< non-owning pointer to the parent const quantum circuit description
        const QCircuit* qc_{nullptr};

        /**
         * \class qpp::QCircuit::iterator::value_type_
         * \brief Value type class for qpp::QCircuit::iterator
         */
        class value_type_ : public IDisplay {
          public:
            ///< non-owning pointer to the grand-parent const quantum circuit
            ///< description
            const QCircuit* value_type_qc_;

            StepType type_{StepType::NONE}; ///< step type
            idx ip_{static_cast<idx>(-1)};  ///< instruction pointer
            std::vector<GateStep>::const_iterator
                gates_ip_{}; ///< gates instruction pointer
            std::vector<MeasureStep>::const_iterator
                measurements_ip_{}; ///< measurements instruction pointer

            /**
             * \brief Constructor
             *
             * \param value_type_qc Pointer to constant quantum circuit
             * description
             */
            explicit value_type_(const QCircuit* value_type_qc)
                : value_type_qc_{value_type_qc} {}

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

          private:
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
                // field spacing for the step number
                idx text_width =
                    std::to_string(value_type_qc_->get_step_count()).size() + 1;

                // gate step
                if (type_ == StepType::GATE) {
                    os << std::left;
                    os << std::setw(text_width) << ip_;
                    os << std::right;
                    idx pos = std::distance(std::begin(value_type_qc_->gates_),
                                            gates_ip_);
                    os << value_type_qc_->get_gates_()[pos];
                }
                // measurement step
                else if (type_ == StepType::MEASUREMENT) {
                    os << std::left;
                    os << std::setw(text_width) << ip_;
                    os << std::right;
                    idx pos =
                        std::distance(std::begin(value_type_qc_->measurements_),
                                      measurements_ip_);

                    switch (
                        value_type_qc_->measurements_[pos].measurement_type_) {
                        case MeasureType::NONE:
                            break;
                        case MeasureType::MEASURE_Z:
                        case MeasureType::MEASURE_Z_MANY:
                        case MeasureType::MEASURE_V:
                        case MeasureType::MEASURE_V_MANY:
                            os << "|> ";
                            break;
                        case MeasureType::MEASURE_Z_ND:
                        case MeasureType::MEASURE_Z_MANY_ND:
                        case MeasureType::MEASURE_V_ND:
                        case MeasureType::MEASURE_V_MANY_ND:
                            os << "|] ";
                            break;
                        case MeasureType::RESET:
                        case MeasureType::RESET_MANY:
                            os << "|* ";
                            break;
                        case MeasureType::DISCARD:
                        case MeasureType::DISCARD_MANY:
                            os << "|x ";
                            break;
                    } /* end switch */

                    os << value_type_qc_->get_measurements_()[pos];
                }
                // no-op
                else if (type_ == StepType::NOP) {
                    os << std::left;
                    os << std::setw(text_width) << ip_;
                    os << std::right;
                    os << "NOP";
                }
                // otherwise
                else {
                }

                return os;
            }
        }; /* class value_type_ */

        value_type_ elem_{nullptr};

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
            if (qc_ == nullptr) {
                throw exception::InvalidIterator(
                    "qpp::QCircuit::iterator::operator++()");
            }

            // protects against incrementing an empty circuit iterator
            if (qc_->get_step_count() == 0) {
                throw exception::InvalidIterator(
                    "qpp::QCircuit::iterator::operator++()");
            }

            // protects against incrementing past the end
            if (elem_.ip_ == qc_->get_step_count()) {
                throw exception::InvalidIterator(
                    "qpp::QCircuit::iterator::operator++()");
            }
            // END EXCEPTION CHECKS

            // gate step
            if (elem_.type_ == StepType::GATE) {
                std::advance(elem_.gates_ip_, 1);
            }
            // measurement step
            else if (elem_.type_ == StepType::MEASUREMENT) {
                std::advance(elem_.measurements_ip_, 1);
            }
            // no-op
            else if (elem_.type_ == StepType::NOP) {
            }
            // otherwise
            else {
            }

            // increment the instruction pointer
            ++elem_.ip_;

            // if we hit the end
            if (elem_.ip_ == qc_->get_step_count()) {
                elem_.type_ = StepType::NONE;
            } else {
                // set the next step type
                elem_.type_ = qc_->step_types_[elem_.ip_];
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
            return std::tie(elem_.type_, elem_.ip_, elem_.gates_ip_,
                            elem_.measurements_ip_) ==
                   std::tie(rhs.elem_.type_, rhs.elem_.ip_, rhs.elem_.gates_ip_,
                            rhs.elem_.measurements_ip_);
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

            // protects against de-referencing past the last element or against
            // de-referencing invalid iterators
            if (qc_ == nullptr || elem_.ip_ == qc_->get_step_count())
                throw exception::InvalidIterator(
                    "qpp::QCircuit::iterator::operator*()");
            // END EXCEPTION CHECKS

            return elem_;
        }

        /**
         * \brief Sets the iterator to std::begin(this)
         *
         * \param qc Pointer to constant quantum circuit description
         */
        void set_begin_(const QCircuit* qc) {
            qc_ = qc;
            elem_ = value_type_{qc_};

            if (qc_ != nullptr) {
                if (qc_->get_step_count() != 0) // non-empty circuit
                {
                    elem_.type_ = qc_->step_types_[0];
                    elem_.ip_ = 0;
                }
                elem_.gates_ip_ = std::begin(qc_->gates_);
                elem_.measurements_ip_ = std::begin(qc_->measurements_);
            }
        }

        /**
         * \brief Sets the iterator to std::begin(this)
         *
         * \param qc Pointer to constant quantum circuit description
         */
        void set_end_(const QCircuit* qc) {
            qc_ = qc;
            elem_ = value_type_{qc_};

            if (qc_ != nullptr) {
                if (qc->get_step_count() != 0) {
                    elem_.ip_ = qc->get_step_count();
                }
                elem_.gates_ip_ = std::end(qc->gates_);
                elem_.measurements_ip_ = std::end(qc->measurements_);
            }
        }

        // iterator traits
        using difference_type = ptrdiff_t;                   ///< iterator trait
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
        it.set_begin_(this);

        return it;
    }

    /**
     * \brief Constant iterator to the first element
     *
     * \return Constant iterator to the first element
     */
    const_iterator begin() const noexcept {
        iterator it;
        it.set_begin_(this);

        return it;
    }

    /**
     * \brief Constant iterator to the first element
     *
     * \return Constant iterator to the first element
     */
    const_iterator cbegin() const noexcept {
        iterator it;
        it.set_begin_(this);

        return it;
    }

    /**
     * \brief Iterator to the next to the last element
     *
     * \return Iterator to the next to the last element
     */
    iterator end() {
        iterator it;
        it.set_end_(this);

        return it;
    }

    /**
     * \brief Constant iterator to the next to the last element
     *
     * \return Constant iterator to the next to the last element
     */
    const_iterator end() const noexcept {
        iterator it;
        it.set_end_(this);

        return it;
    }

    /**
     * \brief Constant iterator to the next to the last element
     *
     * \return Constant iterator to the next to the last element
     */
    const_iterator cend() const noexcept {
        iterator it;
        it.set_end_(this);

        return it;
    }

    /**
     * \brief Constructs a quantum circuit description
     *
     * \note The measurement results can only be stored in the classical dits of
     * which number is specified by \a nc
     *
     * \param nq Number of qudits
     * \param nc Number of classical dits (optional)
     * \param d Subsystem dimensions (optional, default is qubit, i.e. \a d = 2)
     * \param name Circuit name (optional)
     */
    explicit QCircuit(idx nq, idx nc = 0, idx d = 2, std::string name = {})
        : nq_{nq}, nc_{nc}, d_{d}, name_{name}, measured_(nq, false) {
        // EXCEPTION CHECKS

        // if (nq == 0)
        //    throw exception::ZeroSize("qpp::QCircuit::QCircuit()");
        if (d < 2)
            throw exception::OutOfRange("qpp::QCircuit::QCircuit()");
        // END EXCEPTION CHECKS
    }

    /**
     * \brief Default virtual destructor
     */
    virtual ~QCircuit() = default;

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
    std::string get_name() const { return name_; }
    /**
     * \brief Check whether qudit \a i was already measured (destructively)
     * \param i Qudit index
     * \return True if qudit \a i was already measured (destructively), false
     * othwewise
     */
    idx get_measured(idx i) const {
        // EXCEPTION CHECKS

        if (i >= nq_)
            throw exception::OutOfRange("qpp::QCircuit::get_measured()");
        // END EXCEPTION CHECKS

        return measured_[i];
    }

    /**
     * \brief Vector of already measured (destructively) qudit indexes
     *
     * \return Vector of already measured (destructively) qudit indexes
     */
    std::vector<idx> get_measured() const {
        std::vector<idx> result;
        for (idx i = 0; i < nq_; ++i)
            if (get_measured(i))
                result.emplace_back(i);

        return result;
    }

    /**
     * \brief Vector of non-measured (destructively) qudit indexes
     *
     * \return Vector of non-measured (destructively) qudit indexes
     */
    std::vector<idx> get_non_measured() const {
        std::vector<idx> result;
        for (idx i = 0; i < nq_; ++i)
            if (!get_measured(i))
                result.emplace_back(i);

        return result;
    }

    /**
     * \brief Quantum circuit description gate count
     *
     * \param name Gate name
     * \return Gate count
     */
    idx get_gate_count(const std::string& name) const {
        idx result = 0;

        // name not found in the hash table
        try {
            result = count_.at(name);
        } catch (...) {
            return 0;
        }

        return result;
    }

    /**
     * \brief Quantum circuit description total gate count
     *
     * \return Total gate count
     */
    idx get_gate_count() const {
        idx result = 0;

        for (auto&& elem : count_)
            result += elem.second;

        return result;
    }

    /**
     * \brief Quantum circuit description gate depth
     *
     * \param name Gate name
     * \return Gate depth
     */
    idx get_gate_depth(const std::string& name) const {
        bool found = false;
        std::vector<idx> heights(nc_ + nq_, 0);

        // iterate over all steps in the circuit
        for (auto&& step : *this) {
            // gates
            if (step.type_ == StepType::GATE) {
                GateStep gate_step = *step.gates_ip_;

                if (name != __FILE__ "__total_gate_depth__" &&
                    gate_step.name_ != name)
                    continue; // we skip this gate step

                found = true; // gate was found in the circuit

                std::vector<idx> ctrl = gate_step.ctrl_;
                std::vector<idx> target = gate_step.target_;
                std::vector<idx> ctrl_target;
                ctrl_target.reserve(ctrl.size() + target.size());
                ctrl_target.insert(ctrl_target.end(), ctrl.begin(), ctrl.end());
                ctrl_target.insert(ctrl_target.end(), target.begin(),
                                   target.end());

                idx max_height = 0;
                switch (gate_step.gate_type_) {
                    case GateType::NONE:
                    case GateType::SINGLE:
                    case GateType::TWO:
                    case GateType::THREE:
                    case GateType::CUSTOM:
                    case GateType::FAN:
                    case GateType::SINGLE_CTRL_SINGLE_TARGET:
                    case GateType::SINGLE_CTRL_MULTIPLE_TARGET:
                    case GateType::MULTIPLE_CTRL_SINGLE_TARGET:
                    case GateType::MULTIPLE_CTRL_MULTIPLE_TARGET:
                    case GateType::CUSTOM_CTRL:
                        // compute the "height" of the to-be-placed gate
                        for (auto&& i : ctrl_target)
                            if (heights[nc_ + i] > max_height)
                                max_height = heights[nc_ + i];
                        // apply (ctrl) gate
                        for (auto&& i : ctrl_target)
                            heights[nc_ + i] = max_height + 1;
                        break;
                    case GateType::SINGLE_cCTRL_SINGLE_TARGET:
                    case GateType::SINGLE_cCTRL_MULTIPLE_TARGET:
                    case GateType::MULTIPLE_cCTRL_SINGLE_TARGET:
                    case GateType::MULTIPLE_cCTRL_MULTIPLE_TARGET:
                    case GateType::CUSTOM_cCTRL:
                        // compute the "height" of the to-be-placed gate
                        for (auto&& i : ctrl)
                            if (heights[i] > max_height)
                                max_height = heights[i];
                        for (auto&& i : target)
                            if (heights[nc_ + i] > max_height)
                                max_height = heights[nc_ + i];
                        // apply classical ctrl
                        for (auto&& i : ctrl)
                            heights[i] = max_height + 1;
                        // apply gate
                        for (auto&& i : target)
                            heights[nc_ + i] = max_height + 1;
                        break;
                } // end switch
            }     // end if (step.type_ == StepType::GATE)
        }         // end for

        return found ? *std::max_element(std::begin(heights), std::end(heights))
                     : 0;
    }

    /**
     * \brief Quantum circuit description total gate depth
     *
     * \return Total gate depth
     */
    idx get_gate_depth() const {
        return get_gate_depth(__FILE__ "__total_gate_depth__");
    }

    /**
     * \brief Quantum circuit description measurement depth
     *
     * \param name Measurement name
     * \return Measurement depth
     */
    idx get_measurement_depth(const std::string& name) const {
        bool found = false;
        std::vector<idx> heights(nc_ + nq_, 0);

        // iterate over all steps in the circuit
        for (auto&& step : *this) {
            // measurements
            if (step.type_ == StepType::MEASUREMENT) {
                MeasureStep measure_step = *step.measurements_ip_;
                if (name != __FILE__ "__total_measurement_depth__" &&
                    measure_step.name_ != name)
                    continue; // we skip this measurement step

                found = true; // measurement was found in the circuit

                std::vector<idx> target = measure_step.target_;
                idx c_reg = measure_step.c_reg_;

                idx max_height = 0;
                switch (measure_step.measurement_type_) {
                    case MeasureType::NONE:
                    case MeasureType::MEASURE_Z:
                    case MeasureType::MEASURE_Z_MANY:
                    case MeasureType::MEASURE_V:
                    case MeasureType::MEASURE_V_MANY:
                    case MeasureType::MEASURE_Z_ND:
                    case MeasureType::MEASURE_Z_MANY_ND:
                    case MeasureType::MEASURE_V_ND:
                    case MeasureType::MEASURE_V_MANY_ND:
                        // compute the "height" of the to-be-placed measurement
                        if (heights[c_reg] > max_height)
                            max_height = heights[c_reg];
                        for (auto&& i : target)
                            if (heights[nc_ + i] > max_height)
                                max_height = heights[nc_ + i];
                        // apply measurement
                        heights[c_reg] = max_height + 1;
                        for (auto&& i : target)
                            heights[nc_ + i] = max_height + 1;
                        break;
                    case MeasureType::RESET:
                    case MeasureType::RESET_MANY:
                    case MeasureType::DISCARD:
                    case MeasureType::DISCARD_MANY:
                        for (auto&& i : target)
                            if (heights[nc_ + i] > max_height)
                                max_height = heights[nc_ + i];
                        // apply reset/discard
                        for (auto&& i : target)
                            heights[nc_ + i] = max_height + 1;
                        break;
                } // end switch
            }     // if (step.type_ == StepType::MEASUREMENT) }
        }         // end for

        return found ? *std::max_element(std::begin(heights), std::end(heights))
                     : 0;
    }

    /**
     * \brief Quantum circuit description total measurement depth
     *
     * \return Total measurement depth
     */
    idx get_measurement_depth() const {
        return get_measurement_depth(__FILE__ "__total_measurement_depth__");
    }

    // computes the depth greedily, measuring the "height" (depth) of the
    // "pieces" (gates) placed in a Tetris-like style
    /**
     * \brief Quantum circuit description total depth
     *
     * \param name Gate/measurement name (optional)
     * \return Gate/measurement depth
     */
    idx get_depth() const { return get_gate_depth() + get_measurement_depth(); }

    /**
     * \brief Quantum circuit description measurement count
     *
     * \param name Measurement name
     * \return Measurement count
     */
    idx get_measurement_count(const std::string& name) const {
        idx result = 0;

        // name not found in the hash table
        try {
            result = measurement_count_.at(name);
        } catch (...) {
            return 0;
        }

        return result;
    }

    /**
     * \brief Quantum circuit description total measurement count
     *
     * \return Total measurement count
     */
    idx get_measurement_count() const {
        idx result = 0;

        for (auto&& elem : measurement_count_)
            result += elem.second;

        return result;
    }

    /**
     * \brief Quantum circuit total steps count, i.e. the sum of gate count and
     * measurement count
     *
     * \return Total (gates + measurements) count
     */
    idx get_step_count() const noexcept { return step_types_.size(); }

    /**
     * \brief No-op count
     *
     * \return No-op count
     */
    idx get_nop_count() const {
        return std::count_if(
            std::begin(step_types_), std::end(step_types_),
            [](const StepType& s) { return s == StepType::NOP; });
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
     * \brief Adds \a n additional qudits before qudit \a i (by default adds
     * them at the end)
     *
     * \note Qudits with indexes greater or equal than the newly inserted ones
     * have their indexes automatically incremented
     *
     * \param n Number of qudits
     * \param i Qudit index
     * \return Reference to the current instance
     */
    QCircuit& add_qudit(idx n = 1, idx i = -1) {
        // EXCEPTION CHECKS

        if (i == static_cast<idx>(-1)) {
            i = nq_;
            measured_.insert(std::end(measured_), n, false);
        } else if (i > nq_)
            throw exception::OutOfRange("qpp::QCircuit::add_qudit()");
        // END EXCEPTION CHECKS

        nq_ += n;

        // updated the measured qudits
        measured_.insert(std::next(std::begin(measured_), i), n, false);

        // update gate indexes
        for (auto& gate : gates_) {
            for (auto& pos : gate.ctrl_) {
                if (pos >= i)
                    pos += n;
            }
            for (auto& pos : gate.target_) {
                if (pos >= i)
                    pos += n;
            }
        }

        // update measurement indexes
        for (auto& measurement : measurements_) {
            for (auto& pos : measurement.target_) {
                if (pos >= i)
                    pos += n;
            }
        }

        return *this;
    }

    /**
     * \brief Adds \a n additional classical dits before dit \a i (by default
     * adds them at the end)
     *
     * \note Classical dits with indexes greater or equal than the newly
     * inserted ones have their indexes automatically incremented
     *
     * \param n Number of classical dits
     * \param i Classical dit index
     * \return Reference to the current instance
     */
    QCircuit& add_dit(idx n = 1, idx i = -1) {
        // EXCEPTION CHECKS

        if (i == static_cast<idx>(-1))
            i = nc_;
        else if (i > nc_)
            throw exception::OutOfRange("qpp::QCircuit::add_dit()");
        // END EXCEPTION CHECKS

        nc_ += n;

        // update gate indexes
        for (auto& gate : gates_) {
            switch (gate.gate_type_) {
                // classically-controlled gates
                case GateType::SINGLE_cCTRL_SINGLE_TARGET:
                case GateType::SINGLE_cCTRL_MULTIPLE_TARGET:
                case GateType::MULTIPLE_cCTRL_SINGLE_TARGET:
                case GateType::MULTIPLE_cCTRL_MULTIPLE_TARGET:
                case GateType::CUSTOM_cCTRL:
                    for (auto& pos : gate.ctrl_) {
                        if (pos >= i)
                            pos += n;
                    }
                    break;
                default:
                    break;
            }
        }

        // update measurement indexes
        for (auto& measurement : measurements_) {
            if (measurement.c_reg_ >= i)
                measurement.c_reg_ += n;
        }

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
    QCircuit& gate(const cmat& U, idx i, std::string name = {}) {
        // EXCEPTION CHECKS

        try {
            // check valid target
            if (i >= nq_)
                throw exception::OutOfRange("qpp::QCircuit::gate()");
            // check not measured before
            if (get_measured(i))
                throw exception::QuditAlreadyMeasured("qpp::QCircuit::gate()");

            // check square matrix for the gate
            if (!internal::check_square_mat(U))
                throw exception::MatrixNotSquare("qpp::QCircuit::gate()");
            // check correct dimension
            if (static_cast<idx>(U.rows()) != d_)
                throw exception::DimsMismatchMatrix("qpp::QCircuit::gate()");
        } catch (exception::Exception&) {
            std::cerr << "At STEP " << get_step_count() << "\n";
            throw;
        }
        // END EXCEPTION CHECKS

        if (name.empty())
            name = qpp::Gates::get_instance().get_name(U);
        std::size_t hashU = hash_eigen(U);
        add_hash_(U, hashU);
        gates_.emplace_back(GateType::SINGLE, hashU, std::vector<idx>{},
                            std::vector<idx>{i}, std::vector<idx>{}, name);
        step_types_.emplace_back(StepType::GATE);
        ++count_[name];

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
    QCircuit& gate(const cmat& U, idx i, idx j, std::string name = {}) {
        // EXCEPTION CHECKS

        try {
            // check valid target
            if (i >= nq_ || j >= nq_ || i == j)
                throw exception::OutOfRange("qpp::QCircuit::gate()");
            if (get_measured(i) || get_measured(j))
                throw exception::QuditAlreadyMeasured("qpp::QCircuit::gate()");

            // check square matrix for the gate
            if (!internal::check_square_mat(U))
                throw exception::MatrixNotSquare("qpp::QCircuit::gate()");
            // check correct dimension
            if (static_cast<idx>(U.rows()) != d_ * d_)
                throw exception::DimsMismatchMatrix("qpp::QCircuit::gate()");
        } catch (exception::Exception&) {
            std::cerr << "At STEP " << get_step_count() << "\n";
            throw;
        }
        // END EXCEPTION CHECKS

        if (name.empty())
            name = qpp::Gates::get_instance().get_name(U);
        std::size_t hashU = hash_eigen(U);
        add_hash_(U, hashU);
        gates_.emplace_back(GateType::TWO, hashU, std::vector<idx>{},
                            std::vector<idx>{i, j}, std::vector<idx>{}, name);
        step_types_.emplace_back(StepType::GATE);
        ++count_[name];

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
    QCircuit& gate(const cmat& U, idx i, idx j, idx k, std::string name = {}) {
        // EXCEPTION CHECKS

        try {
            // check valid target
            if (i >= nq_ || j >= nq_ || k >= nq_ || (i == j) || (i == k) ||
                (j == k))
                throw exception::OutOfRange("qpp::QCircuit::gate()");
            if (get_measured(i) || get_measured(j) || get_measured(k))
                throw exception::QuditAlreadyMeasured("qpp::QCircuit::gate()");

            // check square matrix for the gate
            if (!internal::check_square_mat(U))
                throw exception::MatrixNotSquare("qpp::QCircuit::gate()");
            // check correct dimension
            if (static_cast<idx>(U.rows()) != d_ * d_ * d_)
                throw exception::DimsMismatchMatrix("qpp::QCircuit::gate()");
        } catch (exception::Exception&) {
            std::cerr << "At STEP " << get_step_count() << "\n";
            throw;
        }
        // END EXCEPTION CHECKS

        if (name.empty())
            name = qpp::Gates::get_instance().get_name(U);
        std::size_t hashU = hash_eigen(U);
        add_hash_(U, hashU);
        gates_.emplace_back(GateType::THREE, hashU, std::vector<idx>{},
                            std::vector<idx>{i, j, k}, std::vector<idx>{},
                            name);
        step_types_.emplace_back(StepType::GATE);
        ++count_[name];

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
                       std::string name = {}) {
        // EXCEPTION CHECKS

        try {
            // check valid target
            if (target.empty())
                throw exception::ZeroSize("qpp::QCircuit::gate_fan()");
            for (auto&& elem : target) {
                if (elem >= nq_)
                    throw exception::OutOfRange("qpp::QCircuit::gate_fan()");
                // check target was not measured before
                if (get_measured(elem))
                    throw exception::QuditAlreadyMeasured(
                        "qpp::QCircuit::gate_fan()");
            }
            // check no duplicates target
            if (!internal::check_no_duplicates(target))
                throw exception::Duplicates("qpp::QCircuit::gate_fan()");

            // check square matrix for the gate
            if (!internal::check_square_mat(U))
                throw exception::MatrixNotSquare("qpp::QCircuit::gate_fan()");
            // check correct dimension
            if (static_cast<idx>(U.rows()) != d_)
                throw exception::DimsMismatchMatrix(
                    "qpp::QCircuit::gate_fan()");
        } catch (exception::Exception&) {
            std::cerr << "At STEP " << get_step_count() << "\n";
            throw;
        }
        // END EXCEPTION CHECKS

        if (name.empty())
            name = qpp::Gates::get_instance().get_name(U);
        std::size_t hashU = hash_eigen(U);
        add_hash_(U, hashU);
        gates_.emplace_back(GateType::FAN, hashU, std::vector<idx>{}, target,
                            std::vector<idx>{}, name);
        step_types_.emplace_back(StepType::GATE);
        count_[name] += target.size();

        return *this;
    }

    // std::initializer_list overload, avoids ambiguity for 2-element lists, see
    // http://stackoverflow.com
    // /questions/26750039/ambiguity-when-using-initializer-list-as-parameter
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
    QCircuit& gate_fan(const cmat& U, const std::initializer_list<idx>& target,
                       std::string name = {}) {
        return gate_fan(U, std::vector<idx>(target), name);
    }

    /**
     * \brief Applies the single qudit gate \a U on all of the remaining
     * non-measured qudits
     *
     * \param U Single qudit quantum gate
     * \param name Optional gate name
     * \return Reference to the current instance
     */
    QCircuit& gate_fan(const cmat& U, std::string name = {}) {
        // EXCEPTION CHECKS

        try {
            // check square matrix for the gate
            if (!internal::check_square_mat(U))
                throw exception::MatrixNotSquare("qpp::QCircuit::gate_fan()");
            // check correct dimension
            if (static_cast<idx>(U.rows()) != d_)
                throw exception::DimsMismatchMatrix(
                    "qpp::QCircuit::gate_fan()");
        } catch (exception::Exception&) {
            std::cerr << "At STEP " << get_step_count() << "\n";
            throw;
        }
        // END EXCEPTION CHECKS

        if (name.empty())
            name = qpp::Gates::get_instance().get_name(U);
        std::size_t hashU = hash_eigen(U);
        add_hash_(U, hashU);
        gates_.emplace_back(GateType::FAN, hashU, std::vector<idx>{},
                            get_non_measured(), std::vector<idx>{}, name);
        step_types_.emplace_back(StepType::GATE);
        count_[name] += get_non_measured().size();

        return *this;
    }

    /**
     * \brief Jointly applies the custom multiple qudit gate \a U on the qudit
     * indexes specified by \a target
     *
     * \param U Multiple qudit quantum gate
     * \param target Subsystem indexes where the gate \a U is applied
     * \param name Optional gate name
     * \return Reference to the current instance
     */
    QCircuit& gate_custom(const cmat& U, const std::vector<idx>& target,
                          std::string name = {}) {
        // EXCEPTION CHECKS

        idx n = static_cast<idx>(target.size());
        idx D = static_cast<idx>(std::llround(std::pow(d_, n)));

        try {
            // check valid target
            if (target.empty())
                throw exception::ZeroSize("qpp::QCircuit::gate_custom()");
            for (auto&& elem : target) {
                if (elem >= nq_)
                    throw exception::OutOfRange("qpp::QCircuit::gate_custom()");
                // check target was not measured before
                if (get_measured(elem))
                    throw exception::QuditAlreadyMeasured(
                        "qpp::QCircuit::gate_custom()");
            }
            // check no duplicates target
            if (!internal::check_no_duplicates(target))
                throw exception::Duplicates("qpp::QCircuit::gate_custom()");

            // check square matrix for the gate
            if (!internal::check_square_mat(U))
                throw exception::MatrixNotSquare(
                    "qpp::QCircuit::gate_custom()");
            // check correct dimension
            if (static_cast<idx>(U.rows()) != D)
                throw exception::DimsMismatchMatrix(
                    "qpp::QCircuit::gate_custom()");
        } catch (exception::Exception&) {
            std::cerr << "At STEP " << get_step_count() << "\n";
            throw;
        }
        // END EXCEPTION CHECKS

        if (name.empty())
            name = qpp::Gates::get_instance().get_name(U);
        std::size_t hashU = hash_eigen(U);
        add_hash_(U, hashU);
        gates_.emplace_back(GateType::CUSTOM, hashU, std::vector<idx>{}, target,
                            std::vector<idx>{}, name);
        step_types_.emplace_back(StepType::GATE);
        ++count_[name];

        return *this;
    }

    /**
     * \brief Applies the quantum Fourier transform (as a series of gates) on
     * the qudit indexes specified by \a target
     *
     * \param target Subsystem indexes where the quantum Fourier transform is
     * applied
     * \param swap Swaps the qubits at the end (true by default)
     * \return Reference to the current instance
     */
    QCircuit& QFT(const std::vector<idx>& target, bool swap = true) {
        // EXCEPTION CHECKS
        try {
            // check valid target
            if (target.empty())
                throw exception::ZeroSize("qpp::QCircuit::QFT()");
            for (auto&& elem : target) {
                if (elem >= nq_)
                    throw exception::OutOfRange("qpp::QCircuit::QFT()");
                // check target was not measured before
                if (get_measured(elem))
                    throw exception::QuditAlreadyMeasured(
                        "qpp::QCircuit::QFT()");
            }
            // check no duplicates target
            if (!internal::check_no_duplicates(target))
                throw exception::Duplicates("qpp::QCircuit::QFT()");
        } catch (exception::Exception&) {
            std::cerr << "At STEP " << get_step_count() << "\n";
            throw;
        }
        // END EXCEPTION CHECKS

        idx n_subsys = target.size();
        if (d_ == 2) // qubits
        {
            for (idx i = 0; i < n_subsys; ++i) {
                // apply Hadamard on qubit i
                gate(Gates::get_instance().H, target[i]);
                // apply controlled rotations
                for (idx j = 2; j <= n_subsys - i; ++j) {
                    // construct Rj
                    cmat Rj(2, 2);
                    Rj << 1, 0, 0, exp(2.0 * pi * 1_i / std::pow(2, j));
                    CTRL(Rj, target[i + j - 1], target[i], {},
                         "CTRL-R" + std::to_string(j));
                }
            }
            if (swap) {
                // we have the qubits in reversed order, we must swap them
                for (idx i = 0; i < n_subsys / 2; ++i) {
                    gate(Gates::get_instance().SWAP, target[i],
                         target[n_subsys - i - 1]);
                }
            }

        } else { // qudits
            for (idx i = 0; i < n_subsys; ++i) {
                // apply qudit Fourier on qudit i
                gate(Gates::get_instance().Fd(d_), target[i], "Fd");
                // apply controlled rotations
                for (idx j = 2; j <= n_subsys - i; ++j) {
                    // construct Rj
                    cmat Rj = cmat::Zero(d_, d_);
                    for (idx m = 0; m < d_; ++m) {
                        Rj(m, m) = exp(2.0 * pi * m * 1_i / std::pow(d_, j));
                    }
                    CTRL(Rj, target[i + j - 1], target[i], {},
                         "CTRL-R" + std::to_string(j) + "d");
                }
            }
            if (swap) {
                // we have the qudits in reversed order, we must swap them
                for (idx i = 0; i < n_subsys / 2; ++i) {
                    gate(Gates::get_instance().SWAPd(d_), target[i],
                         target[n_subsys - i - 1], "SWAPd");
                }
            }
        }

        return *this;
    }

    // std::initializer_list overload, avoids ambiguity for {idx} -> bool
    /**
     * \brief Applies the quantum Fourier transform (as a series of gates) on
     * the qudit indexes specified by \a target
     *
     * \param target Subsystem indexes where the quantum Fourier transform is
     * applied
     * \param swap Swaps the qubits at the end (true by default)
     * \return Reference to the current instance
     */
    QCircuit& QFT(const std::initializer_list<idx>& target, bool swap = true) {
        return QFT(std::vector<idx>(target), swap);
    }

    /**
     * \brief Applies the quantum Fourier transform (as a series of gates) on
     * all of remaining non-measured qudits
     *
     * \param swap Swaps the qubits at the end (true by default)
     * \return Reference to the current instance
     */
    QCircuit& QFT(bool swap = true) { return QFT(get_non_measured(), swap); }

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
                  bool swap QPP_UNUSED_ = true) {
        // EXCEPTION CHECKS
        try {
            // check valid target
            if (target.empty())
                throw exception::ZeroSize("qpp::QCircuit::TFQ()");
            for (auto&& elem : target) {
                if (elem >= nq_)
                    throw exception::OutOfRange("qpp::QCircuit::TFQ()");
                // check target was not measured before
                if (get_measured(elem))
                    throw exception::QuditAlreadyMeasured(
                        "qpp::QCircuit::TFQ()");
            }
            // check no duplicates target
            if (!internal::check_no_duplicates(target))
                throw exception::Duplicates("qpp::QCircuit::TFQ()");
        } catch (exception::Exception&) {
            std::cerr << "At STEP " << get_step_count() << "\n";
            throw;
        }
        // END EXCEPTION CHECKS

        idx n_subsys = target.size();
        if (d_ == 2) // qubits
        {
            if (swap) {
                // we have the qubits in reversed order, we must swap them
                for (idx i = n_subsys / 2; i-- > 0;) {
                    gate(Gates::get_instance().SWAP, target[i],
                         target[n_subsys - i - 1]);
                }
            }
            for (idx i = n_subsys; i-- > 0;) {
                // apply controlled rotations
                for (idx j = n_subsys - i + 1; j-- > 2;) {
                    // construct Rj
                    cmat Rj(2, 2);
                    Rj << 1, 0, 0, exp(-2.0 * pi * 1_i / std::pow(2, j));
                    CTRL(Rj, target[i + j - 1], target[i], {},
                         "CTRL-R" + std::to_string(j) + "+");
                }
                // apply Hadamard on qubit i
                gate(Gates::get_instance().H, target[i]);
            }
        } else { // qudits
            if (swap) {
                // we have the qudits in reversed order, we must swap them
                for (idx i = n_subsys / 2; i-- > 0;) {
                    gate(Gates::get_instance().SWAPd(d_), target[i],
                         target[n_subsys - i - 1], "SWAPd");
                }
            }
            for (idx i = n_subsys; i-- > 0;) {
                // apply controlled rotations
                for (idx j = n_subsys - i + 1; j-- > 2;) {
                    // construct Rj
                    cmat Rj = cmat::Zero(d_, d_);
                    for (idx m = 0; m < d_; ++m) {
                        Rj(m, m) = exp(-2.0 * pi * m * 1_i / std::pow(d_, j));
                    }
                    CTRL(Rj, target[i + j - 1], target[i], {},
                         "CTRL-R" + std::to_string(j) + "d+");
                }
                // apply qudit Fourier on qudit i
                gate(qpp::adjoint(Gates::get_instance().Fd(d_)), target[i],
                     "Fd+");
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
    QCircuit& TFQ(bool swap = true) { return TFQ(get_non_measured(), swap); }

    // single ctrl single target
    /**
     * \brief Applies the single qudit controlled gate \a U with control qudit
     * \a ctrl and target qudit \a target
     *
     * \param U Single qudit quantum gate
     * \param ctrl Control qudit index
     * \param target Target qudit index
     * \param shift Performs the control as if the \a ctrl qudit state was
     * \f$ X\f$-incremented by \a shift
     * \param name Optional gate name
     * \return Reference to the current instance
     */
    QCircuit& CTRL(const cmat& U, idx ctrl, idx target, idx shift = 0,
                   std::string name = {}) {
        // EXCEPTION CHECKS

        try {
            // check valid ctrl and target
            if (ctrl >= nq_ || target >= nq_ || ctrl == target)
                throw exception::OutOfRange("qpp::QCircuit::CTRL()");
            if (get_measured(ctrl) || get_measured(target))
                throw exception::QuditAlreadyMeasured("qpp::QCircuit::CTRL()");

            // check square matrix for the gate
            if (!internal::check_square_mat(U))
                throw exception::MatrixNotSquare("qpp::QCircuit::CTRL()");
            // check correct dimension
            if (static_cast<idx>(U.rows()) != d_)
                throw exception::DimsMismatchMatrix("qpp::QCircuit::CTRL()");

            // check shift
            if (shift >= d_)
                throw exception::OutOfRange("qpp::QCircuit::CTRL()");
        } catch (exception::Exception&) {
            std::cerr << "At STEP " << get_step_count() << "\n";
            throw;
        }
        // END EXCEPTION CHECKS

        if (name.empty()) {
            std::string gate_name = qpp::Gates::get_instance().get_name(U);
            name = gate_name.empty() ? "CTRL" : "CTRL-" + gate_name;
        }
        std::size_t hashU = hash_eigen(U);
        add_hash_(U, hashU);
        gates_.emplace_back(GateType::SINGLE_CTRL_SINGLE_TARGET, hashU,
                            std::vector<idx>{ctrl}, std::vector<idx>{target},
                            std::vector<idx>{shift}, name);
        step_types_.emplace_back(StepType::GATE);
        ++count_[name];

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
     * \param shift Performs the control as if the \a ctrl qudit state was
     * \f$ X\f$-incremented by \a shift
     * \param name Optional gate name
     * \return Reference to the current instance
     */
    QCircuit& CTRL(const cmat& U, idx ctrl, const std::vector<idx>& target,
                   idx shift = 0, std::string name = {}) {
        // EXCEPTION CHECKS

        try {
            // check valid ctrl
            if (ctrl >= nq_)
                throw exception::OutOfRange("qpp::QCircuit::CTRL()");
            if (get_measured(ctrl))
                throw exception::QuditAlreadyMeasured("qpp::QCircuit::CTRL()");

            // check valid target
            if (target.empty())
                throw exception::ZeroSize("qpp::QCircuit::CTRL()");
            for (auto&& elem : target) {
                if (elem >= nq_)
                    throw exception::OutOfRange("qpp::QCircuit::CTRL()");
                // check target was not measured before
                if (get_measured(elem))
                    throw exception::QuditAlreadyMeasured(
                        "qpp::QCircuit::CTRL()");
            }
            // check no duplicates target
            if (!internal::check_no_duplicates(target))
                throw exception::Duplicates("qpp::QCircuit::CTRL()");

            // check ctrl and target don't share common elements
            for (auto&& elem : target)
                if (elem == ctrl)
                    throw exception::OutOfRange("qpp::QCircuit::CTRL()");

            // check square matrix for the gate
            if (!internal::check_square_mat(U))
                throw exception::MatrixNotSquare("qpp::QCircuit::CTRL()");
            // check correct dimension
            if (static_cast<idx>(U.rows()) != d_)
                throw exception::DimsMismatchMatrix("qpp::QCircuit::CTRL()");

            // check shift
            if (shift >= d_)
                throw exception::OutOfRange("qpp::QCircuit::CTRL()");
        } catch (exception::Exception&) {
            std::cerr << "At STEP " << get_step_count() << "\n";
            throw;
        }
        // END EXCEPTION CHECKS

        if (name.empty()) {
            std::string gate_name = qpp::Gates::get_instance().get_name(U);
            name = gate_name.empty() ? "CTRL" : "CTRL-" + gate_name;
        }
        std::size_t hashU = hash_eigen(U);
        add_hash_(U, hashU);
        gates_.emplace_back(GateType::SINGLE_CTRL_MULTIPLE_TARGET, hashU,
                            std::vector<idx>{ctrl}, target,
                            std::vector<idx>{shift}, name);
        step_types_.emplace_back(StepType::GATE);
        ++count_[name];

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
     * \param shift Performs the control as if the \a ctrl qudit states were
     * \f$ X\f$-incremented component-wise by \a shift. If non-empty (default),
     * the size of \a shift must be the same as the size of \a ctrl.
     * \param name Optional gate name
     * \return Reference to the current instance
     */
    QCircuit& CTRL(const cmat& U, const std::vector<idx>& ctrl, idx target,
                   const std::vector<idx>& shift = {}, std::string name = {}) {
        // EXCEPTION CHECKS

        try {
            // check valid ctrl
            for (auto&& elem : ctrl) {
                if (elem >= nq_)
                    throw exception::OutOfRange("qpp::QCircuit::CTRL()");
                // check ctrl was not measured before
                if (get_measured(elem))
                    throw exception::QuditAlreadyMeasured(
                        "qpp::QCircuit::CTRL()");
            }
            // check no duplicates ctrl
            if (!internal::check_no_duplicates(ctrl))
                throw exception::Duplicates("qpp::QCircuit::CTRL()");

            // check valid target
            if (target >= nq_)
                throw exception::OutOfRange("qpp::QCircuit::CTRL()");
            if (get_measured(target))
                throw exception::QuditAlreadyMeasured("qpp::QCircuit::CTRL()");

            // check ctrl and target don't share common elements
            for (auto&& elem : ctrl)
                if (elem == target)
                    throw exception::OutOfRange("qpp::QCircuit::CTRL()");

            // check square matrix for the gate
            if (!internal::check_square_mat(U))
                throw exception::MatrixNotSquare("qpp::QCircuit::CTRL()");
            // check correct dimension
            if (static_cast<idx>(U.rows()) != d_)
                throw exception::DimsMismatchMatrix("qpp::QCircuit::CTRL()");

            // check shift
            if (!shift.empty() && (shift.size() != ctrl.size()))
                throw exception::SizeMismatch("qpp::QCircuit::CTRL()");
            if (!shift.empty())
                for (auto&& elem : shift)
                    if (elem >= d_)
                        throw exception::OutOfRange("qpp::QCircuit::CTRL()");
        } catch (exception::Exception&) {
            std::cerr << "At STEP " << get_step_count() << "\n";
            throw;
        }
        // END EXCEPTION CHECKS

        if (name.empty()) {
            std::string gate_name = qpp::Gates::get_instance().get_name(U);
            name = gate_name.empty() ? "CTRL" : "CTRL-" + gate_name;
        }
        std::size_t hashU = hash_eigen(U);
        add_hash_(U, hashU);
        gates_.emplace_back(GateType::MULTIPLE_CTRL_SINGLE_TARGET, hashU, ctrl,
                            std::vector<idx>{target}, shift, name);
        step_types_.emplace_back(StepType::GATE);
        ++count_[name];

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
     * \param shift Performs the control as if the \a ctrl qudit states were
     * \f$ X\f$-incremented component-wise by \a shift. If non-empty (default),
     * the size of \a shift must be the same as the size of \a ctrl.
     * \param name Optional gate name
     * \return Reference to the current instance
     */
    QCircuit& CTRL(const cmat& U, const std::vector<idx>& ctrl,
                   const std::vector<idx>& target,
                   const std::vector<idx>& shift = {}, std::string name = {}) {
        // EXCEPTION CHECKS

        try {
            // check valid ctrl
            for (auto&& elem : ctrl) {
                if (elem >= nq_)
                    throw exception::OutOfRange("qpp::QCircuit::CTRL()");
                // check ctrl was not measured before
                if (get_measured(elem))
                    throw exception::QuditAlreadyMeasured(
                        "qpp::QCircuit::CTRL()");
            }
            // check no duplicates ctrl
            if (!internal::check_no_duplicates(ctrl))
                throw exception::Duplicates("qpp::QCircuit::CTRL()");

            // check valid target
            if (target.empty())
                throw exception::ZeroSize("qpp::QCircuit::CTRL()");
            for (auto&& elem : target) {
                if (elem >= nq_)
                    throw exception::OutOfRange("qpp::QCircuit::CTRL()");
                // check target was not measured before
                if (get_measured(elem))
                    throw exception::QuditAlreadyMeasured(
                        "qpp::QCircuit::CTRL()");
            }
            // check no duplicates target
            if (!internal::check_no_duplicates(target))
                throw exception::Duplicates("qpp::QCircuit::CTRL()");

            // check ctrl and target don't share common elements
            for (auto&& elem_ctrl : ctrl)
                for (auto&& elem_target : target)
                    if (elem_ctrl == elem_target)
                        throw exception::OutOfRange("qpp::QCircuit::CTRL()");

            // check square matrix for the gate
            if (!internal::check_square_mat(U))
                throw exception::MatrixNotSquare("qpp::QCircuit::CTRL()");
            // check correct dimension
            if (static_cast<idx>(U.rows()) != d_)
                throw exception::DimsMismatchMatrix("qpp::QCircuit::CTRL()");

            // check shift
            if (!shift.empty() && (shift.size() != ctrl.size()))
                throw exception::SizeMismatch("qpp::QCircuit::CTRL()");
            if (!shift.empty())
                for (auto&& elem : shift)
                    if (elem >= d_)
                        throw exception::OutOfRange("qpp::QCircuit::CTRL()");
        } catch (exception::Exception&) {
            std::cerr << "At STEP " << get_step_count() << "\n";
            throw;
        }
        // END EXCEPTION CHECKS

        if (name.empty()) {
            std::string gate_name = qpp::Gates::get_instance().get_name(U);
            name = gate_name.empty() ? "CTRL" : "CTRL-" + gate_name;
        }
        std::size_t hashU = hash_eigen(U);
        add_hash_(U, hashU);
        gates_.emplace_back(GateType::MULTIPLE_CTRL_MULTIPLE_TARGET, hashU,
                            ctrl, std::vector<idx>{target}, shift, name);
        step_types_.emplace_back(StepType::GATE);
        ++count_[name];

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
     * \param shift Performs the control as if the \a ctrl qudit states were
     * \f$ X\f$-incremented component-wise by \a shift. If non-empty (default),
     * the size of \a shift must be the same as the size of \a ctrl.
     * \param name Optional gate name
     * \return Reference to the current instance
     */
    QCircuit& CTRL_custom(const cmat& U, const std::vector<idx>& ctrl,
                          const std::vector<idx>& target,
                          const std::vector<idx>& shift = {},
                          std::string name = {}) {
        // EXCEPTION CHECKS

        idx n = static_cast<idx>(target.size());
        idx D = static_cast<idx>(std::llround(std::pow(d_, n)));

        try {
            // check valid ctrl
            for (auto&& elem : ctrl) {
                if (elem >= nq_)
                    throw exception::OutOfRange("qpp::QCircuit::CTRL_custom()");
                // check ctrl was not measured before
                if (get_measured(elem))
                    throw exception::QuditAlreadyMeasured(
                        "qpp::QCircuit::CTRL_custom()");
            }
            // check no duplicates ctrl
            if (!internal::check_no_duplicates(ctrl))
                throw exception::Duplicates("qpp::QCircuit::CTRL_custom()");

            // check valid target
            if (target.empty())
                throw exception::ZeroSize("qpp::QCircuit::CTRL_custom()");
            for (auto&& elem : target) {
                if (elem >= nq_)
                    throw exception::OutOfRange("qpp::QCircuit::CTRL_custom()");
                // check target was not measured before
                if (get_measured(elem))
                    throw exception::QuditAlreadyMeasured(
                        "qpp::QCircuit::CTRL_custom()");
            }

            // check ctrl and target don't share common elements
            for (auto&& elem_ctrl : ctrl)
                for (auto&& elem_target : target)
                    if (elem_ctrl == elem_target)
                        throw exception::OutOfRange(
                            "qpp::QCircuit::CTRL_custom()");

            // check square matrix for the gate
            if (!internal::check_square_mat(U))
                throw exception::MatrixNotSquare(
                    "qpp::QCircuit::CTRL_custom()");
            // check correct dimension
            if (static_cast<idx>(U.rows()) != D)
                throw exception::DimsMismatchMatrix(
                    "qpp::QCircuit::CTRL_custom()");

            // check shift
            if (!shift.empty() && (shift.size() != ctrl.size()))
                throw exception::SizeMismatch("qpp::QCircuit::CTRL_custom()");
            if (!shift.empty())
                for (auto&& elem : shift)
                    if (elem >= d_)
                        throw exception::OutOfRange(
                            "qpp::QCircuit::CTRL_custom()");
        } catch (exception::Exception&) {
            std::cerr << "At STEP " << get_step_count() << "\n";
            throw;
        }
        // END EXCEPTION CHECKS

        if (name.empty()) {
            std::string gate_name = qpp::Gates::get_instance().get_name(U);
            name = gate_name.empty() ? "CTRL" : "CTRL-" + gate_name;
        }
        std::size_t hashU = hash_eigen(U);
        add_hash_(U, hashU);
        gates_.emplace_back(GateType::CUSTOM_CTRL, hashU, ctrl, target, shift,
                            name);
        step_types_.emplace_back(StepType::GATE);
        ++count_[name];

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
     * \param shift Performs the control as if the \a ctrl_dit classical dit was
     * incremented by \a shift
     * \param name Optional gate name
     * \return Reference to the current instance
     */
    QCircuit& cCTRL(const cmat& U, idx ctrl_dit, idx target, idx shift = 0,
                    std::string name = {}) {
        // EXCEPTION CHECKS

        try {
            // check valid ctrl_dit and target
            if (ctrl_dit >= nc_ || target >= nq_)
                throw exception::OutOfRange("qpp::QCircuit::cCTRL()");
            if (get_measured(target))
                throw exception::QuditAlreadyMeasured("qpp::QCircuit::cCTRL()");

            // check square matrix for the gate
            if (!internal::check_square_mat(U))
                throw exception::MatrixNotSquare("qpp::QCircuit::cCTRL()");
            // check correct dimension
            if (static_cast<idx>(U.rows()) != d_)
                throw exception::DimsMismatchMatrix("qpp::QCircuit::cCTRL()");

            // check shift
            if (shift >= d_)
                throw exception::OutOfRange("qpp::QCircuit::cCTRL()");
        } catch (exception::Exception&) {
            std::cerr << "At STEP " << get_step_count() << "\n";
            throw;
        }
        // END EXCEPTION CHECKS

        if (name.empty()) {
            std::string gate_name = qpp::Gates::get_instance().get_name(U);
            name = gate_name.empty() ? "cCTRL" : "cCTRL-" + gate_name;
        }

        std::size_t hashU = hash_eigen(U);
        add_hash_(U, hashU);
        gates_.emplace_back(GateType::SINGLE_cCTRL_SINGLE_TARGET, hashU,
                            std::vector<idx>{ctrl_dit},
                            std::vector<idx>{target}, std::vector<idx>{shift},
                            name);
        step_types_.emplace_back(StepType::GATE);
        ++count_[name];

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
     * \param shift Performs the control as if the \a ctrl_dit classical dit was
     * incremented by \a shift
     * \param name Optional gate name
     * \return Reference to the current instance
     */
    QCircuit& cCTRL(const cmat& U, idx ctrl_dit, const std::vector<idx>& target,
                    idx shift = 0, std::string name = {}) {
        // EXCEPTION CHECKS

        try {
            // check valid ctrl_dit
            if (ctrl_dit >= nc_)
                throw exception::OutOfRange("qpp::QCircuit::cCTRL()");

            // check valid target
            if (target.empty())
                throw exception::ZeroSize("qpp::QCircuit::cCTRL()");
            for (auto&& elem : target) {
                if (elem >= nq_)
                    throw exception::OutOfRange("qpp::QCircuit::cCTRL()");
                // check target was not measured before
                if (get_measured(elem))
                    throw exception::QuditAlreadyMeasured(
                        "qpp::QCircuit::cCTRL()");
            }
            // check no duplicates target
            if (!internal::check_no_duplicates(target))
                throw exception::Duplicates("qpp::QCircuit::cCTRL()");

            // check square matrix for the gate
            if (!internal::check_square_mat(U))
                throw exception::MatrixNotSquare("qpp::QCircuit::cCTRL()");
            // check correct dimension
            if (static_cast<idx>(U.rows()) != d_)
                throw exception::DimsMismatchMatrix("qpp::QCircuit::cCTRL()");

            // check shift
            if (shift >= d_)
                throw exception::OutOfRange("qpp::QCircuit::cCTRL()");
        } catch (exception::Exception&) {
            std::cerr << "At STEP " << get_step_count() << "\n";
            throw;
        }
        // END EXCEPTION CHECKS

        if (name.empty()) {
            std::string gate_name = qpp::Gates::get_instance().get_name(U);
            name = gate_name.empty() ? "cCTRL" : "cCTRL-" + gate_name;
        }
        std::size_t hashU = hash_eigen(U);
        add_hash_(U, hashU);
        gates_.emplace_back(GateType::SINGLE_cCTRL_MULTIPLE_TARGET, hashU,
                            std::vector<idx>{ctrl_dit}, target,
                            std::vector<idx>{shift}, name);
        step_types_.emplace_back(StepType::GATE);
        ++count_[name];

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
     * \param shift Performs the control as if the \a ctrl_dits classical dits
     * were incremented component-wise by \a shift. If non-empty (default), the
     * size of \a shift must be the same as the size of \a ctrl_dits.
     * \param name Optional gate name
     * \return Reference to the current instance
     */
    QCircuit& cCTRL(const cmat& U, const std::vector<idx>& ctrl_dits,
                    idx target, const std::vector<idx>& shift = {},
                    std::string name = {}) {
        // EXCEPTION CHECKS

        try {
            // check valid ctrl_dits
            for (auto&& elem : ctrl_dits) {
                if (elem >= nc_)
                    throw exception::OutOfRange("qpp::QCircuit::cCTRL()");
            }
            // check no duplicates ctrl_dits
            if (!internal::check_no_duplicates(ctrl_dits))
                throw exception::Duplicates("qpp::QCircuit::cCTRL()");

            // check valid target
            if (target >= nq_)
                throw exception::OutOfRange("qpp::QCircuit::cCTRL()");
            // check target was not measured before
            if (get_measured(target))
                throw exception::QuditAlreadyMeasured("qpp::QCircuit::cCTRL()");

            // check square matrix for the gate
            if (!internal::check_square_mat(U))
                throw exception::MatrixNotSquare("qpp::QCircuit::cCTRL()");
            // check correct dimension
            if (static_cast<idx>(U.rows()) != d_)
                throw exception::DimsMismatchMatrix("qpp::QCircuit::cCTRL()");

            // check shift
            if (!shift.empty() && (shift.size() != ctrl_dits.size()))
                throw exception::SizeMismatch("qpp::QCircuit::cCTRL()");
            if (!shift.empty())
                for (auto&& elem : shift)
                    if (elem >= d_)
                        throw exception::OutOfRange("qpp::QCircuit::cCTRL()");

        } catch (exception::Exception&) {
            std::cerr << "At STEP " << get_step_count() << "\n";
            throw;
        }
        // END EXCEPTION CHECKS

        if (name.empty()) {
            std::string gate_name = qpp::Gates::get_instance().get_name(U);
            name = gate_name.empty() ? "cCTRL" : "cCTRL-" + gate_name;
        }
        std::size_t hashU = hash_eigen(U);
        add_hash_(U, hashU);
        gates_.emplace_back(GateType::MULTIPLE_cCTRL_SINGLE_TARGET, hashU,
                            ctrl_dits, std::vector<idx>{target}, shift, name);
        step_types_.emplace_back(StepType::GATE);
        ++count_[name];

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
     * \param shift Performs the control as if the \a ctrl_dits classical dits
     * were incremented component-wise by \a shift. If non-empty (default), the
     * size of \a shift must be the same as the size of \a ctrl_dits.
     * \param name Optional gate name
     * \return Reference to the current instance
     */
    QCircuit& cCTRL(const cmat& U, const std::vector<idx>& ctrl_dits,
                    const std::vector<idx>& target,
                    const std::vector<idx>& shift = {}, std::string name = {}) {
        // EXCEPTION CHECKS

        try {
            // check valid ctrl_dits
            for (auto&& elem : ctrl_dits) {
                if (elem >= nc_)
                    throw exception::OutOfRange("qpp::QCircuit::cCTRL()");
            }
            // check no duplicates ctrl_dits
            if (!internal::check_no_duplicates(ctrl_dits))
                throw exception::Duplicates("qpp::QCircuit::cCTRL()");

            // check valid target
            if (target.empty())
                throw exception::ZeroSize("qpp::QCircuit::cCTRL()");
            for (auto&& elem : target) {
                if (elem >= nq_)
                    throw exception::OutOfRange("qpp::QCircuit::cCTRL()");
                // check target was not measured before
                if (get_measured(elem))
                    throw exception::QuditAlreadyMeasured(
                        "qpp::QCircuit::cCTRL()");
            }
            // check no duplicates target
            if (!internal::check_no_duplicates(target))
                throw exception::Duplicates("qpp::QCircuit::cCTRL()");

            // check square matrix for the gate
            if (!internal::check_square_mat(U))
                throw exception::MatrixNotSquare("qpp::QCircuit::cCTRL()");
            // check correct dimension
            if (static_cast<idx>(U.rows()) != d_)
                throw exception::DimsMismatchMatrix("qpp::QCircuit::cCTRL()");

            // check shift
            if (!shift.empty() && (shift.size() != ctrl_dits.size()))
                throw exception::SizeMismatch("qpp::QCircuit::cCTRL()");
            if (!shift.empty())
                for (auto&& elem : shift)
                    if (elem >= d_)
                        throw exception::OutOfRange("qpp::QCircuit::cCTRL()");
        } catch (exception::Exception&) {
            std::cerr << "At STEP " << get_step_count() << "\n";
            throw;
        }
        // END EXCEPTION CHECKS

        if (name.empty()) {
            std::string gate_name = qpp::Gates::get_instance().get_name(U);
            name = gate_name.empty() ? "cCTRL" : "cCTRL-" + gate_name;
        }
        std::size_t hashU = hash_eigen(U);
        add_hash_(U, hashU);
        gates_.emplace_back(GateType::MULTIPLE_cCTRL_MULTIPLE_TARGET, hashU,
                            ctrl_dits, std::vector<idx>{target}, shift, name);
        step_types_.emplace_back(StepType::GATE);
        ++count_[name];

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
     * \param shift Performs the control as if the \a ctrl_dits classical dits
     * were incremented component-wise by \a shift. If non-empty (default), the
     * size of \a shift must be the same as the size of \a ctrl_dits.
     * \param name Optional gate name
     * \return Reference to the current instance
     */
    QCircuit& cCTRL_custom(const cmat& U, const std::vector<idx>& ctrl_dits,
                           const std::vector<idx>& target,
                           const std::vector<idx>& shift = {},
                           std::string name = {}) {
        // EXCEPTION CHECKS

        idx n = static_cast<idx>(target.size());
        idx D = static_cast<idx>(std::llround(std::pow(d_, n)));

        try {
            // check valid ctrl_dits
            for (auto&& elem : ctrl_dits) {
                if (elem >= nc_)
                    throw exception::OutOfRange(
                        "qpp::QCircuit::cCTRL_custom()");
            }
            // check no duplicates ctrl_dits
            if (!internal::check_no_duplicates(ctrl_dits))
                throw exception::Duplicates("qpp::QCircuit::cCTRL_custom()");

            // check valid target
            if (target.empty())
                throw exception::ZeroSize("qpp::QCircuit::cCTRL_custom()");
            for (auto&& elem : target) {
                if (elem >= nq_)
                    throw exception::OutOfRange(
                        "qpp::QCircuit::cCTRL_custom()");
                // check target was not measured before
                if (get_measured(elem))
                    throw exception::QuditAlreadyMeasured(
                        "qpp::QCircuit::cCTRL_custom()");
            }
            // check no duplicates target
            if (!internal::check_no_duplicates(target))
                throw exception::Duplicates("qpp::QCircuit::cCTRL_custom()");

            // check square matrix for the gate
            if (!internal::check_square_mat(U))
                throw exception::MatrixNotSquare(
                    "qpp::QCircuit::cCTRL_custom()");
            // check correct dimension
            if (static_cast<idx>(U.rows()) != D)
                throw exception::DimsMismatchMatrix(
                    "qpp::QCircuit::cCTRL_custom()");

            // check shift
            if (!shift.empty() && (shift.size() != ctrl_dits.size()))
                throw exception::SizeMismatch("qpp::QCircuit::cCTRL()");
            if (!shift.empty())
                for (auto&& elem : shift)
                    if (elem >= d_)
                        throw exception::OutOfRange("qpp::QCircuit::cCTRL()");
        } catch (exception::Exception&) {
            std::cerr << "At STEP " << get_step_count() << "\n";
            throw;
        }
        // END EXCEPTION CHECKS

        if (name.empty()) {
            std::string gate_name = qpp::Gates::get_instance().get_name(U);
            name = gate_name.empty() ? "cCTRL" : "cCTRL-" + gate_name;
        }
        std::size_t hashU = hash_eigen(U);
        add_hash_(U, hashU);
        gates_.emplace_back(GateType::CUSTOM_cCTRL, hashU, ctrl_dits, target,
                            shift, name);
        step_types_.emplace_back(StepType::GATE);
        ++count_[name];

        return *this;
    }

    // Z measurement of single qudit
    /**
     * \brief Measurement of single qudit in the computational basis (Z-basis)
     *
     * \param target Target qudit index that is measured
     * \param c_reg Classical register where the value of the measurement is
     * being stored
     * \param destructive Destructive measurement, true by default
     * \param name Optional measurement name, default is "mZ"
     * \return Reference to the current instance
     */
    QCircuit& measureZ(idx target, idx c_reg, bool destructive = true,
                       std::string name = {}) {
        // EXCEPTION CHECKS

        try {
            // measuring non-existing qudit
            if (target >= nq_)
                throw exception::OutOfRange("qpp::QCircuit::measureZ()");
            // trying to put the result into a non-existing classical slot
            if (c_reg >= nc_)
                throw exception::OutOfRange("qpp::QCircuit::measureZ()");
            // qudit was measured before
            if (get_measured(target))
                throw exception::QuditAlreadyMeasured(
                    "qpp:QCircuit::measureZ()");
        } catch (exception::Exception&) {
            std::cerr << "At STEP " << get_step_count() << "\n";
            throw;
        }
        // END EXCEPTION CHECKS

        if (name.empty())
            name = "mZ";
        if (destructive) {
            measured_[target] = true;
            measurements_.emplace_back(MeasureType::MEASURE_Z,
                                       std::vector<std::size_t>{},
                                       std::vector<idx>{target}, c_reg, name);
        } else {
            measurements_.emplace_back(MeasureType::MEASURE_Z_ND,
                                       std::vector<std::size_t>{},
                                       std::vector<idx>{target}, c_reg, name);
        }
        step_types_.emplace_back(StepType::MEASUREMENT);
        ++measurement_count_[name];

        return *this;
    }

    // Z measurement of multiple qudits
    /**
     * \brief Measurement of multiple qudits in the computational basis
     * (Z-basis)
     *
     * \param target Target qudit indexes that are measured
     * \param c_reg Classical register where the value of the measurement is
     * being stored, as a decimal representation of the binary string
     * representing the measurement, with the most significant dit on the left
     * (corresponding to the first qudit that is being measured)
     * \param destructive Destructive measurement, true by default
     * \param name Optional measurement name, default is "mZ"
     * \return Reference to the current instance
     */
    QCircuit& measureZ(const std::vector<idx>& target, idx c_reg,
                       bool destructive = true, std::string name = {}) {
        // EXCEPTION CHECKS

        try {
            // check valid target
            if (target.empty())
                throw exception::ZeroSize("qpp::QCircuit::measureZ()");
            for (auto&& elem : target) {
                // measuring non-existing qudit
                if (elem >= nq_)
                    throw exception::OutOfRange("qpp::QCircuit::measureZ()");
                // qudit was measured before
                if (get_measured(elem))
                    throw exception::QuditAlreadyMeasured(
                        "qpp:QCircuit::measureZ()");
            }
            // trying to put the result into a non-existing classical slot
            if (c_reg >= nc_)
                throw exception::OutOfRange("qpp::QCircuit::measureZ()");
        } catch (exception::Exception&) {
            std::cerr << "At STEP " << get_step_count() << "\n";
            throw;
        }
        // END EXCEPTION CHECKS

        if (name.empty())
            name = "mZ";

        if (destructive) {
            for (auto&& elem : target) {
                measured_[elem] = true;
            }
            measurements_.emplace_back(MeasureType::MEASURE_Z_MANY,
                                       std::vector<std::size_t>{}, target,
                                       c_reg, name);
        } else {
            measurements_.emplace_back(MeasureType::MEASURE_Z_MANY_ND,
                                       std::vector<std::size_t>{}, target,
                                       c_reg, name);
        }
        step_types_.emplace_back(StepType::MEASUREMENT);
        ++measurement_count_[name];

        return *this;
    }

    // measurement of single qudit in the orthonormal basis or rank-1 projectors
    // specified by the columns of matrix V
    /**
     * \brief Measurement of single qudit in the orthonormal basis or rank-1
     * projectors specified by the columns of matrix \a V
     *
     * \param V Orthonormal basis or rank-1 projectors specified by the columns
     * of matrix V
     * \param target Target qudit index that is measured
     * \param c_reg Classical register where the value of the measurement is
     * stored
     * \param destructive Destructive measurement, true by default
     * \param name Optional measurement name
     * \return Reference to the current instance
     */
    QCircuit& measureV(const cmat& V, idx target, idx c_reg,
                       bool destructive = true, std::string name = {}) {
        // EXCEPTION CHECKS

        try {
            // measuring non-existing qudit
            if (target >= nq_)
                throw exception::OutOfRange("qpp::QCircuit::measureV()");
            // trying to put the result into a non-existing classical slot
            if (c_reg >= nc_)
                throw exception::OutOfRange("qpp::QCircuit::measureV()");
            // qudit was measured before
            if (get_measured(target))
                throw exception::QuditAlreadyMeasured(
                    "qpp:QCircuit::measureV()");
        } catch (exception::Exception&) {
            std::cerr << "At STEP " << get_step_count() << "\n";
            throw;
        }
        // END EXCEPTION CHECKS

        if (name.empty())
            name = "m" + qpp::Gates::get_instance().get_name(V);

        if (destructive) {
            measured_[target] = true;
            measurements_.emplace_back(MeasureType::MEASURE_V,
                                       std::vector<std::size_t>{hash_eigen(V)},
                                       std::vector<idx>{target}, c_reg, name);
        } else {
            measurements_.emplace_back(MeasureType::MEASURE_V_ND,
                                       std::vector<std::size_t>{hash_eigen(V)},
                                       std::vector<idx>{target}, c_reg, name);
        }
        step_types_.emplace_back(StepType::MEASUREMENT);
        ++measurement_count_[name];

        return *this;
    }

    // measurement of multiple qudits in the orthonormal basis or rank-1
    // projectors specified by the columns of matrix V
    /**
     * \brief Joint measurement of multiple qudits in the orthonormal basis or
     * rank-1 projectors specified by the columns of matrix \a V
     *
     * \param V Orthonormal basis or rank-1 projectors specified by the columns
     * of matrix V
     * \param target Target qudit indexes that are jointly measured
     * \param c_reg Classical register where the value of the measurement is
     * stored
     * \param destructive Destructive measurement, true by default
     * \param name Optional measurement name
     * \return Reference to the current instance
     */
    QCircuit& measureV(const cmat& V, const std::vector<idx>& target, idx c_reg,
                       bool destructive = true, std::string name = {}) {
        // EXCEPTION CHECKS

        try {
            // check valid target
            if (target.empty())
                throw exception::ZeroSize("qpp::QCircuit::measureV()");
            for (auto&& elem : target) {
                if (elem >= nq_)
                    throw exception::OutOfRange("qpp::QCircuit::measureV()");
                // check target was not measured before
                if (get_measured(elem))
                    throw exception::QuditAlreadyMeasured(
                        "qpp::QCircuit::measureV()");
            }
            // check no duplicates target
            if (!internal::check_no_duplicates(target))
                throw exception::Duplicates("qpp::QCircuit::measureV()");

            // trying to put the result into a non-existing classical slot
            if (c_reg >= nc_)
                throw exception::OutOfRange("qpp::QCircuit::measureV()");
            // qudit was measured before
            for (auto&& elem : target)
                if (get_measured(elem))
                    throw exception::QuditAlreadyMeasured(
                        "qpp::QCircuit::measureV()");
        } catch (exception::Exception&) {
            std::cerr << "At STEP " << get_step_count() << "\n";
            throw;
        }
        // END EXCEPTION CHECKS

        if (name.empty())
            name = "m" + qpp::Gates::get_instance().get_name(V);

        if (destructive) {
            for (auto&& elem : target)
                measured_[elem] = true;
            measurements_.emplace_back(MeasureType::MEASURE_V_MANY,
                                       std::vector<std::size_t>{hash_eigen(V)},
                                       target, c_reg, name);
        } else {
            measurements_.emplace_back(MeasureType::MEASURE_V_MANY_ND,
                                       std::vector<std::size_t>{hash_eigen(V)},
                                       target, c_reg, name);
        }
        step_types_.emplace_back(StepType::MEASUREMENT);
        ++measurement_count_[name];

        return *this;
    }

    /**
     * \brief Discard single qudit by measuring it destructively in the
     * computational basis (Z-basis) and discarding the measurement result
     *
     * \param target Target qudit index that is discarded
     * \param name Optional discard operation name, default is "discard"
     * \return Reference to the current instance
     */
    QCircuit& discard(idx target, std::string name = {}) {
        // EXCEPTION CHECKS

        try {
            // discarding non-existing qudit
            if (target >= nq_)
                throw exception::OutOfRange("qpp::QCircuit::discard()");
            // qudit was measured before
            if (get_measured(target))
                throw exception::QuditAlreadyMeasured(
                    "qpp:QCircuit::discard()");
        } catch (exception::Exception&) {
            std::cerr << "At STEP " << get_step_count() << "\n";
            throw;
        }
        // END EXCEPTION CHECKS

        if (name.empty())
            name = "discard";

        measured_[target] = true;
        measurements_.emplace_back(MeasureType::DISCARD,
                                   std::vector<std::size_t>{},
                                   std::vector<idx>{target}, -1, name);

        step_types_.emplace_back(StepType::MEASUREMENT);
        ++measurement_count_[name];

        return *this;
    }

    /**
     * \brief Discard multiple qudits by measuring them destructively in the
     * computational basis (Z-basis) and discarding the measurement result
     *
     * \param target Target qudit indexes that are discarded
     * \param name Optional discard operation name, default is "discard"
     * \return Reference to the current instance
     */
    QCircuit& discard(const std::vector<idx>& target, std::string name = {}) {
        // EXCEPTION CHECKS

        try {
            // check valid target
            if (target.empty())
                throw exception::ZeroSize("qpp::QCircuit::discard()");
            for (auto&& elem : target) {
                // discarding non-existing qudit
                if (elem >= nq_)
                    throw exception::OutOfRange("qpp::QCircuit::discard()");
                // qudit was measured before
                if (get_measured(elem))
                    throw exception::QuditAlreadyMeasured(
                        "qpp:QCircuit::discard()");
            }
        } catch (exception::Exception&) {
            std::cerr << "At STEP " << get_step_count() << "\n";
            throw;
        }
        // END EXCEPTION CHECKS

        if (name.empty())
            name = "discard";

        for (auto&& elem : target) {
            measured_[elem] = true;
        }
        measurements_.emplace_back(MeasureType::DISCARD_MANY,
                                   std::vector<std::size_t>{}, target, -1,
                                   name);
        step_types_.emplace_back(StepType::MEASUREMENT);
        ++measurement_count_[name];

        return *this;
    }

    /**
     * \brief No operation (no-op)
     *
     * \note If the underlying step is executed on a noisy engine, then noise
     * acts before it
     *
     * \return Reference to the current instance
     */
    QCircuit& nop() {
        step_types_.emplace_back(StepType::NOP);

        return *this;
    }

    // reset single qudit
    /**
     * \brief Resets single qudit by first measuring it non-destructively in the
     * computational basis and discarding the measurement result, followed by
     * shifting it back to the \f$|0\rangle\f$ state
     *
     * \param target Target qudit index that is reset
     * \param name Optional name, default is "reset"
     * \return Reference to the current instance
     */
    QCircuit& reset(idx target, std::string name = {}) {
        // EXCEPTION CHECKS

        try {
            // resetting non-existing qudit
            if (target >= nq_)
                throw exception::OutOfRange("qpp::QCircuit::reset()");
            // qudit was measured before
            if (get_measured(target))
                throw exception::QuditAlreadyMeasured("qpp:QCircuit::reset()");
        } catch (exception::Exception&) {
            std::cerr << "At STEP " << get_step_count() << "\n";
            throw;
        }
        // END EXCEPTION CHECKS

        if (name.empty())
            name = "reset";
        measurements_.emplace_back(MeasureType::RESET,
                                   std::vector<std::size_t>{},
                                   std::vector<idx>{target}, -1, name);
        step_types_.emplace_back(StepType::MEASUREMENT);
        ++measurement_count_[name];

        return *this;
    }

    // reset multiple qudits
    /**
     * \brief Resets multiple qudits by first measuring them non-destructively
     * in the computational basis and discarding the measurement results,
     * followed by shifting them back to the \f$|0\cdots 0\rangle\f$ state
     *
     * \param target Target qudit indexes that are reset
     * \param name Optional measurement name, default is "reset"
     * \return Reference to the current instance
     */
    QCircuit& reset(const std::vector<idx>& target, std::string name = {}) {
        // EXCEPTION CHECKS

        try {
            // check valid target
            if (target.empty())
                throw exception::ZeroSize("qpp::QCircuit::reset()");
            for (auto&& elem : target) {
                // resetting non-existing qudit
                if (elem >= nq_)
                    throw exception::OutOfRange("qpp::QCircuit::reset()");
                // qudit was measured before
                if (get_measured(elem))
                    throw exception::QuditAlreadyMeasured(
                        "qpp:QCircuit::reset()");
            }
        } catch (exception::Exception&) {
            std::cerr << "At STEP " << get_step_count() << "\n";
            throw;
        }
        // END EXCEPTION CHECKS

        if (name.empty())
            name = "reset";
        measurements_.emplace_back(MeasureType::RESET_MANY,
                                   std::vector<std::size_t>{}, target, -1,
                                   name);
        step_types_.emplace_back(StepType::MEASUREMENT);
        ++measurement_count_[name];

        return *this;
    }

    /**
     * \brief Replicates the circuit, in place
     * \note The circuit should not contain any measurements when invoking this
     * member function
     *
     * \param n Number of repetitions. If \a n == 1, returns the original
     * circuit.
     * \return Reference to the current instance
     */
    QCircuit& replicate(idx n) {
        // EXCEPTION CHECKS

        if (n == 0)
            throw exception::OutOfRange("qpp::QCircuit::replicate()");
        if (!get_measured().empty())
            throw exception::QuditAlreadyMeasured("qpp::QCircuit::replicate()");
        if (n == 1)
            return *this;
        // END EXCEPTION CHECKS

        auto gates_copy = gates_;
        idx gates_size = gates_.size();

        auto step_types_copy = step_types_;
        idx step_types_size = step_types_.size();

        gates_.resize(gates_.size() * n);
        step_types_.resize(step_types_size * n);

        for (idx i = 0; i < n - 1; ++i) {
            std::copy(std::begin(gates_copy), std::end(gates_copy),
                      std::next(std::begin(gates_), (i + 1) * gates_size));
            std::copy(
                std::begin(step_types_copy), std::end(step_types_copy),
                std::next(std::begin(step_types_), (i + 1) * step_types_size));
        }

        for (auto& elem : count_)
            elem.second *= n;

        return *this;
    }

    /**
     * \brief Appends a quantum circuit description to the current one
     *
     * \note If the qudit indexes of the added quantum circuit description do
     * not totally overlap with the indexes of the current quantum circuit
     * description, then the required number of additional qudits are
     * automatically added to the current quantum circuit description
     *
     * \param other Quantum circuit description
     * \param pos_qudit The index of the first qudit of \a other quantum circuit
     * description relative to the index of the first qudit of the current
     * quantum circuit description, with the rest following in order. If
     * negative or greater than the total number of qudits of the current
     * quantum circuit description, then the required number of additional
     * qudits are automatically added to the current quantum circuit
     * description.
     * \param pos_dit The first classical dit of \a other is inserted before the
     * \a pos_dit classical dit index of the current quantum circuit description
     * (in the classical dits array), the rest following in order. By default,
     * insertion is performed at the end.
     * \return Reference to the current instance
     */
    QCircuit& add_circuit(QCircuit other, bigint pos_qudit, idx pos_dit = -1) {
        // EXCEPTION CHECKS

        // check equal dimensions
        if (other.d_ != d_)
            throw exception::DimsNotEqual("qpp::QCircuit::add_circuit()");
        // check classical dits
        if (pos_dit == static_cast<idx>(-1))
            pos_dit = nc_;
        else if (pos_dit > nc_)
            throw exception::OutOfRange("qpp::QCircuit::add_circuit()");
        // check overlapping qudits (in the current instance) were not already
        // destructively measured
        if (pos_qudit >= 0 && static_cast<idx>(pos_qudit) < nq_) {
            for (idx i = 0;
                 i < std::min(static_cast<idx>(nq_ - pos_qudit), other.nq_);
                 ++i)
                if (get_measured(pos_qudit + i))
                    throw exception::QuditAlreadyMeasured(
                        "qpp::QCircuit::add_circuit()");
        }
        // END EXCEPTION CHECKS

        // STEP 0: add additional qudits (if needed)
        idx extra_qudits = 0;
        // add qudits before beginning
        if (pos_qudit < 0) {
            extra_qudits = std::abs(pos_qudit);
            add_qudit(extra_qudits, 0);
        } else {
            idx tmp = pos_qudit + other.nq_;
            if (tmp > nq_) {
                extra_qudits = tmp - nq_;
                add_qudit(extra_qudits);
            }

            // STEP 1: modify the copy of the to-be-added circuit

            // update gate indexes
            for (auto& gate : other.gates_) {
                switch (gate.gate_type_) {
                    // classically-controlled gates
                    case GateType::SINGLE_cCTRL_SINGLE_TARGET:
                    case GateType::SINGLE_cCTRL_MULTIPLE_TARGET:
                    case GateType::MULTIPLE_cCTRL_SINGLE_TARGET:
                    case GateType::MULTIPLE_cCTRL_MULTIPLE_TARGET:
                    case GateType::CUSTOM_cCTRL:
                        for (auto& pos : gate.ctrl_) {
                            pos += pos_dit;
                        }
                        break;
                    // quantumly-controlled gates
                    default:
                        for (auto& pos : gate.ctrl_) {
                            pos += pos_qudit;
                        }
                        break;
                }

                // non-controlled gates
                for (auto& pos : gate.target_) {
                    pos += pos_qudit;
                }
            }

            // update measurement indexes
            for (auto& measurement : other.measurements_) {
                for (auto& pos : measurement.target_) {
                    pos += pos_qudit;
                }
                measurement.c_reg_ += pos_dit;
            }
        } // end else

        // STEP 2: append the copy of the to-be-added circuit to the current
        // instance
        // insert classical dits from the to-be-added circuit
        add_dit(other.nc_, pos_dit);

        // insert the measured vector
        measured_.insert(std::next(std::begin(measured_), pos_dit),
                         std::begin(other.measured_),
                         std::end(other.measured_));

        // append gate steps vector
        gates_.insert(std::end(gates_), std::begin(other.gates_),
                      std::end(other.gates_));

        // append measurement steps vector
        measurements_.insert(std::end(measurements_),
                             std::begin(other.measurements_),
                             std::end(other.measurements_));

        // append step types vector
        step_types_.insert(std::end(step_types_), std::begin(other.step_types_),
                           std::end(other.step_types_));

        // STEP 3: modify gate counts, hash tables etc accordingly
        // update matrix hash table
        for (auto& elem : other.cmat_hash_tbl_)
            cmat_hash_tbl_[elem.first] = elem.second;
        // update gate counts
        for (auto& elem : other.count_)
            count_[elem.first] += elem.second;
        // update measurement counts
        for (auto& elem : other.measurement_count_)
            measurement_count_[elem.first] += elem.second;

        return *this;
    }

    /**
     * \brief Kronecker product with another quantum circuit description, in
     * place
     *
     * \param qc Quantum circuit description
     * \return Reference to the current instance
     */
    QCircuit& kron(const QCircuit& qc) {
        add_circuit(qc, nq_);

        return *this;
    }

    /**
     * \brief Adjoint quantum circuit description, in place
     *
     * \return Reference to the current instance
     */
    QCircuit& adjoint() {
        // EXCEPTION CHECKS

        if (get_measured().size() > 0)
            throw exception::QuditAlreadyMeasured("qpp::QCircuit()::adjoint()");
        // END EXCEPTION CHECKS

        auto htbl = cmat_hash_tbl_; // copy the gate hash table of other
        cmat_hash_tbl_.clear();

        std::reverse(std::begin(gates_), std::end(gates_));
        std::reverse(std::begin(step_types_), std::end(step_types_));

        for (auto& elem : gates_) {
            // get the gate and its corresponding hash
            std::size_t hashU = elem.gate_hash_;
            cmat U = htbl[hashU];

            // compute the adjoints
            cmat Udagger = qpp::adjoint(U);
            std::size_t hashUdagger = hash_eigen(Udagger);

            // modify and add hash
            elem.gate_hash_ = hashUdagger;
            if (!elem.name_.empty())
                elem.name_ += "+";
            add_hash_(Udagger, hashUdagger);
        }

        return *this;
    }

    /**
     * \brief qpp::IJSON::to_JSON() override
     *
     * Displays the quantum circuit description in JSON format
     *
     * \param enclosed_in_curly_brackets If true, encloses the result in curly
     * brackets
     * \return String containing the JSON representation of
     * the quantum circuit description
     */
    std::string to_JSON(bool enclosed_in_curly_brackets = true) const override {
        std::string result;

        if (enclosed_in_curly_brackets)
            result += "{";

        result += "\"nq\" : " + std::to_string(nq_);
        result += ", \"nc\" : " + std::to_string(nc_);
        result += ", \"d\" : " + std::to_string(d_);
        result += ", \"name\" : \"" + name_ + '\"';

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
            // gate step
            if (elem.type_ == StepType::GATE) {
                idx pos = std::distance(std::begin(elem.value_type_qc_->gates_),
                                        elem.gates_ip_);
                ss.str("");
                ss.clear();
                ss << gates_[pos].gate_type_;
                result += '\"' + ss.str() + "\", ";
                if (!gates_[pos].ctrl_.empty()) {
                    ss.str("");
                    ss.clear();
                    ss << disp(gates_[pos].ctrl_, ", ");
                    if (gates_[pos].gate_type_ >=
                        GateType::SINGLE_cCTRL_SINGLE_TARGET)
                        result += "\"c_ctrl\" : " + ss.str() + ", ";
                    else
                        result += "\"ctrl\" : " + ss.str() + ", ";
                }
                ss.str("");
                ss.clear();
                ss << disp(gates_[pos].target_, ", ");
                result += "\"target\" : " + ss.str() + ", ";

                if (!gates_[pos].shift_.empty()) {
                    ss.str("");
                    ss.clear();
                    ss << disp(gates_[pos].shift_, ", ");
                    result += "\"shift\" : " + ss.str() + ", ";
                }

                result += "\"name\" : ";
                result += '\"' + gates_[pos].name_ + "\"}";
            }
            // measurement step
            else if (elem.type_ == StepType::MEASUREMENT) {
                idx pos = std::distance(
                    std::begin(elem.value_type_qc_->measurements_),
                    elem.measurements_ip_);
                ss.str("");
                ss.clear();
                ss << measurements_[pos].measurement_type_;
                result += '\"' + ss.str() + "\", ";
                ss.str("");
                ss.clear();
                ss << disp(measurements_[pos].target_, ", ");
                result += "\"target\" : " + ss.str() + ", ";

                if (measurements_[pos].measurement_type_ !=
                        MeasureType::RESET &&
                    measurements_[pos].measurement_type_ !=
                        MeasureType::RESET_MANY &&
                    measurements_[pos].measurement_type_ !=
                        MeasureType::DISCARD &&
                    measurements_[pos].measurement_type_ !=
                        MeasureType::DISCARD_MANY)
                    result += "\"c_reg\" : " +
                              std::to_string(measurements_[pos].c_reg_) + ", ";

                result += "\"name\" : ";
                result += '\"' + measurements_[pos].name_ + "\"}";

            }
            // no-op
            else if (elem.type_ == StepType::NOP) {
                result += std::string{"\"NOP\""} + "}";
            }
            // otherwise
            else {
            }
        }                // end for
        result += "], "; // end steps

        result += "\"step count\" : " + std::to_string(get_step_count()) + ", ";
        result +=
            "\"total gate count\" : " + std::to_string(get_gate_count()) + ", ";
        result +=
            "\"total gate depth\" : " + std::to_string(get_gate_depth()) + ", ";
        result += "\"total measurement count\" : " +
                  std::to_string(get_measurement_count()) + ", ";
        result += "\"total measurement depth\" : " +
                  std::to_string(get_measurement_depth()) + ", ";
        result += "\"total depth\" : " + std::to_string(get_depth()) + ", ";

        ss.str("");
        ss.clear();
        ss << disp(get_measured(), ", ");
        result += "\"measured/discarded (destructive)\" : " + ss.str() + ", ";

        ss.str("");
        ss.clear();
        ss << disp(get_non_measured(), ", ");
        result += "\"non-measured/non-discarded\" : " + ss.str();

        if (enclosed_in_curly_brackets)
            result += "}";

        return result;
    } /* to_JSON() */

    /**
     * \brief Adjoint quantum circuit description
     *
     * \param qc Quantum circuit description
     * \return Adjoint quantum circuit description
     */
    friend QCircuit adjoint(QCircuit qc) {
        // EXCEPTION CHECKS

        if (qc.get_measured().size() > 0)
            throw exception::QuditAlreadyMeasured("qpp::adjoint()");
        // END EXCEPTION CHECKS

        return qc.adjoint();
    }

    /**
     * \brief Kronecker product between two quantum circuit descriptions
     *
     * \param qc1 Quantum circuit description
     * \param qc2 Quantum circuit description
     * \return Quantum circuit description of the Kronecker product of \a qc1
     * with \a qc2
     */
    friend QCircuit kron(const QCircuit& qc1, const QCircuit& qc2) {
        QCircuit qc{qc1}; // create a copy

        return qc.kron(qc2);
    }

    /**
     * \brief Replicates the circuit
     * \note The circuit should not contain any measurements when invoking this
     * function
     *
     * \param qc Quantum circuit description
     * \param n Number of repetitions. If \a n == 1, returns the original
     * circuit.
     * \return Reference to the current instance
     */
    friend QCircuit replicate(QCircuit qc, idx n) {
        // EXCEPTION CHECKS

        if (n == 0)
            throw exception::OutOfRange("qpp::replicate()");
        if (!qc.get_measured().empty())
            throw exception::QuditAlreadyMeasured("qpp::replicate()");
        if (n == 1)
            return qc;
        // END EXCEPTION CHECKS

        return qc.replicate(n);
    }

    /**
     * \brief Appends a quantum circuit description to another one
     *
     * \note If qudit indexes of the second quantum circuit description do not
     * totally overlap with the indexes of the first quantum circuit
     * description, then the required number of additional qudits are
     * automatically added to the output quantum circuit description
     *
     * \param qc1 Quantum circuit description
     * \param qc2 Quantum circuit description
     * \param pos_qudit The index of the first qudit of \a qc2 quantum circuit
     * description relative to the index of the first qudit of the \a qc1
     * quantum circuit description, with the rest following in order. If
     * negative or greater than the total number of qudits of \a qc1, then the
     * required number of additional qudits are automatically added to the
     * output quantum circuit description.
     * \param pos_dit The first classical dit of \a qc2 quantum circuit
     * description is inserted before the \a pos_dit classical dit index of
     * \a qc1 quantum circuit description (in the classical dits array), the
     * rest following in order. By default, insertion is performed at the end.
     * \return Combined quantum circuit description
     */
    friend QCircuit add_circuit(QCircuit qc1, const QCircuit& qc2,
                                bigint pos_qudit, idx pos_dit = -1) {
        return qc1.add_circuit(qc2, pos_qudit, pos_dit);
    }

  private:
    /**
     * \brief qpp::IDisplay::display() override
     *
     * Writes to the output stream a textual representation of the quantum
     * circuit
     *
     * \param os Output stream passed by reference
     * \return Reference to the output stream
     */
    std::ostream& display(std::ostream& os) const override {
        os << "nq = " << nq_ << ", nc = " << nc_ << ", d = " << d_;

        if (!name_.empty()) // if the circuit is named
            os << ", name = \"" << name_ << "\"\n";
        else
            os << ", name = \"\"\n";

        for (auto&& elem : *this) {
            os << elem << '\n';
        }

        os << "step count: " << get_step_count() << '\n';
        os << "total gate count: " << get_gate_count() << '\n';
        os << "total gate depth: " << get_gate_depth() << '\n';
        os << "total measurement count: " << get_measurement_count() << '\n';
        os << "total measurement depth: " << get_measurement_depth() << '\n';
        os << "total depth: " << get_depth() << '\n';
        os << "measured/discarded (destructive): " << disp(get_measured(), ", ")
           << '\n';
        os << "non-measured/non-discarded: " << disp(get_non_measured(), ", ");

        return os;
    }
}; /* class QCircuit */

} /* namespace qpp */

#endif /* CLASSES_CIRCUITS_CIRCUITS_H_ */
