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
 * \file qpp/internal/classes/qcircuit_gate_step.hpp
 * \brief qpp::internal::QCircuitGateStep
 */

#ifndef QPP_INTERNAL_CLASSES_QCIRCUIT_GATE_STEP_HPP_
#define QPP_INTERNAL_CLASSES_QCIRCUIT_GATE_STEP_HPP_

#include "qpp/input_output.hpp"

#include "qpp/classes/idisplay.hpp"

namespace qpp {
namespace internal {
/**
 * \brief One step consisting only of gates/operators in the circuit
 */
struct QCircuitGateStep : IDisplay {
    /**
     * \brief Type of gate being executed in a gate step
     */
    enum class Type {
        NONE, ///< represents no gate

        SINGLE, ///< unitary gate on a single qudit

        TWO, ///< unitary gate on 2 qudits

        THREE, ///< unitary gate on 3 qudits

        JOINT, ///< joint gate on multiple qudits

        FAN, ///< the same unitary gate on multiple qudits

        CTRL, ///< controlled unitary gate with joint target

        CTRL_FAN, ///< controlled unitary gate with multiple targets

        cCTRL, ///< classically controlled unitary gate with joint target

        cCTRL_FAN, ///< classically controlled unitary gate with multiple
                   ///< targets
    };
    /**
     * \brief Extraction operator overload for
     * qpp::internal::QCircuitGateStep::Type enum class
     *
     * \param os Output stream passed by reference
     * \param gate_type qpp::internal::QCircuitGateStep::Type enum class
     * \return Reference to the output stream
     */
    friend std::ostream& operator<<(std::ostream& os, const Type& gate_type) {
        switch (gate_type) {
            case Type::NONE:
                os << "GATE NONE";
                break;
            case Type::SINGLE:
                os << "SINGLE";
                break;
            case Type::TWO:
                os << "TWO";
                break;
            case Type::THREE:
                os << "THREE";
                break;
            case Type::FAN:
                os << "FAN";
                break;
            case Type::JOINT:
                os << "JOINT";
                break;
            case Type::CTRL:
                os << "CTRL";
                break;
            case Type::CTRL_FAN:
                os << "CTRL_FAN";
                break;
            case Type::cCTRL:
                os << "cCTRL";
                break;
            case Type::cCTRL_FAN:
                os << "cCTRL_FAN";
                break;
        }

        return os;
    }

    Type gate_type_ = Type::NONE;            ///< gate type
    std::size_t gate_hash_{};                ///< gate hash
    std::optional<std::vector<idx>> ctrl_{}; ///< control
    std::vector<idx> target_{}; ///< target where the gate is applied
    std::optional<std::vector<idx>> shift_{}; ///< shifts in CTRL gates
    std::optional<std::string> name_{};       ///< custom name of the gate(s)

    /**
     * \brief Default constructor
     */
    QCircuitGateStep() = default;

    /**
     * \brief Constructs a gate step instance
     *
     * \param gate_type Gate type
     * \param gate_hash Hash of the quantum gate
     * \param ctrl Optional control qudit indexes
     * \param target Target qudit indexes
     * \param shift Optional gate shifts (for CTRL gates)
     * \param name Optional gate name
     */
    explicit QCircuitGateStep(
        Type gate_type, std::size_t gate_hash,
        std::optional<std::vector<idx>> ctrl, std::vector<idx> target,
        std::optional<std::vector<idx>> shift = std::nullopt,
        std::optional<std::string> name = std::nullopt)
        : gate_type_{gate_type}, gate_hash_{gate_hash}, ctrl_{std::move(ctrl)},
          target_{std::move(target)}, shift_{std::move(shift)},
          name_{std::move(name)} {}

    /**
     * \brief Equality operator
     * \note Ignores gate names
     *
     * \param rhs qpp::internal::QCircuitGateStep against which the equality is
     * being tested
     * \return True if the qpp::internal::QCircuitGateStep(s) are equal, false
     * otherwise
     */
    bool operator==(const QCircuitGateStep& rhs) const noexcept {
        return std::tie(rhs.target_, rhs.shift_, rhs.gate_type_, rhs.gate_hash_,
                        rhs.ctrl_) ==
               std::tie(target_, shift_, gate_type_, gate_hash_, ctrl_);
    }

    /**
     * \brief Inequality operator
     * \note Ignores gate names
     *
     * \param rhs qpp::internal::QCircuitGateStep against which the inequality
     * is being tested
     * \return True if the qpp::internal::QCircuitGateStep(s) are not equal,
     * false otherwise
     */
    bool operator!=(const QCircuitGateStep& rhs) const noexcept {
        return !(*this == rhs);
    }

  private:
    /**
     * \brief qpp::IDisplay::display() override
     *
     * Writes to the output stream a textual representation of the
     * \a qpp::internal::QCircuitGateStep instance
     *
     * \param os Output stream passed by reference
     * \return Reference to the output stream
     */
    std::ostream& display(std::ostream& os) const override {
        os << gate_type_ << ", ";
        if (ctrl_.has_value()) {
            if (gate_type_ == Type::CTRL || gate_type_ == Type::CTRL_FAN) {
                os << "ctrl = ";
                os << disp(ctrl_.value(), IOManipContainerOpts{}.set_sep(", "))
                   << ", ";
            } else if (gate_type_ == Type::cCTRL ||
                       gate_type_ == Type::cCTRL_FAN) {
                os << "c_ctrl = "
                   << disp(ctrl_.value(), IOManipContainerOpts{}.set_sep(", "))
                   << ", ";
            }
        }
        os << "target = "
           << disp(target_, IOManipContainerOpts{}.set_sep(", "));
        // os << "shift = ";
        if (shift_.has_value()) {
            os << ", shift = ";
            os << disp(shift_.value(), IOManipContainerOpts{}.set_sep(", "));
        }
        if (name_.has_value()) {
            os << ", name = " << std::quoted(name_.value());
        }

        return os;
    }
};

} /* namespace internal */
}; /* namespace qpp */

#endif /* QPP_INTERNAL_CLASSES_QCIRCUIT_GATE_STEP_HPP_ */
