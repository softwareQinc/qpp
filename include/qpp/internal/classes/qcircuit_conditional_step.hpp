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
 * \file qpp/internal/classes/qcircuit_conditional_step.hpp
 * \brief qpp::internal::QCircuitConditionalStep
 */

#ifndef QPP_INTERNAL_CLASSES_QCIRCUIT_CONDITIONAL_STEP_HPP_
#define QPP_INTERNAL_CLASSES_QCIRCUIT_CONDITIONAL_STEP_HPP_

#include <optional>

#include "qpp/input_output.hpp"

#include "qpp/classes/idisplay.hpp"

namespace qpp {
namespace internal {
/**
 * \brief One step consisting only of measurements in the circuit
 */
struct QCircuitConditionalStep : IDisplay {
    /**
     * \brief Type of measurement being executed in a measurement step
     */
    enum class Type {
        NONE, ///< no condition

        IF, ///< if branch statement

        ELSE, ///< else branch statement

        ENDIF, ///< end if statement
    };

    /**
     * \brief Conditional functor type in qpp::QCircuit conditional statements
     */
    using cond_func_t = std::function<bool(std::vector<idx>)>;

    /**
     * \class qpp::internal::QCircuitConditionalStep::IfElseLoc
     * \brief Stores the location of conditional statements
     */
    struct Context {
        std::optional<std::pair<idx, cond_func_t>>
            if_expr; ///< location of if statement and corresponding condition
                     ///< function
        std::optional<idx> else_expr;  ///< location of else statement
        std::optional<idx> endif_expr; ///< location of endif statement
    };

    /**
     * \brief Extraction operator overload for
     * qpp::internal::QCircuitConditionalStep::Type enum class
     *
     * \param os Output stream passed by reference
     * \param measure_type qpp::internal::QCircuitConditionalStep::Type enum
     * class
     * \return Reference to the output stream
     */
    friend std::ostream& operator<<(std::ostream& os,
                                    const Type& condition_type) {
        switch (condition_type) {
            case Type::NONE:
                os << "CONDITIONAL NONE";
                break;
            case Type::IF:
                os << "IF";
                break;
            case Type::ELSE:
                os << "ELSE";
                break;
            case Type::ENDIF:
                os << "ENDIF";
                break;
        }

        return os;
    }

    Type condition_type_ = Type::NONE; ///< condition type
    Context ctx_;

    /**
     * \brief Default constructor
     */
    QCircuitConditionalStep() = default;

    /**
     * \brief Constructs a conditional step instance
     *
     * \tparam Func Conditional functor
     * \param condition_type Measurement type
     */
    explicit QCircuitConditionalStep(Type condition_type, Context ctx)
        : condition_type_{condition_type}, ctx_(std::move(ctx)) {}

    /**
     * \brief Equality operator
     *
     * \param rhs qpp::internal::QCircuitConditionalStep against which the
     * equality is being tested
     * \return True if the qpp::internal::QCircuitConditionalStep(s) are equal,
     * false otherwise
     */
    // FIXME:
    bool operator==(const QCircuitConditionalStep& rhs) const noexcept {
        return std::tie(rhs.condition_type_) == std::tie(condition_type_);
    }

    /**
     * \brief Inequality operator
     *
     * \param rhs qpp::internal::QCircuitConditionalStep against which the
     * inequality is being tested
     * \return True if the qpp::internal::QCircuitConditionalStep(s) are not
     * equal, false otherwise
     */
    bool operator!=(const QCircuitConditionalStep& rhs) const noexcept {
        return !(*this == rhs);
    }

  private:
    /**
     * \brief qpp::IDisplay::display() override
     *
     * Writes to the output stream a textual representation of the
     * \a qpp::internal::QCircuitConditionalStep instance
     *
     * \param os Output stream passed by reference
     * \return Reference to the output stream
     */
    std::ostream& display(std::ostream& os) const override {
        os << condition_type_;
        if (condition_type_ == Type::IF) {
            os << ", <FUNC " << std::addressof(ctx_.if_expr->second) << ">";
        }

        return os;
    }
};

} /* namespace internal */
} /* namespace qpp */

#endif /* QPP_INTERNAL_CLASSES_QCIRCUIT_CONDITIONAL_STEP_HPP_ */
