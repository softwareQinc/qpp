/*
 * This file is part of Quantum++.
 *
 * Copyright (c) 2017 - 2025 softwareQ Inc. All rights reserved.
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
 * \file qpp/internal/classes/qcircuit_runtime_step.hpp
 * \brief qpp::internal::QCircuitRuntimeStep
 */

#ifndef QPP_INTERNAL_CLASSES_QCIRCUIT_RUNTIME_STEP_HPP_
#define QPP_INTERNAL_CLASSES_QCIRCUIT_RUNTIME_STEP_HPP_

#include <functional>
#include <iterator>
#include <optional>
#include <ostream>
#include <tuple>
#include <utility>
#include <vector>

#include "qpp/classes/idisplay.hpp"
#include "qpp/internal/classes/labelled_vector_proxy.hpp"

namespace qpp {
/**
 * \brief Conditional functor type (boolean predicate) in qpp::QCircuit
 * conditional (runtime) statements
 */
using cond_pred_t = std::function<bool(internal::const_proxy_to_engine_dits_t)>;
// using cond_pred_t = std::function<bool(std::vector<idx>&)>;

/**
 * \brief Functor type in qpp::QCircuit used for overwriting quantum engine
 * dits at runtime
 */
using mutable_dits_functor_t =
    std::function<void(internal::proxy_to_engine_dits_t)>;

namespace internal {
/**
 * \brief One step consisting only of runtime statements such as
 * if/while/set_dits_runtime etc.
 */
struct QCircuitRuntimeStep : IDisplay {
    /**
     * \brief Type of runtime statement being executed in a runtime step
     */
    enum class Type {
        NONE, ///< no statement

        IF, ///< if branch statement

        ELSE, ///< else branch statement

        ENDIF, ///< end if block statement

        WHILE, ///< while branch statement

        ENDWHILE, ///< end while branch statement

        SET_DITS_RUNTIME, ///< overwrite quantum engine dits at runtime
    };

    /**
     * \brief Extraction operator overload for
     * qpp::internal::QCircuitRuntimeStep::Type enum class
     *
     * \param os Output stream passed by reference
     * \param runtime_type qpp::internal::QCircuitRuntimeStep::Type enum
     * class
     * \return Reference to the output stream
     */
    friend std::ostream& operator<<(std::ostream& os,
                                    const Type& runtime_type) {
        switch (runtime_type) {
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
            case Type::WHILE:
                os << "WHILE";
                break;
            case Type::ENDWHILE:
                os << "ENDWHILE";
                break;
            case Type::SET_DITS_RUNTIME:
                os << "SET_DITS_RUNTIME";
                break;
        }

        return os;
    }

    /**
     * \class qpp::internal::QCircuitRuntimeStep::Context
     * \brief Keeps track of the location of runtime statements
     */
    struct Context : IDisplay {
        /**
         * \brief Vector of pairs (offset, length) that keep track of where
         * classical dits were added by QCircuit::add_dit() or when composing
         * quantum circuit descriptions
         */
        std::vector<std::pair<idx, idx>> dit_ctx{};
        /**
         * \brief location of if/while statement and corresponding condition
         * functor (boolean predicate)
         */
        std::optional<std::pair<idx, cond_pred_t>> start_expr{};
        std::optional<idx> else_expr{}; ///< location of else statement
        std::optional<idx> end_expr{}; ///< location of endif/endwhile statement
        std::optional<mutable_dits_functor_t>
            set_dits_func{}; ///< functor that overwrites engine dits at runtime

        /**
         * \brief Given a vector of pairs (offset, length) that keep tracks of
         * where additional classical dits were added, restore the vector of
         * dits to its original form
         *
         * \param dits Vector of classical dits
         * \return Vector of classical dits, with the additional dits listed in
         * \a dit_ctx removed
         */
        std::vector<idx>
        restore_dits_from_dit_ctx(std::vector<idx> dits) const {
            if (dit_ctx.empty()) {
                return dits;
            }
            for (auto rit = dit_ctx.rbegin(); rit != dit_ctx.rend(); ++rit) {
                auto elem = *rit;
                auto start = std::next(dits.begin(), elem.first);
                auto finish = std::next(dits.begin(), elem.first + elem.second);
                dits.erase(start, finish);
            }

            return dits;
        }

        /**
         * \brief Increments each offset in the \a dit_ctx by \a i
         *
         * \param i Non-negative integer
         */
        void inc_offset_dit_ctx(idx i) {
            for (auto& elem : dit_ctx) {
                elem.first += i;
            }
        }

        /**
         * \brief Increments conditional locations by \a i
         *
         * \param i Non-negative integer
         */
        void inc_cond_locs(idx i) {
            if (start_expr.has_value()) {
                start_expr.value().first += i;
            }
            if (else_expr.has_value()) {
                else_expr.value() += i;
            }
            if (end_expr.has_value()) {
                end_expr.value() += i;
            }
        }

        /**
         * \brief Equality operator
         * \note Does not test equality of the comparison functions
         *
         * \param rhs Context against which the equality is being tested
         * \return True if the contexts are equal, false otherwise
         */
        bool operator==(const Context& rhs) const {
            if (dit_ctx != rhs.dit_ctx) {
                return false;
            }
            if (else_expr != rhs.else_expr) {
                return false;
            }
            if (end_expr != rhs.end_expr) {
                return false;
            }
            if (start_expr.has_value() && !rhs.start_expr.has_value()) {
                return false;
            }
            if (!start_expr.has_value() && rhs.start_expr.has_value()) {
                return false;
            }
            if (!start_expr.has_value() && !rhs.start_expr.has_value()) {
                return true;
            }

            auto start_expr_val = start_expr.value();
            auto rhs_start_expr_val = rhs.start_expr.value();
            return start_expr_val.first == rhs_start_expr_val.first;
        }

        /**
         * \brief qpp::IDisplay::display() override
         *
         * \param os Output stream passed by reference
         * \return Reference to the output stream
         */
        std::ostream& display(std::ostream& os) const override {
            if (start_expr.has_value()) {
                os << "START: " << start_expr.value().first;
                // os << " (addr. " << std::addressof(start_expr.value().second)
                //    << ")";
            }
            if (else_expr.has_value()) {
                os << ", ELSE: " << else_expr.value();
            }
            if (end_expr.has_value()) {
                os << ", END: " << end_expr.value();
            }

            // NOTE: comment in production
            /*
            if (!dit_ctx.empty()) {
                os << ", dit_ctx: [";
                bool first = true;
                for (auto&& elem : dit_ctx) {
                    if (first) {
                        first = false;
                    } else {
                        os << ", ";
                    }
                    os << '(' << elem.first << ',' << elem.second << ')';
                }
                os << ']';
            }
            */

            return os;
        }
    }; /* struct QCircuitRuntimeStep::Context */

    Type runtime_type_ = Type::NONE; ///< runtime type
    Context ctx_{};                  ///< additional context

    /**
     * \brief Default constructor
     */
    QCircuitRuntimeStep() = default;

    /**
     * \brief Constructs a runtime step instance
     *
     * \param runtime_type Runtime type
     * \param ctx Context
     */
    explicit QCircuitRuntimeStep(Type runtime_type, Context ctx)
        : runtime_type_{runtime_type}, ctx_(std::move(ctx)) {}

    /**
     * \brief Equality operator
     *
     * \param rhs qpp::internal::QCircuitRuntimeStep against which the
     * equality is being tested
     * \return True if the qpp::internal::QCircuitRuntimeStep(s) are
     * equal, false otherwise
     */
    bool operator==(const QCircuitRuntimeStep& rhs) const noexcept {
        return std::tie(rhs.runtime_type_, rhs.ctx_) ==
               std::tie(runtime_type_, ctx_);
    }

    /**
     * \brief Inequality operator
     *
     * \param rhs qpp::internal::QCircuitRuntimeStep against which the
     * inequality is being tested
     * \return True if the qpp::internal::QCircuitRuntimeStep(s) are not
     * equal, false otherwise
     */
    bool operator!=(const QCircuitRuntimeStep& rhs) const noexcept {
        return !(*this == rhs);
    }

  private:
    /**
     * \brief qpp::IDisplay::display() override
     *
     * Writes to the output stream a textual representation of the
     * \a qpp::internal::QCircuitRuntimeStep instance
     *
     * \param os Output stream passed by reference
     * \return Reference to the output stream
     */
    std::ostream& display(std::ostream& os) const override {
        os << runtime_type_;
        if (runtime_type_ != Type::NONE &&
            runtime_type_ != Type::SET_DITS_RUNTIME) {
            os << "-> Context [" << ctx_ << ']';
        }

        return os;
    }
};

} /* namespace internal */
} /* namespace qpp */

#endif /* QPP_INTERNAL_CLASSES_QCIRCUIT_RUNTIME_STEP_HPP_ */
