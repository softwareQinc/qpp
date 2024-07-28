
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
 * \file qpp/internal/classes/qcircuit_measurement_step.hpp
 * \brief qpp::internal::QCircuitMeasurementStep
 */

#ifndef QPP_INTERNAL_CLASSES_QCIRCUIT_MEASUREMENT_STEP_HPP_
#define QPP_INTERNAL_CLASSES_QCIRCUIT_MEASUREMENT_STEP_HPP_

#include "qpp/input_output.hpp"

#include "qpp/classes/idisplay.hpp"

namespace qpp {
namespace internal {
/**
 * \brief One step consisting only of measurements in the circuit
 */
struct QCircuitMeasurementStep : IDisplay {
    /**
     * \brief Type of measurement being executed in a measurement step
     */
    enum class Type {
        NONE, ///< represents no measurement

        MEASURE, ///< Z measurement of single qudit

        MEASURE_MANY, ///< Z measurement of multiple qudits

        MEASURE_V, ///< measurement of single qudit in the orthonormal
                   ///< basis or rank-1 projectors specified by the
                   ///< columns of matrix \a V

        MEASURE_V_JOINT, ///< joint measurement of multiple qudits in
                         ///< the orthonormal basis or rank-1 projectors
                         ///< specified by the columns of the matrix \a
                         ///< V

        MEASURE_ND, ///< Z measurement of single qudit, non-destructive

        MEASURE_MANY_ND, ///< Z measurement of multiple qudits,
                         ///< non-destructive

        MEASURE_V_ND, ///< measurement of single qudit in the
                      ///< orthonormal basis or rank-1 projectors
                      ///< specified by the columns of matrix \a V,
                      ///< non-destructive

        MEASURE_V_JOINT_ND, ///< joint measurement of multiple qudits in
                            ///< the orthonormal basis or rank-1
                            ///< projectors specified by the columns of
                            ///< the matrix \a V, non-destructive

        RESET, ///< resets single qudit

        RESET_MANY, ///< resets multiple qudits

        DISCARD, ///< discards single qudit

        DISCARD_MANY, ///< discards multiple qudits

        POST_SELECT, ///< post-selects single qudit

        POST_SELECT_MANY, ///< post-selects multiple qudits

        POST_SELECT_V, ///< post-selection of single qudit in the
                       ///< orthonormal basis or rank-1 projectors specified
                       ///< by the columns of matrix \a V

        POST_SELECT_V_JOINT, ///< joint post-selection of multiple qudits in
                             ///< the orthonormal basis or rank-1 projectors
                             ///< specified by the columns of the matrix \a
                             ///< V

        POST_SELECT_ND, ///< post-selects single qudit, non-destructive

        POST_SELECT_MANY_ND, ///< post-selects multiple qudits,
                             ///< non-destructive

        POST_SELECT_V_ND, ///< post-selection of single qudit in the
                          ///< orthonormal basis or rank-1 projectors
                          ///< specified by the columns of matrix \a V,
                          ///< non-destructive

        POST_SELECT_V_JOINT_ND, ///< joint post-selection of multiple qudits
                                ///< in the orthonormal basis or rank-1
                                ///< projectors specified by the columns of
                                ///< the matrix \a V, non-destructive
    };

    /**
     * \brief Extraction operator overload for
     * qpp::internal::QCircuitMeasurementStep::Type enum class
     *
     * \param os Output stream passed by reference
     * \param measure_type qpp::internal::QCircuitMeasurementStep::Type enum
     * class
     * \return Reference to the output stream
     */
    friend std::ostream& operator<<(std::ostream& os,
                                    const Type& measure_type) {
        switch (measure_type) {
            case Type::NONE:
                os << "MEASURE NONE";
                break;
            case Type::MEASURE:
                os << "MEASURE";
                break;
            case Type::MEASURE_MANY:
                os << "MEASURE_MANY";
                break;
            case Type::MEASURE_V:
                os << "MEASURE_V";
                break;
            case Type::MEASURE_V_JOINT:
                os << "MEASURE_V_JOINT";
                break;
            case Type::MEASURE_ND:
                os << "MEASURE_ND";
                break;
            case Type::MEASURE_MANY_ND:
                os << "MEASURE_MANY_ND";
                break;
            case Type::MEASURE_V_ND:
                os << "MEASURE_V_ND";
                break;
            case Type::MEASURE_V_JOINT_ND:
                os << "MEASURE_V_JOINT_ND";
                break;
            case Type::RESET:
                os << "RESET";
                break;
            case Type::RESET_MANY:
                os << "RESET_MANY";
                break;
            case Type::DISCARD:
                os << "DISCARD";
                break;
            case Type::DISCARD_MANY:
                os << "DISCARD_MANY";
                break;
            case Type::POST_SELECT:
                os << "POST_SELECT";
                break;
            case Type::POST_SELECT_MANY:
                os << "POST_SELECT_MANY";
                break;
            case Type::POST_SELECT_V:
                os << "POST_SELECT_V";
                break;
            case Type::POST_SELECT_V_JOINT:
                os << "POST_SELECT_V_JOINT";
                break;
            case Type::POST_SELECT_ND:
                os << "POST_SELECT_ND";
                break;
            case Type::POST_SELECT_MANY_ND:
                os << "POST_SELECT_MANY_ND";
                break;
            case Type::POST_SELECT_V_ND:
                os << "POST_SELECT_V_ND";
                break;
            case Type::POST_SELECT_V_JOINT_ND:
                os << "POST_SELECT_V_JOINT_ND";
                break;
        }

        return os;
    }

    Type measurement_type_ = Type::NONE;   ///< measurement type
    std::vector<std::size_t> mats_hash_{}; ///< hashes of measurement
                                           ///< matrix/matrices
    std::vector<idx> target_{}; ///< target where the measurement is applied
    idx c_reg_{};               ///< index of the classical register where the
                                ///< measurement result is being stored
    std::optional<std::vector<idx>> ps_vals_{}; ///< post-selection values
    std::optional<std::string> name_{}; ///< custom name of the measurement(s)

    /**
     * \brief Default constructor
     */
    QCircuitMeasurementStep() = default;

    /**
     * \brief Constructs a measurement step instance
     *
     * \param measurement_type Measurement type
     * \param mats_hash Vector of hashes of the measurement matrix/matrices
     * \param target Target qudit indexes
     * \param c_reg Classical register where the value of the measurement is
     * stored
     * \param name Optional gate name
     * \param ps_vals Optional post-selection values
     */
    explicit QCircuitMeasurementStep(
        Type measurement_type, std::vector<std::size_t> mats_hash,
        std::vector<idx> target, idx c_reg,
        std::optional<std::string> name = std::nullopt,
        std::optional<std::vector<idx>> ps_vals = std::nullopt)
        : measurement_type_{measurement_type}, mats_hash_{std::move(mats_hash)},
          target_{std::move(target)}, c_reg_{c_reg},
          ps_vals_{std::move(ps_vals)}, name_{std::move(name)} {}

    /**
     * \brief Equality operator
     * \note Ignores measurement names
     *
     * \param rhs qpp::internal::QCircuitMeasurementStep against which the
     * equality is being tested
     * \return True if the qpp::internal::QCircuitMeasurementStep(s) are equal,
     * false otherwise
     */
    bool operator==(const QCircuitMeasurementStep& rhs) const noexcept {
        return std::tie(rhs.target_, rhs.measurement_type_, rhs.mats_hash_,
                        rhs.c_reg_,
                        rhs.ps_vals_) == std::tie(target_, measurement_type_,
                                                  mats_hash_, c_reg_, ps_vals_);
    }

    /**
     * \brief Inequality operator
     * \note Ignores measurement names
     *
     * \param rhs qpp::internal::QCircuitMeasurementStep against which the
     * inequality is being tested
     * \return True if the qpp::internal::QCircuitMeasurementStep(s) are not
     * equal, false otherwise
     */
    bool operator!=(const QCircuitMeasurementStep& rhs) const noexcept {
        return !(*this == rhs);
    }

  private:
    /**
     * \brief qpp::IDisplay::display() override
     *
     * Writes to the output stream a textual representation of the
     * \a qpp::internal::QCircuitMeasurementStep instance
     *
     * \param os Output stream passed by reference
     * \return Reference to the output stream
     */
    std::ostream& display(std::ostream& os) const override {
        os << measurement_type_ << ", ";
        os << "target = "
           << disp(target_, IOManipContainerOpts{}.set_sep(", "));
        if (ps_vals_.has_value()) {
            os << ", ps_vals = "
               << disp(ps_vals_.value(), IOManipContainerOpts{}.set_sep(", "));
        }
        if (measurement_type_ != Type::RESET &&
            measurement_type_ != Type::RESET_MANY &&
            measurement_type_ != Type::DISCARD &&
            measurement_type_ != Type::DISCARD_MANY) {
            os << ", c_reg = " << c_reg_;
        }
        if (name_.has_value()) {
            os << ", name = " << std::quoted(name_.value());
        }

        return os;
    }
};

} /* namespace internal */
} /* namespace qpp */

#endif /* QPP_INTERNAL_CLASSES_QCIRCUIT_MEASUREMENT_STEP_HPP_ */
