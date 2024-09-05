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
 * \file qpp/internal/classes/qengine_statistics.hpp
 * \brief Quantum engine statistics
 */

#ifndef QPP_INTERNAL_CLASSES_QENGINE_STATISTICS_HPP_
#define QPP_INTERNAL_CLASSES_QENGINE_STATISTICS_HPP_

#include <map>
#include <ostream>
#include <string>
#include <vector>

#include "qpp/input_output.hpp"

#include "qpp/classes/idisplay.hpp"
#include "qpp/classes/ijson.hpp"

namespace qpp {
namespace internal {
/**
 * \class qpp::internal::QEngineStatistics
 * \brief Sampling/measurement statistics
 */
class QEngineStatistics : public IDisplay, public IJSON {
    /**
     * \brief Measurement/sampling statistics
     */
    using stats_t_ = std::map<std::vector<idx>, idx>;
    mutable stats_t_ stats_data_{}; ///< statistics data

  public:
    /**
     * \brief Default constructor
     */
    QEngineStatistics() = default;

    /**
     * \brief Constructor
     *
     * \param stats Instance of qpp::QEngine::Statistics
     */
    explicit QEngineStatistics(stats_t_ stats)
        : stats_data_{std::move(stats)} {}

    /**
     * \brief Number of samples
     *
     * \return Number of samples
     */
    idx get_num_reps() const {
        idx result = 0;
        for (auto&& [_, val] : stats_data_) {
            result += val;
        }
        return result;
    }

    /**
     * \brief Number of distinct outcomes
     *
     * \return Number of distinct outcomes
     */
    idx get_num_outcomes() const { return stats_data_.size(); }

    /**
     * \brief Raw data structure representing the statistics
     *
     * \return Raw data structure representing the statistics
     */
    stats_t_& data() const& { return stats_data_; }

    /**
     * \brief qpp::IJSON::to_JSON() override
     *
     * Displays the statistics in JSON format
     *
     * \param enclosed_in_curly_brackets If true, encloses the result in
     * curly brackets
     * \return String containing the JSON representation of the statistics
     */
    std::string to_JSON(bool enclosed_in_curly_brackets = true) const override {
        std::string result;

        if (enclosed_in_curly_brackets) {
            result += "{";
        }

        std::ostringstream ss;
        ss << "\"num_reps\": " << get_num_reps() << ", ";
        ss << "\"num_outcomes\": " << get_num_outcomes() << ", ";

        bool is_first = true;
        std::string sep{};
        ss << "\"outcomes\": ";
        ss << "{";
        for (auto&& [key, val] : stats_data_) {
            ss << sep << "\"" << disp(key, IOManipContainerOpts{}.set_sep(", "))
               << "\": " << val;
            if (is_first) {
                is_first = false;
                sep = ", ";
            }
        }
        ss << "}";
        result += ss.str();

        if (enclosed_in_curly_brackets) {
            result += "}";
        }

        return result;
    }

  private:
    /**
     * \brief qpp::IDisplay::display() override
     *
     * Writes to the output stream a textual representation of the
     * statistics
     *
     * \param os Output stream passed by reference
     * \return Reference to the output stream
     */
    std::ostream& display(std::ostream& os) const override {
        os << "[Statistics]\n";
        os << "\tnum_reps: " << get_num_reps() << '\n';
        os << "\tnum_outcomes: " << get_num_outcomes() << '\n';
        bool is_first = true;
        std::string sep{};
        for (auto&& [key, val] : stats_data_) {
            os << sep << '\t' << disp(key, IOManipContainerOpts{}.set_sep(" "))
               << ": " << val;
            if (is_first) {
                is_first = false;
                sep = '\n';
            }
        }

        return os;
    };
}; /* class QEngineStatistics */

} /* namespace internal */
} /* namespace qpp */

#endif /* QPP_INTERNAL_CLASSES_QENGINE_STATISTICS_HPP_ */
