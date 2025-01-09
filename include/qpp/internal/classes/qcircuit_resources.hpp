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
 * \file qpp/internal/classes/qcircuit_resources.hpp
 * \brief qpp::internal::QCircuitResources
 */

#ifndef QPP_INTERNAL_CLASSES_QCIRCUIT_RESOURCES_HPP_
#define QPP_INTERNAL_CLASSES_QCIRCUIT_RESOURCES_HPP_

#include "qpp/classes/idisplay.hpp"
#include "qpp/classes/ijson.hpp"

namespace qpp {
namespace internal {
/**
 * \brief Quantum circuit resources
 */
struct QCircuitResources : IDisplay, IJSON {
    idx nq{}, nc{}, d{};
    std::optional<std::string> name{};
    idx step_count{};
    idx gate_count{};
    idx gate_depth{};
    idx measurement_count{};
    idx measurement_depth{};
    idx total_depth{};

    /**
     * \brief qpp::IJSON::to_JSON() override
     *
     * Displays the quantum circuit resources in JSON format
     *
     * \param enclosed_in_curly_brackets If true, encloses the result in
     * curly brackets
     * \return String containing the JSON representation of the quantum
     * circuit resources
     */
    std::string to_JSON(bool enclosed_in_curly_brackets = true) const override {
        std::string result;

        if (enclosed_in_curly_brackets) {
            result += "{";
        }

        result += "\"nq\": " + std::to_string(nq) + ", ";
        result += "\"nc\": " + std::to_string(nc) + ", ";
        result += "\"d\": " + std::to_string(d) + ", ";

        result += "\"step count\": " + std::to_string(step_count) + ", ";
        result += "\"total gate count\": " + std::to_string(gate_count) + ", ";
        result += "\"total gate depth\": " + std::to_string(gate_depth) + ", ";
        result += "\"total measurement count\": " +
                  std::to_string(measurement_count) + ", ";
        result += "\"total measurement depth\": " +
                  std::to_string(measurement_depth) + ", ";
        result += "\"total depth\": " + std::to_string(total_depth);

        if (enclosed_in_curly_brackets) {
            result += "}";
        }

        return result;
    }

  private:
    /**
     * \brief qpp::IDisplay::display() override
     *
     * Writes to the output stream a textual representation of the quantum
     * circuit resources
     *
     * \param os Output stream passed by reference
     * \return Reference to the output stream
     */
    std::ostream& display(std::ostream& os) const override {
        os << "[QCircuitResources]\n";

        os << "<QCircuit " << "nq: " << nq << ", nc: " << nc << ", d: " << d;

        if (name.has_value()) {
            os << ", name: " << std::quoted(name.value());
        }
        os << ">\n";

        os << "step count: " << step_count << '\n';
        os << "total gate count: " << gate_count << '\n';
        os << "total gate depth: " << gate_depth << '\n';
        os << "total measurement count: " << measurement_count << '\n';
        os << "total measurement depth: " << measurement_depth << '\n';
        os << "total depth: " << total_depth;

        return os;
    }
};

} /* namespace internal */
} /* namespace qpp */

#endif /* QPP_INTERNAL_CLASSES_QCIRCUIT_RESOURCES_HPP_ */
