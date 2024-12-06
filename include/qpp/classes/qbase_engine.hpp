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
 * \file qpp/classes/qbase_engine.hpp
 * \brief Base class for quantum engines
 */

#ifndef QPP_CLASSES_QBASE_ENGINE_HPP_
#define QPP_CLASSES_QBASE_ENGINE_HPP_

#include <memory>
#include <ostream>
#include <sstream>
#include <string>
#include <type_traits>

#include "qpp/types.hpp"

#include "qpp/classes/idisplay.hpp"
#include "qpp/classes/ijson.hpp"
#include "qpp/classes/qcircuit_traits.hpp"
#include "qpp/classes/qengine_traits.hpp"

namespace qpp {
/**
 * \class QBaseEngine
 * \brief Base class for all quantum engines
 *
 * \tparam T Engine's state underlying type
 * \tparam QCT Circuit underlying type
 *
 */
template <typename T, typename QCT>
class QBaseEngine : public IQEngineTraits, public IDisplay, public IJSON {
  protected:
    const QCT* qc_ptr_; ///< pointer to constant quantum circuit description
  public:
    /**
     * \brief Constructs a quantum engine out of a quantum circuit description
     *
     * \note The quantum circuit description must be an lvalue
     * \see qpp::QBaseEngineT(QCT&&)
     *
     * \param qc Quantum circuit description
     */
    explicit QBaseEngine(const QCT& qc) : qc_ptr_{std::addressof(qc)} {}

    // silence -Weffc++ class has pointer data members
    /**
     * \brief Default copy constructor
     */
    QBaseEngine(const QBaseEngine&) = default;

    // silence -Weffc++ class has pointer data members
    /**
     * \brief Default copy assignment operator
     *
     * \return Reference to the current instance
     */
    QBaseEngine& operator=(const QBaseEngine&) = default;

    /**
     * \brief Disables rvalue QCT
     */
    QBaseEngine(QCT&&) = delete;

    /**
     * \brief Default virtual destructor
     */
    virtual ~QBaseEngine() override = default;

    // traits
    /**
     * \brief qpp::IQEngineTraits::traits_get_name() override
     */
    std::string traits_get_name() const override { return "QBaseEngine"; }

    /**
     * \brief qpp::IQEngineTraits::traits_is_noisy() override
     */
    bool traits_is_noisy() const override { return false; }

    /**
     * \brief qpp::IQEngineTraits::traits_is_pure() override
     */
    bool traits_is_pure() const override {
        if (std::is_same_v<T, cmat>) {
            return false;
        } else if (std::is_same_v<T, ket>) {
            return true;
        }
        // default, we assume pure states
        return true;
    }
    // end traits

    /**
     * \brief Executes one step in the quantum circuit description
     *
     * \note This function is a no-op in qpp::QBaseEngine; override it in all
     * derived classes to achieve the desired behaviour
     *
     * \param it Iterator to the step to be executed
     * \return Reference to the current instance
     */
    virtual QBaseEngine&
    execute([[maybe_unused]] typename QCircuitTraits<QCT>::iterator_type it) {
        return *this;
    };

    QBaseEngine& execute([[maybe_unused]]
                         typename QCircuitTraits<QCT>::value_type& step) {}

    /**
     * \brief Executes the entire quantum circuit description
     * \note This function is a no-op in qpp::QBaseEngine; override it in all
     * derived classes to achieve the desired behaviour
     *
     * \param reps Number of repetitions
     * \return Reference to the current instance
     */
    virtual QBaseEngine& execute([[maybe_unused]] idx reps = 1) {
        return *this;
    };

    /**
     * \brief Underlying quantum state
     * \note This function is a no-op in qpp::QBaseEngine; override it in all
     * derived classes to achieve the desired behaviour
     *
     * \return Underlying quantum state
     */
    virtual T get_state() const { return T{}; }

    /**
     * \brief Sets the underlying quantum state to \a state
     * \note This function is a no-op in qpp::QBaseEngine; override it in all
     * derived classes to achieve the desired behaviour
     *
     * \param state Quantum state
     * \return Reference to the current instance
     */
    virtual QBaseEngine& set_state([[maybe_unused]] const T& state) {
        return *this;
    }

    /**
     * \brief Quantum circuit description, lvalue ref qualifier
     *
     * \return Const reference to the underlying quantum circuit description
     */
    const QCT& get_circuit() const& noexcept { return *qc_ptr_; }

    /**
     * \brief Quantum circuit description, rvalue ref qualifier
     *
     * \return Copy of the underlying quantum circuit description
     */
    QCT get_circuit() const&& noexcept { return *qc_ptr_; }

    std::string to_JSON(bool enclosed_in_curly_brackets = true) const override {
        std::string result;

        if (enclosed_in_curly_brackets) {
            result += "{";
        }

        std::ostringstream ss;
        ss << "\"Engine\": ";
        ss << "{\"name\": " << std::quoted(traits_get_name()) << ", ";
        ss << "\"is_noisy\": " << (traits_is_noisy() ? "true" : "false")
           << ", ";
        ss << "\"is_pure\": " << (traits_is_pure() ? "true" : "false");
        ss << "}, ";
        ss << "\"Circuit\": ";
        ss << "{\"nq\": " << get_circuit().get_nq() << ", ";
        ss << "\"nc\": " << get_circuit().get_nc() << ", ";
        ss << "\"d\": " << get_circuit().get_d() << ", ";
        ss << "\"name\": ";
        if (get_circuit().get_name().has_value()) {
            ss << "\"" << get_circuit().get_name().value() << "\"";
        } else {
            ss << "null";
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
     * Writes to the output stream a textual representation of the state of
     * the engine
     *
     * \param os Output stream passed by reference
     * \return Reference to the output stream
     */
    std::ostream& display(std::ostream& os) const override {
        os << '[' << traits_get_name() << ']';
        os << '\n';

        os << "<Circuit nq: " << get_circuit().get_nq()
           << ", nc: " << get_circuit().get_nc()
           << ", d: " << get_circuit().get_d();

        if (get_circuit().get_name().has_value()) {
            os << ", name: ";
            os << "\"" << get_circuit().get_name().value() << "\"";
        }
        os << ">\n";

        return os;
    }
};
} /* namespace qpp */

#endif /* QPP_CLASSES_QBASE_ENGINE_HPP_ */
