/*
 * This file is part of qasmtools.
 *
 * Copyright (c) 2019 - 2023 softwareQ Inc. All rights reserved.
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
 *
 * Adapted from Bruno Schmitt's Tweedledum library
 */

/**
 * \file qasmtools/utils/angle.hpp
 * \brief Either symbolic or concrete representation of rotation angles
 */

#ifndef QASMTOOLS_UTILS_ANGLE_HPP_
#define QASMTOOLS_UTILS_ANGLE_HPP_

#include <cmath>
#include <iostream>
#include <numeric>
#include <optional>
#include <variant>

#include "templates.hpp"

namespace qasmtools {
namespace utils {

/**
 * \brief \f$ \pi \f$
 */
constexpr double pi = 3.141592653589793238462643383279502884;

/**
 * \brief Simple class to represent rotation angles
 *
 * A angle can be defined symbolically or numerically. When defined symbolic the
 * angle is always a multiple of pi, i.e., the symbolic value will always
 * multiplied by pi.
 *
 * The numeric value of a rotation angle is given in radians (rad).
 */
class Angle {
    using fraction = std::pair<int, int>;

    std::variant<fraction, double> value_;

  public:
    constexpr Angle(int n, int d) : value_(std::make_pair(n, d)) {
        if (d == 0) {
            throw std::invalid_argument(
                "Trying to construct angle with denominator 0");
        }

        normalize();
    }
    constexpr Angle(fraction angle) : value_(angle) {
        if (angle.second == 0) {
            throw std::invalid_argument(
                "Trying to construct angle with denominator 0");
        }

        normalize();
    }
    constexpr Angle(double angle) : value_(angle) {}

    /** \brief Returns true if this angle is symbolically defined. */
    constexpr bool is_symbolic() const {
        return std::holds_alternative<fraction>(value_);
    }

    /*! \brief Returns true if this angle is symbolically defined. */
    constexpr bool is_numeric() const {
        return std::holds_alternative<double>(value_);
    }

    /*! \brief Returns the symbolic value of this angle. */
    constexpr std::optional<fraction> symbolic_value() const {
        if (std::holds_alternative<fraction>(value_)) {
            return std::get<fraction>(value_);
        } else {
            return std::nullopt;
        }
    }

    /*! \brief Returns the numeric value of this angle. */
    constexpr double numeric_value() const {
        if (std::holds_alternative<fraction>(value_)) {
            auto frac = std::get<fraction>(value_);
            return ((double)frac.first * pi) / (double)frac.second;
        } else {
            return std::get<double>(value_);
        }
    }

    constexpr Angle operator-() const {
        if (std::holds_alternative<fraction>(value_)) {
            auto frac = std::get<fraction>(value_);
            return Angle(2 * frac.second - frac.first, frac.second);
        } else {
            return Angle(-std::get<double>(value_));
        }
    }

    bool operator==(const Angle& other) const {
        return numeric_value() == other.numeric_value();
    }

    bool operator!=(const Angle& other) const {
        return numeric_value() != other.numeric_value();
    }

    Angle& operator+=(const Angle& rhs) {
        if (is_symbolic() && rhs.is_symbolic()) {
            auto [a, b] = std::get<fraction>(value_);
            auto [c, d] = std::get<fraction>(rhs.value_);
            value_ = fraction(a * d + c * b, b * d);
            normalize();
        } else {
            auto a = numeric_value();
            auto b = rhs.numeric_value();
            value_ = a + b;
        }

        return *this;
    }

    Angle& operator-=(const Angle& rhs) {
        *this += -rhs;
        return *this;
    }

    Angle& operator*=(int fac) {
        if (is_symbolic()) {
            auto [a, b] = std::get<fraction>(value_);
            value_ = fraction(a * fac, b);
            normalize();
        } else {
            auto a = numeric_value();
            value_ = a * (double)fac;
        }

        return *this;
    }

    Angle& operator/=(int div) {
        if (is_symbolic()) {
            auto [a, b] = std::get<fraction>(value_);
            value_ = fraction(a, b * div);
            normalize();
        } else {
            auto a = numeric_value();
            value_ = a / (double)div;
        }

        return *this;
    }

    friend Angle operator+(const Angle& lhs, const Angle& rhs) {
        auto tmp = lhs;
        tmp += rhs;
        return tmp;
    }

    friend Angle operator-(const Angle& lhs, const Angle& rhs) {
        return lhs + (-rhs);
    }

    friend Angle operator*(const Angle& lhs, int rhs) {
        auto tmp = lhs;
        tmp *= rhs;
        return tmp;
    }

    friend Angle operator/(const Angle& lhs, int rhs) {
        auto tmp = lhs;
        tmp /= rhs;
        return tmp;
    }

    friend std::ostream& operator<<(std::ostream& os, const Angle& angle) {
        if (angle.is_symbolic()) {
            auto [a, b] = std::get<fraction>(angle.value_);
            switch (a) {
                case 1:
                    break;
                case -1:
                    os << '-';
                    break;
                default:
                    os << a << '*';
                    break;
            }
            os << "pi";
            if (b != 1) {
                os << "/" << b;
            }
        } else {
            os << std::get<double>(angle.value_);
        }

        return os;
    }

  private:
    constexpr void normalize() {
        if (is_numeric()) {
            return;
        }

        auto& [n, d] = std::get<fraction>(value_);
        if (n == 0) {
            d = 1;
        } else {
            // calculate the sign & get absolute values
            int sgn = 1;
            if (n < 0) {
                sgn = -1;
                n = std::abs(n);
            }
            if (d < 0) {
                sgn *= -1;
                d = std::abs(d);
            }

            // bring into lowest form
            auto tmp = std::gcd(n, d);
            n = n / tmp;
            d = d / tmp;

            // bring into [0, 2pi) range
            n = n % (2 * d);
            if (sgn == -1) {
                n = (2 * d) - n;
            }
        }
    }
};

namespace angles {
/*! \brief identity */
constexpr Angle zero(0, 1);
/*! \brief rotation angle of a T gate */
constexpr Angle pi_quarter(1, 4);
/*! \brief rotation angle of a S gate (phase gate) */
constexpr Angle pi_half(1, 2);
/*! \brief rotation angle of a Pauli-Z gate, Pauli-X (NOT) */
constexpr Angle pi(1, 1);
} /* namespace angles */

} /* namespace utils */
} /* namespace qasmtools */

#endif /* QASMTOOLS_UTILS_ANGLE_HPP_ */
