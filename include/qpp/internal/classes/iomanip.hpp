/*
 * This file is part of Quantum++.
 *
 * Copyright (c) 2013 - 2023 softwareQ Inc. All rights reserved.
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
 * \file internal/classes/iomanip.hpp
 * \brief Input/output manipulators
 */

#ifndef QPP_INTERNAL_CLASSES_IOMANIP_HPP_
#define QPP_INTERNAL_CLASSES_IOMANIP_HPP_

#include <cmath>
#include <string>
#include <utility>

#include "qpp/constants.hpp"
#include "qpp/types.hpp"

#include "qpp/classes/idisplay.hpp"
#include "qpp/internal/util.hpp"

namespace qpp::internal {
// ostream manipulators for nice formatting of
// Eigen matrices and STL/C-style containers/vectors

template <typename InputIterator>
class IOManipRange : public IDisplay {
    InputIterator first_, last_;
    std::string separator_, start_, end_;
    realT chop_;

  public:
    explicit IOManipRange(InputIterator first, InputIterator last,
                          std::string separator, std::string start = "[",
                          std::string end = "]", realT chop = qpp::chop)
        : first_{first}, last_{last}, separator_{std::move(separator)},
          start_{std::move(start)}, end_{std::move(end)}, chop_{chop} {}

    // to silence -Weffc++ warnings for classes that have pointer members
    // (whenever we have a pointer instantiation, i.e., iterator is a raw
    // pointer)
    IOManipRange(const IOManipRange&) = default;

    IOManipRange& operator=(const IOManipRange&) = default;

  private:
    std::ostream& display(std::ostream& os) const override {
        os << start_;

        std::string sep;
        for (InputIterator it = first_; it != last_; ++it) {
            os << sep << abs_chop(*it, chop_);
            sep = separator_;
        }
        os << end_;

        return os;
    }
}; /* class IOManipRange */

template <typename PointerType>
class IOManipPointer : public IDisplay {
    const PointerType* p_;
    idx N_;
    std::string separator_, start_, end_;
    realT chop_;

  public:
    explicit IOManipPointer(const PointerType* p, idx N, std::string separator,
                            std::string start = "[", std::string end = "]",
                            realT chop = qpp::chop)
        : p_{p}, N_{N}, separator_{std::move(separator)},
          start_{std::move(start)}, end_{std::move(end)}, chop_{chop} {}

    // to silence -Weffc++ warnings for classes that have pointer members
    IOManipPointer(const IOManipPointer&) = default;

    IOManipPointer& operator=(const IOManipPointer&) = default;

  private:
    std::ostream& display(std::ostream& os) const override {
        os << start_;

        for (idx i = 0; i < N_ - 1; ++i)
            os << abs_chop(p_[i], chop_) << separator_;
        if (N_ > 0)
            os << abs_chop(p_[N_ - 1], chop_);

        os << end_;

        return os;
    }
}; /* class IOManipPointer */

// silence g++4.8.x bogus warning -Wnon-virtual-dtor for
// qpp::internal::Display_impl_ when class qpp::internal::IOManipEigen
// privately inherits from it
#if defined(__GNUC__) && !defined(__clang__) && !defined(__INTEL_COMPILER) &&  \
    (__GNUC__ == 4) && (__GNUC_MINOR__ == 8)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Weffc++"
#endif
class IOManipEigen : public IDisplay, private Display_Impl_ {
#if defined(__GNUC__) && !defined(__clang__) && !defined(__INTEL_COMPILER) &&  \
    (__GNUC__ == 4) && (__GNUC_MINOR__ == 8)
#pragma GCC diagnostic pop
#endif
    cmat A_;
    realT chop_;

  public:
    // Eigen matrices
    template <typename Derived>
    explicit IOManipEigen(const Eigen::MatrixBase<Derived>& A,
                          realT chop = qpp::chop)
        : A_{A.template cast<cplx>()}, // copy, so we can bind rvalues safely
          chop_{chop} {}

    // Complex numbers
    explicit IOManipEigen(const cplx z, realT chop = qpp::chop)
        : A_{cmat::Zero(1, 1)}, chop_{chop} {
        // put the complex number inside an Eigen matrix
        A_(0, 0) = z;
    }

  private:
    std::ostream& display(std::ostream& os) const override {
        return display_impl_(A_, os, chop_);
    }
}; /* class IOManipEigen */

template <typename Scalar>
class IOManipDirac : public IDisplay {
    io_braket<Scalar> A_;
    bool normal_form_;
    std::string add_op_;
    std::string mult_op_;
    realT chop_;

  public:
    explicit IOManipDirac(const io_braket<Scalar>& A, bool normal_form,
                          const std::string& add_op, const std::string& mult_op,
                          realT chop = qpp::chop)
        : A_{A}, normal_form_{normal_form}, add_op_{add_op}, mult_op_{mult_op},
          chop_{chop} {}

  private:
    void display_ket_dits_(std::ostream& os,
                           const std::vector<idx>& dits) const {
        os << "|";
        os << IOManipRange(dits.begin(), dits.end(), "", "", "");
        os << ">";
    }

    void display_bra_dits_(std::ostream& os,
                           const std::vector<idx>& dits) const {
        os << "<";
        os << IOManipRange(dits.begin(), dits.end(), "", "", "");
        os << "|";
    }

    void display_mat_dits_(std::ostream& os, const std::vector<idx>& dits,
                           idx n_subsys_rows) const {
        auto split_it = std::next(dits.begin(), n_subsys_rows);
        std::vector<idx> row_dits(dits.begin(), split_it);
        std::vector<idx> col_dits(split_it, dits.end());
        display_ket_dits_(os, row_dits);
        display_bra_dits_(os, col_dits);
    }

    std::ostream& display(std::ostream& os) const override {
        if (A_.states.empty())
            return os << "0";

        idx D_rows = prod(A_.dims_rows);
        idx D_cols = prod(A_.dims_cols);

        bool is_scalar = D_rows == 1 && D_cols == 1;
        bool is_col_vec = D_rows > 1 && D_cols == 1;
        bool is_row_vec = D_rows == 1 && D_cols > 1;
        bool is_matrix = D_rows > 1 && D_cols > 1;

        idx n_subsys_rows = A_.dims_rows.size();

        // display the coefficient
        auto display_coeff = [&](Scalar coeff) {
            if (std::abs(std::real(coeff)) > chop_ &&
                std::abs(std::imag(coeff)) > chop_) {
                os << '(' << IOManipEigen(coeff, chop_) << ')';
            } else {
                os << IOManipEigen(coeff, chop_);
            }
        };

        // display the dits
        auto display_dits = [&](const std::vector<idx>& dits) {
            // ket
            if (is_col_vec) {
                display_ket_dits_(os, dits);
            }
            // bra
            else if (is_row_vec) {
                display_bra_dits_(os, dits);
            }
            // density matrix
            else if (is_matrix) {
                display_mat_dits_(os, dits, n_subsys_rows);
            }
        };

        if (is_scalar) {
            display_coeff(A_.states[0].first);
            return os;
        }

        bool first = true;
        for (auto&& elem : A_.states) {
            auto coeff = elem.first;
            auto dits = elem.second;
            if (!first) {
                os << add_op_;
            } else {
                first = false;
            }
            if (normal_form_) {
                display_coeff(coeff);
                os << mult_op_;
                display_dits(dits);
            } else {
                display_dits(dits);
                os << mult_op_;
                display_coeff(coeff);
            }
        }

        return os;
    }
}; /* class IOManipDirac */

} /* namespace qpp::internal */

#endif /* QPP_INTERNAL_CLASSES_IOMANIP_HPP_ */
