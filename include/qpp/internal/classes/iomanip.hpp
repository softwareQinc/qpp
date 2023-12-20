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
#include <complex>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

#include "qpp/options.hpp"
#include "qpp/types.hpp"

#include "qpp/classes/idisplay.hpp"
#include "qpp/internal/util.hpp"

namespace qpp::internal {

// std::ostream manipulators for formatting std::complex<> numbers
template <typename Scalar>
class IOManipComplex : public IDisplay {
    std::complex<Scalar> z_;
    IOManipComplexOpts opts_;

  public:
    IOManipComplex(std::complex<Scalar> z, IOManipComplexOpts opts)
        : z_{z}, opts_{std::move(opts)} {}

  private:
    std::ostream& display(std::ostream& os) const override {
        realT re = std::real(z_);
        realT im = std::imag(z_);
        realT abs_re = std::abs(re);
        realT abs_im = std::abs(im);

        os << opts_.left;
        // zero
        if (abs_re < opts_.chop && abs_im < opts_.chop) {
            os << '0';
        }

        // pure imaginary
        if (abs_re < opts_.chop && !(abs_im < opts_.chop)) {
            os << im << opts_.im_suffix;
        }

        // pure real
        if (!(abs_re < opts_.chop) && abs_im < opts_.chop) {
            os << re;
        }

        // complex
        if (!(abs_re < opts_.chop) && !(abs_im < opts_.chop)) {
            os << re;
            os << (im < 0 ? opts_.minus_op : opts_.plus_op);
            os << abs_im << opts_.im_suffix;
        }
        os << opts_.right;

        return os;
    }
}; /* class IOManipComplex */

// std::ostream manipulator for formatting STL/C-style containers/vectors
template <typename InputIterator>
class IOManipRange : public IDisplay {
    InputIterator first_, last_;
    IOManipRangeOpts opts_;

  public:
    explicit IOManipRange(InputIterator first, InputIterator last,
                          IOManipRangeOpts opts)
        : first_{first}, last_{last}, opts_{std::move(opts)} {}

    // to silence -Weffc++ warnings for classes that have pointer members
    // (whenever we have a pointer instantiation, i.e., iterator is a raw
    // pointer)
    IOManipRange(const IOManipRange&) = default;

    IOManipRange& operator=(const IOManipRange&) = default;

  private:
    std::ostream& display(std::ostream& os) const override {
        os << opts_.left;

        std::string sep;
        for (InputIterator it = first_; it != last_; ++it) {
            os << sep << abs_chop(*it, opts_.chop);
            sep = opts_.sep;
        }
        os << opts_.right;

        return os;
    }
}; /* class IOManipRange */

// std::ostream manipulator for formatting C-style buffers
template <typename PointerType>
class IOManipPointer : public IDisplay {
    const PointerType* p_;
    idx N_;
    IOManipPointerOpts opts_;

  public:
    explicit IOManipPointer(const PointerType* p, idx N,
                            IOManipPointerOpts opts)
        : p_{p}, N_{N}, opts_{std::move(opts)} {}

    // to silence -Weffc++ warnings for classes that have pointer members
    IOManipPointer(const IOManipPointer&) = default;

    IOManipPointer& operator=(const IOManipPointer&) = default;

  private:
    std::ostream& display(std::ostream& os) const override {
        os << opts_.left;

        for (idx i = 0; i < N_ - 1; ++i) {
            os << abs_chop(p_[i], opts_.chop) << opts_.sep;
        }
        if (N_ > 0) {
            os << abs_chop(p_[N_ - 1], opts_.chop);
        }

        os << opts_.right;

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
// std::ostream manipulator for formatting Eigen expressions
class IOManipEigen : public IDisplay, private Display_Impl_ {
#if defined(__GNUC__) && !defined(__clang__) && !defined(__INTEL_COMPILER) &&  \
    (__GNUC__ == 4) && (__GNUC_MINOR__ == 8)
#pragma GCC diagnostic pop
#endif
    cmat A_;
    IOManipEigenOpts opts_;

  public:
    // Eigen matrices
    template <typename Derived>
    explicit IOManipEigen(const Eigen::MatrixBase<Derived>& A,
                          IOManipEigenOpts opts)
        : A_{A.template cast<cplx>()}, // copy, so we can bind rvalues safely
          opts_{std::move(opts)} {}

  private:
    std::ostream& display(std::ostream& os) const override {
        return display_impl_(A_, os, opts_);
    }
}; /* class IOManipEigen */

// std::ostream manipulator for formatting expressions in Dirac notation
template <typename Scalar>
class IOManipDirac : public IDisplay {
    dirac_t<Scalar> A_;
    IOManipDiracOpts opts_{};

  public:
    explicit IOManipDirac(const dirac_t<Scalar>& A, IOManipDiracOpts opts)
        : A_{A}, opts_{std::move(opts)} {}

  private:
    void display_ket_dits_(std::ostream& os,
                           const std::vector<idx>& dits) const {
        os << "|";
        os << IOManipRange(
            dits.begin(), dits.end(),
            IOManipRangeOpts{}.set_sep("").set_left("").set_right(""));
        os << ">";
    }

    void display_bra_dits_(std::ostream& os,
                           const std::vector<idx>& dits) const {
        os << "<";
        os << IOManipRange(
            dits.begin(), dits.end(),
            IOManipRangeOpts{}.set_sep("").set_left("").set_right(""));
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
        if (A_.states.empty()) {
            return os << "0";
        }

        idx D_rows = prod(A_.dims_rows);
        idx D_cols = prod(A_.dims_cols);

        bool is_scalar = D_rows == 1 && D_cols == 1;
        bool is_col_vec = D_rows > 1 && D_cols == 1;
        bool is_row_vec = D_rows == 1 && D_cols > 1;
        bool is_matrix = D_rows > 1 && D_cols > 1;

        idx n_subsys_rows = A_.dims_rows.size();

        // display the coefficient
        auto display_coeff = [&](Scalar coeff) {
            // display parenthesis when both real and imag parts are present
            if (std::abs(std::real(coeff)) > opts_.cplx_opts.chop &&
                std::abs(std::imag(coeff)) > opts_.cplx_opts.chop) {
                os << IOManipComplex(
                    coeff,
                    IOManipComplexOpts{opts_.cplx_opts}.set_left("(").set_right(
                        ")"));
            }
            // don't display parenthesis
            else {
                os << IOManipComplex(
                    coeff,
                    IOManipComplexOpts{opts_.cplx_opts}.set_left("").set_right(
                        ""));
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
            auto abs_re = std::abs(std::real(coeff));
            auto abs_im = std::abs(std::imag(coeff));
            if (abs_re < opts_.cplx_opts.chop &&
                abs_im < opts_.cplx_opts.chop && opts_.discard_zeros) {
                continue;
            }
            auto dits = elem.second;
            if (!first) {
                os << opts_.plus_op;
            } else {
                first = false;
            }
            if (opts_.amplitudes_after) {
                display_dits(dits);
                os << opts_.mul_op;
                display_coeff(coeff);
            } else {
                display_coeff(coeff);
                os << opts_.mul_op;
                display_dits(dits);
            }
        }
        if (first) {
            os << 0;
        }

        return os;
    }
}; /* class IOManipDirac */

} /* namespace qpp::internal */

#endif /* QPP_INTERNAL_CLASSES_IOMANIP_HPP_ */
