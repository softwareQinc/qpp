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
 * \file options.hpp
 * \brief Options
 */

#ifndef QPP_OPTIONS_HPP_
#define QPP_OPTIONS_HPP_

#include <string>

#include "qpp/types.hpp"

namespace qpp {
/**
 * \brief  Used in qpp::disp() for setting to zero  numbers that have their
 * absolute value smaller than qpp::chop
 */
constexpr realT chop = 1e-14;

/*
 * \brief Formatting options for scalars
 * \see qpp::disp()
 */
struct IOManipScalarOpts {
    std::string left = "";  ///< left delimiter
    std::string right = ""; ///< right delimiter
    realT chop = qpp::chop; ///<  sets to zero coefficients smaller than this
                            ///<  w.r.t. some norm
    // setters

    /* \brief Sets left delimiter
     *
     * \param left_delimiter Left delimiter
     * \return Reference to the current instance
     */
    IOManipScalarOpts& set_left(std::string left_delimiter) {
        left = std::move(left_delimiter);
        return *this;
    }

    /* \brief Sets right delimiter
     *
     * \param right_delimiter Right delimiter
     * \return Reference to the current instance
     */
    IOManipScalarOpts& set_right(std::string right_delimiter) {
        right = std::move(right_delimiter);
        return *this;
    }

    /* \brief Sets chopping threshold
     *
     * \param chop_at Chopping threshold, sets to zero coefficients (real or
     * imaginary) smaller than this
     * \return Reference to the current instance
     */
    IOManipScalarOpts& set_chop(realT chop_at) {
        chop = chop_at;
        return *this;
    }
};

/*
 * \brief Formatting options for std::complex<T>
 * \see qpp::disp()
 */
struct IOManipComplexOpts {
    std::string im_suffix =
        "i"; ///< imaginary symbol for the complex number "i"
    std::string plus_op = " + ";  ///< binary addition operator
    std::string minus_op = " - "; ///< binary subtraction operator
    std::string left = "";        ///< left delimiter
    std::string right = "";       ///< right delimiter
    realT chop = qpp::chop;       ///<  sets to zero coefficients (real or
                                  ///<  imaginary) smaller than this
    // setters

    /* \brief Sets imaginary symbol
     *
     * \param imaginary_suffix Imaginary symbol for the complex number "i"
     * \return Reference to the current instance
     */
    IOManipComplexOpts& set_im_suffix(std::string imaginary_suffix) {
        im_suffix = std::move(imaginary_suffix);
        return *this;
    }

    /* \brief Sets plus operator
     *
     * \param plus_operator Plus operator
     * \return Reference to the current instance
     */
    IOManipComplexOpts& set_plus_op(std::string plus_operator) {
        plus_op = std::move(plus_operator);
        return *this;
    }

    /* \brief Sets minus operator
     *
     * \param minus_operator Minus operator
     * \return Reference to the current instance
     */
    IOManipComplexOpts& set_minus_op(std::string minus_operator) {
        minus_op = std::move(minus_operator);
        return *this;
    }

    /* \brief Sets left delimiter
     *
     * \param left_delimiter Left delimiter
     * \return Reference to the current instance
     */
    IOManipComplexOpts& set_left(std::string left_delimiter) {
        left = std::move(left_delimiter);
        return *this;
    }

    /* \brief Sets right delimiter
     *
     * \param right_delimiter Right delimiter
     * \return Reference to the current instance
     */
    IOManipComplexOpts& set_right(std::string right_delimiter) {
        right = std::move(right_delimiter);
        return *this;
    }

    /* \brief Sets chopping threshold
     *
     * \param chop_at Chopping threshold, sets to zero coefficients (real or
     * imaginary) smaller than this
     * \return Reference to the current instance
     */
    IOManipComplexOpts& set_chop(realT chop_at) {
        chop = chop_at;
        return *this;
    }
};

/*
 * \brief Formatting options for Eigen::MatrixBase<Derived>
 * \see qpp::disp()
 */
struct IOManipEigenOpts {
    IOManipComplexOpts cplx_opts{};

    // setters

    /* \brief Sets std::complex<T> formatting options
     *
     * \param complex_opts Instance of qpp::IOManipComplexOpts
     * \return Reference to the current instance
     */
    IOManipEigenOpts& set_complex_opts(IOManipComplexOpts complex_opts) {
        cplx_opts = std::move(complex_opts);
        return *this;
    }
};

/*
 * \brief Formatting options for iterable objects
 * \see qpp::disp()
 */
struct IOManipRangeOpts {
    std::string sep = " ";   ///< separator
    std::string left = "[";  ///< left marking
    std::string right = "]"; ///< right marking
    realT chop = qpp::chop;  ///<  sets to zero coefficients (real or imaginary)
                             ///<  smaller than this
    // setters

    /* \brief Sets separator
     *
     * \param separator Separator
     * \return Reference to the current instance
     */
    IOManipRangeOpts& set_sep(std::string separator) {
        sep = std::move(separator);
        return *this;
    }

    /* \brief Sets left delimiter
     *
     * \param left_delimiter Left delimiter
     * \return Reference to the current instance
     */
    IOManipRangeOpts& set_left(std::string left_delimiter) {
        left = std::move(left_delimiter);
        return *this;
    }

    /* \brief Sets right delimiter
     *
     * \param right_delimiter Right delimiter
     * \return Reference to the current instance
     */
    IOManipRangeOpts& set_right(std::string right_delimiter) {
        right = std::move(right_delimiter);
        return *this;
    }

    /* \brief Sets chopping threshold
     *
     * \param chop_at Chopping threshold, sets to zero coefficients (real or
     * imaginary) smaller than this
     * \return Reference to the current instance
     */
    IOManipRangeOpts& set_chop(realT chop_at) {
        chop = chop_at;
        return *this;
    }
};

/*
 * \brief Formatting options for containers
 * \see qpp::disp()
 */
struct IOManipContainerOpts {
    std::string sep = " ";   ///< separator
    std::string left = "[";  ///< left marking
    std::string right = "]"; ///< right marking
    realT chop = qpp::chop;  ///<  sets to zero coefficients (real or imaginary)
                             ///<  smaller than this

    // setters

    /* \brief Sets separator
     *
     * \param separator Separator
     * \return Reference to the current instance
     */
    IOManipContainerOpts& set_sep(std::string separator) {
        sep = std::move(separator);
        return *this;
    }

    /* \brief Sets left delimiter
     *
     * \param left_delimiter Left delimiter
     * \return Reference to the current instance
     */
    IOManipContainerOpts& set_left(std::string left_delimiter) {
        left = std::move(left_delimiter);
        return *this;
    }

    /* \brief Sets right delimiter
     *
     * \param right_delimiter Right delimiter
     * \return Reference to the current instance
     */
    IOManipContainerOpts& set_right(std::string right_delimiter) {
        right = std::move(right_delimiter);
        return *this;
    }

    /* \brief Sets chopping threshold
     *
     * \param chop_at Chopping threshold, sets to zero coefficients (real or
     * imaginary) smaller than this
     * \return Reference to the current instance
     */
    IOManipContainerOpts& set_chop(realT chop_at) {
        chop = chop_at;
        return *this;
    }

    /*
     * \brief Conversion operator to qpp::IOManipRange
     */
    operator IOManipRangeOpts() {
        IOManipRangeOpts range_opts;

        range_opts.sep = sep;
        range_opts.left = left;
        range_opts.right = right;
        range_opts.chop = chop;
        return range_opts;
    }
};

/*
 * \brief Formatting options for C-style pointers
 * \see qpp::disp()
 */
struct IOManipPointerOpts {
    std::string sep = " ";   ///< separator
    std::string left = "[";  ///< left marking
    std::string right = "]"; ///< right marking
    realT chop = qpp::chop;  ///<  sets to zero coefficients (real or imaginary)
                             ///<  smaller than this
    // setters

    /* \brief Sets separator
     *
     * \param separator Separator
     * \return Reference to the current instance
     */
    IOManipPointerOpts& set_sep(std::string separator) {
        sep = std::move(separator);
        return *this;
    }

    /* \brief Sets left delimiter
     *
     * \param left_delimiter Left delimiter
     * \return Reference to the current instance
     */
    IOManipPointerOpts& set_left(std::string left_delimiter) {
        left = std::move(left_delimiter);
        return *this;
    }

    /* \brief Sets right delimiter
     *
     * \param right_delimiter Right delimiter
     * \return Reference to the current instance
     */
    IOManipPointerOpts& set_right(std::string right_delimiter) {
        right = std::move(right_delimiter);
        return *this;
    }

    /* \brief Sets chopping threshold
     *
     * \param chop_at Chopping threshold, sets to zero coefficients (real or
     * imaginary) smaller than this
     * \return Reference to the current instance
     */
    IOManipPointerOpts& set_chop(realT chop_at) {
        chop = chop_at;
        return *this;
    }
};

/**
 * \brief Formatting options for qpp::dirac_t objects
 * \see qpp::disp()
 */
struct IOManipDiracOpts {
    IOManipComplexOpts cplx_opts{};
    std::string plus_op = "\n";   ///< addition operator
    std::string mul_op = " * ";   ///< multiplication operator
                                  ///< absolute value
    bool amplitudes_after = true; ///< amplitudes are displayed after bra/kets
    bool discard_zeros = true;    ///< do not display zero coefficients

    // setters

    /* \brief Sets std::complex<T> formatting options
     *
     * \param complex_opts Instance of qpp::IOManipComplexOpts
     * \return Reference to the current instance
     */
    IOManipDiracOpts& set_complex_opts(IOManipComplexOpts complex_opts) {
        cplx_opts = std::move(complex_opts);
        return *this;
    }

    /* \brief Sets plus operator
     *
     * \param plus_operator Plus operator
     * \return Reference to the current instance
     */
    IOManipDiracOpts& set_plus_op(std::string plus_operator) {
        plus_op = std::move(plus_operator);
        return *this;
    }

    /* \brief Sets multiplication operator
     *
     * \param mul_operator Multiplication operator
     * \return Reference to the current instance
     */
    IOManipDiracOpts& set_mul_op(std::string mul_operator) {
        mul_op = std::move(mul_operator);
        return *this;
    }

    /* \brief Sets amplitudes after
     *
     * \param show_amplitudes_after If true, amplitudes are displayed after
     * bra/kets \return Reference to the current instance
     */

    IOManipDiracOpts& set_amplitudes_after(bool show_amplitudes_after) {
        amplitudes_after = show_amplitudes_after;
        return *this;
    }

    /* \brief Sets discarding zeros
     *
     * \param discarding_zeros If true, do not display real/imag parts of the
     * coefficients that are smaller than this->cplx_opts.chop
     * \return Reference to the current instance
     */
    IOManipDiracOpts& set_discard_zeros(bool discarding_zeros) {
        discard_zeros = discarding_zeros;
        return *this;
    }
};

namespace internal {
/**
 * \brief Maximum number of allowed qubits/qudits (subsystems)
 *
 * Used internally to allocate arrays on the stack (for performance reasons)
 */
constexpr idx maxn = 64;
} /* namespace internal */

} /* namespace qpp */

#endif /* QPP_OPTIONS_HPP_ */
