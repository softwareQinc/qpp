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
 * \file classes/idisplay.hpp
 * \brief Textual representation interface and very basic JSON serialization
 * interface
 */

#ifndef QPP_CLASSES_IDISPLAY_HPP_
#define QPP_CLASSES_IDISPLAY_HPP_

#include <cmath>
#include <iomanip>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>

#include "qpp/options.hpp"
#include "qpp/types.hpp"

namespace qpp {
/**
 * \class qpp::IDisplay
 * \brief Abstract class (interface) that mandates the definition of
 * virtual std::ostream& display(std::ostream& os) const
 *
 * This class defines
 * friend std::ostream& operator<<(std::ostream& os, const qpp::IDisplay& rhs).
 * The latter delegates the work to the pure private virtual function
 * qpp::IDisplay::display() which has to be overridden by all derived classes.
 */
class IDisplay {
    /**
     * \brief Must be overridden by all derived classes
     *
     * The actual stream extraction processing is performed by the overridden
     * member function in the derived class. This function is automatically
     * invoked by
     * friend std::ostream& operator<<(std::ostream& os, const IDisplay& rhs).
     */
    virtual std::ostream& display(std::ostream& os) const = 0;

  public:
    /**
     * \brief Default virtual destructor
     */
    virtual ~IDisplay() = default;

    /**
     * \brief Overloads the extraction operator
     *
     * Delegates the work to the virtual function qpp::IDisplay::display()
     */
    friend std::ostream& operator<<(std::ostream& os, const IDisplay& rhs) {
        return rhs.display(os);
    }
}; /* class IDisplay */

namespace internal {
// implementation details for pretty formatting
struct Display_Impl_ {
    template <typename T>
    // T must support rows(), cols(), operator()(idx, idx) const
    std::ostream& display_impl_(const T& A, std::ostream& os,
                                IOManipEigenOpts opts) const {
        std::ostringstream ostr;
        ostr.copyfmt(os); // copy os' state

        std::vector<std::string> vstr;
        std::string str;

        for (idx i = 0; i < static_cast<idx>(A.rows()); ++i) {
            for (idx j = 0; j < static_cast<idx>(A.cols()); ++j) {
                str.clear(); // clear the temporary string
                ostr.clear();
                ostr.str(std::string{}); // clear the ostringstream

                // convert to complex
                realT re = static_cast<cplx>(A(i, j)).real();
                realT im = static_cast<cplx>(A(i, j)).imag();

                // zero
                if (std::abs(re) < opts.chop && std::abs(im) < opts.chop) {
                    ostr << "0"; // otherwise, segfault on destruction
                    // if using only vstr.emplace_back("0 ");
                    // bug in MATLAB libmx
                    vstr.emplace_back(ostr.str());
                }
                // pure imag
                else if (std::abs(re) < opts.chop) {
                    ostr << im;
                    vstr.emplace_back(ostr.str() + "i");
                }
                // real
                else if (std::abs(im) < opts.chop) {
                    ostr << re;
                    vstr.emplace_back(ostr.str());
                }
                // full complex
                else {
                    ostr << re;
                    str = ostr.str();

                    str += (im > 0 ? opts.plus_op : opts.minus_op);
                    ostr.clear();
                    ostr.str(std::string()); // clear
                    ostr << std::abs(im);
                    str += ostr.str();
                    str += "i";
                    vstr.emplace_back(str);
                }
            }
        }

        // determine the maximum lenght of the entries in each column
        std::vector<idx> maxlengthcols(A.cols(), 0);

        for (idx i = 0; i < static_cast<idx>(A.rows()); ++i) {
            for (idx j = 0; j < static_cast<idx>(A.cols()); ++j) {
                if (static_cast<idx>(vstr[i * A.cols() + j].size()) >
                    maxlengthcols[j]) {
                    maxlengthcols[j] = vstr[i * A.cols() + j].size();
                }
            }
        }

        // finally, display it!
        for (idx i = 0; i < static_cast<idx>(A.rows()); ++i) {
            os << std::setw(static_cast<int>(maxlengthcols[0])) << std::right
               << vstr[i * A.cols()]; // display first column
            // then the rest
            idx spacer = 2;
            for (idx j = 1; j < static_cast<idx>(A.cols()); ++j) {
                os << std::setw(static_cast<int>(maxlengthcols[j] + spacer))
                   << std::right << vstr[i * A.cols() + j];
            }

            if (i < static_cast<idx>(A.rows()) - 1) {
                os << '\n';
            }
        }

        return os;
    }
};
} /* namespace internal */

/**
 * \class qpp::IJSON
 * \brief Abstract class (interface) that mandates the definition of
 * very basic JSON serialization support
 */
class IJSON {
  public:
    /**
     * \brief Default virtual destructor
     */
    virtual ~IJSON() = default;

    /**
     * \brief JSON representation of the derived instance, must be overridden by
     * all derived classes
     *
     * \param enclosed_in_curly_brackets If true, encloses the result in curly
     * brackets
     */
    virtual std::string
    to_JSON(bool enclosed_in_curly_brackets = true) const = 0;

}; /* class IJSON */

} /* namespace qpp */

#endif /* QPP_CLASSES_IDISPLAY_HPP_ */
