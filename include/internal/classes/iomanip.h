/*
 * This file is part of Quantum++.
 *
 * MIT License
 *
 * Copyright (c) 2013 - 2018 Vlad Gheorghiu (vgheorgh@gmail.com)
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
* \file internal/classes/iomanip.h
* \brief Input/output manipulators
*/

#ifndef INTERNAL_CLASSES_IOMANIP_H_
#define INTERNAL_CLASSES_IOMANIP_H_

namespace qpp {
namespace internal {
// ostream manipulators for nice formatting of
// Eigen matrices and STL/C-style containers/vectors

template <typename InputIterator>
class IOManipRange : public IDisplay {
    InputIterator first_, last_;
    std::string separator_, start_, end_;

  public:
    explicit IOManipRange(InputIterator first, InputIterator last,
                          const std::string& separator,
                          const std::string& start = "[",
                          const std::string& end = "]")
        : first_{first}, last_{last}, separator_{separator}, start_{start},
          end_{end} {}

    // to silence -Weffc++ warnings for classes that have pointer members
    // (whenever we have a pointer instantiation,
    // i.e. iterator is a raw pointer)
    IOManipRange(const IOManipRange&) = default;

    IOManipRange& operator=(const IOManipRange&) = default;

  private:
    std::ostream& display(std::ostream& os) const override {
        os << start_;

        bool first = true;
        for (InputIterator it = first_; it != last_; ++it) {
            if (!first)
                os << separator_;
            first = false;
            os << *it;
        }
        os << end_;

        return os;
    }
}; // class IOManipRange

template <typename PointerType>
class IOManipPointer : public IDisplay {
    const PointerType* p_;
    idx N_;
    std::string separator_, start_, end_;

  public:
    explicit IOManipPointer(const PointerType* p, idx N,
                            const std::string& separator,
                            const std::string& start = "[",
                            const std::string& end = "]")
        : p_{p}, N_{N}, separator_{separator}, start_{start}, end_{end} {}

    // to silence -Weffc++ warnings for classes that have pointer members
    IOManipPointer(const IOManipPointer&) = default;

    IOManipPointer& operator=(const IOManipPointer&) = default;

  private:
    std::ostream& display(std::ostream& os) const override {
        os << start_;

        for (idx i = 0; i < N_ - 1; ++i)
            os << p_[i] << separator_;
        if (N_ > 0)
            os << p_[N_ - 1];

        os << end_;

        return os;
    }
}; // class IOManipPointer

// silence g++4.8.x bogus warning -Wnon-virtual-dtor for
// qpp::internal::Display_impl_ when class qpp::internal::IOManipEigen
// privately inherits from it
#if ((__GNUC__ == 4) && (__GNUC_MINOR__ == 8) && !__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Weffc++"
#endif
class IOManipEigen : public IDisplay, private Display_Impl_ {
#if ((__GNUC__ == 4) && (__GNUC_MINOR__ == 8) && !__clang__)
#pragma GCC diagnostic pop
#endif
    cmat A_;
    double chop_;

  public:
    // Eigen matrices
    template <typename Derived>
    explicit IOManipEigen(const Eigen::MatrixBase<Derived>& A,
                          double chop = qpp::chop)
        : A_{A.template cast<cplx>()}, // copy, so we can bind rvalues safely
          chop_{chop} {}

    // Complex numbers
    explicit IOManipEigen(const cplx z, double chop = qpp::chop)
        : A_{cmat::Zero(1, 1)}, chop_{chop} {
        // put the complex number inside an Eigen matrix
        A_(0, 0) = z;
    }

  private:
    std::ostream& display(std::ostream& os) const override {
        return display_impl_(A_, os, chop_);
    }
}; // class IOManipEigen

} /* namespace internal */
} /* namespace qpp */

#endif /* INTERNAL_CLASSES_IOMANIP_H_ */
