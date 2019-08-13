/*
 * This file is part of Quantum++.
 *
 * MIT License
 *
 * Copyright (c) 2013 - 2019 Vlad Gheorghiu (vgheorgh@gmail.com)
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
 * \file classes/layouts.h
 * \brief Various physical qudit layouts, all must implement
 * template <class... Ts> idx operator()(Ts... xs)
 */

#ifndef CLASSES_LAYOUTS_H_
#define CLASSES_LAYOUTS_H_

namespace qpp {
/**
 * \class Lattice
 * \brief N-dimensional orthogonal lattice coordinate system
 */
class Lattice {
    std::vector<idx> dims_; ///< lattice dimensions
    /**
     * \brief Computes the product of the first \a n dimensions
     *
     * \param n Index
     * \return Product of the first \a n dimensions
     */
    idx prod_dims(idx n) {
        return prod(std::begin(dims_), std::next(std::begin(dims_), n));
    }

    /**
     * \brief Computes the index of the point represented by \a xs in the
     * lattice coordinate system
     *
     * \param xs Coordinates in the lattice coordinate system
     * \param n Recursion index
     * \return Index of the point represented by \a xs in the lattice coordinate
     * system
     */
    idx compute_(const std::vector<idx>& xs, idx n) {
        if (n == 0)
            return 0;
        else {
            if (xs[n - 1] >= dims_[n - 1])
                throw exception::OutOfRange("qpp::Lattice::operator()()");
            return xs[n - 1] * prod_dims(n - 1) + compute_(xs, n - 1);
        }
    }

    /**
     * \brief Constructor (private)
     * \param dims Vector of lattice dimensions
     */
    Lattice(const std::vector<idx>& dims) : dims_{dims} {
        if (dims.empty())
            throw exception::ZeroSize("qpp::Lattice::Lattice()");
    }

  public:
    /**
     * \brief Constructor
     *
     * \tparam Ts Type list
     * \param ds Lattice dimensions
     */
    template <class... Ts>
    Lattice(Ts... ds) : Lattice(std::vector<idx>{static_cast<idx>(ds)...}) {}

    /**
     * \brief Computes the index of the point represented by \a xs in the
     * lattice coordinate system

     * \tparam Ts Type list
     * \param xs Coordinates in the lattice coordinate system
     * \return Index of the point represented by \a xs in the lattice coordinate
     * system
     */
    template <class... Ts>
    idx operator()(Ts... xs) {
        if (sizeof...(xs) != dims_.size())
            throw exception::DimsNotEqual("qpp::Lattice::operator()()");
        return compute_(std::vector<idx>{static_cast<idx>(xs)...},
                        dims_.size());
    }
};

} /* namespace qpp */

#endif /* CLASSES_LAYOUTS_H_ */