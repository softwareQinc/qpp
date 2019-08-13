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
 * \see qpp::ILayout
 * \brief Various qudit placement layouts, all must implement the interface 
 * \a ILayout
 *
 */

#ifndef CLASSES_LAYOUTS_H_
#define CLASSES_LAYOUTS_H_

namespace qpp {

/**
 * \class qpp::ILayout
 * \brief Mandatory interface for layouts
 */
class ILayout {
  public:
    /**
     * \brief Computes the index of the point represented by \a xs in the
     * lattice coordinate system
     *
     * \param xs Vector of coordinates in the lattice coordinate system
     * \return Index of the point represented by \a xs in the lattice coordinate
     * system
     */
    virtual idx operator()(const std::vector<idx>& xs) const = 0;

    /**
     * \brief Converts index to coordinates
     *
     * \param i Index
     * \return coordinates of the point with index \a i
     */
    virtual std::vector<idx> to_coordinates(idx i) const = 0;

    /**
     * \brief Coordinate system dimensions
     *
     * \return Coordinate system dimensions
     */
    virtual std::vector<idx> get_dims() const = 0;

    /**
     * \brief Default virtual destructor
     */
    virtual ~ILayout() = default;
}; /* class ILayout */

/**
 * \class qpp::Lattice
 * \brief N-dimensional orthogonal lattice coordinate system
 */
class Lattice : public ILayout {
    std::vector<idx> dims_; ///< lattice dimensions
    /**
     * \brief Computes the product of the first \a n dimensions
     *
     * \param n Index
     * \return Product of the first \a n dimensions
     */
    idx prod_dims(idx n) const {
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
    idx compute_(const std::vector<idx>& xs, idx n) const {
        if (n == 0)
            return 0;
        else {
            // EXCEPTION CHECKS

            if (xs[n - 1] >= dims_[n - 1])
                throw exception::OutOfRange("qpp::Lattice::operator()");
            // END EXCEPTION CHECKS

            return xs[n - 1] * prod_dims(n - 1) + compute_(xs, n - 1);
        }
    }

  public:
    /**
     * \brief Constructor
     * \param dims Vector of lattice dimensions
     */
    Lattice(const std::vector<idx>& dims) : dims_{dims} {
        // EXCEPTION CHECKS

        if (dims.empty())
            throw exception::ZeroSize("qpp::Lattice::Lattice()");
        // END EXCEPTION CHECKS
    }

    /**
     * \brief Variadic constructor
     *
     * \tparam Ts Variadic type list
     * \param ds Lattice dimensions
     */
    template <class... Ts>
    Lattice(Ts... ds) : Lattice(std::vector<idx>{static_cast<idx>(ds)...}) {}

    /**
     * \brief Computes the index of the point represented by \a xs in the
     * lattice coordinate system
     *
     * \param xs Vector of coordinates in the lattice coordinate system
     * \return Index of the point represented by \a xs in the lattice coordinate
     * system
     */
    idx operator()(const std::vector<idx>& xs) const override {
        // EXCEPTION CHECKS

        if (xs.size() != dims_.size())
            throw exception::DimsNotEqual("qpp::Lattice::operator()");
        // END EXCEPTION CHECKS

        return compute_(xs, dims_.size());
    }

    /**
     * \brief Computes the index of the point represented by \a xs in the
     * lattice coordinate system

     * \tparam Ts Variadic yype list
     * \param xs Coordinates in the lattice coordinate system
     * \return Index of the point represented by \a xs in the lattice coordinate
     * system
     */
    template <class... Ts>
    idx operator()(Ts... xs) const {
        return operator()(std::vector<idx>{static_cast<idx>(xs)...});
    }

    /**
     * \brief Converts index to lattice coordinates
     *
     * \param i Index
     * \return Lattice coordinates of the point with index \a i
     */
    std::vector<idx> to_coordinates(idx i) const override {
        // EXCEPTION CHECKS

        if (i >= prod(dims_))
            throw exception::OutOfRange("qpp::Lattice::to_coordinates()");
        // END EXCEPTION CHECKS

        std::vector<idx> result = dims_;
        std::reverse(std::begin(result), std::end(result));
        result = n2multiidx(i, result);
        std::reverse(std::begin(result), std::end(result));

        return result;
    }

    /**
     * \brief Lattice dimensions
     *
     * \return Lattice dimensions
     */
    std::vector<idx> get_dims() const override { return dims_; }
}; /* class Lattice */

} /* namespace qpp */

#endif /* CLASSES_LAYOUTS_H_ */
