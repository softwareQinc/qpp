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
 * \file classes/circuits/layouts.h
 * \brief Various physical qudit layouts, all must implement
 * idx operator(idx, idx, ...)
 */

#ifndef CLASSES_CIRCUITS_LAYOUTS_H_
#define CLASSES_CIRCUITS_LAYOUTS_H_

namespace qpp {
/**
 * \class qpp::GridView
 * \brief Square grid qudit layout
 */
class GridView {
  protected:
    idx dim_x_, dim_y_;

  public:
    /**
     * \brief Constructs a GridView (planar orthogonal grid)
     * \param dim_x Number of elements on the X axis
     * \param dim_y Number of elements on the Y axis
     */
    GridView(idx dim_x, idx dim_y) : dim_x_{dim_x}, dim_y_{dim_y} {}

    /**
     * \brief Converts GridView coordinates to linear coordinates
     *
     * \param x X coordinate
     * \param y Y coordinate
     * \return Linear coordinate of the (x, y) pair
     */
    idx operator()(idx x, idx y) const {
        // EXCEPTION CHECKS

        if (x >= dim_x_ || y >= dim_y_)
            throw exception::OutOfRange("qpp::GridView::operator()");
        // END EXCEPTION CHECKS

        return dim_x_ * y + x;
    }
};

/**
 * \class qpp::3DGridView
 * \brief 3D-grid qudit layout
 */
class GridView3D : private GridView {
    idx dim_z_;

  public:
    /**
     * \brief Constructs a GridView3D (3 dimensional orthogonal grid)
     * \param dim_x Number of elements on the X axis
     * \param dim_y Number of elements on the Y axis
     * \param dim_z Number of elements on the Z axis
     */
    GridView3D(idx dim_x, idx dim_y, idx dim_z)
        : GridView(dim_x, dim_y), dim_z_{dim_z} {}

    /**
     * \brief Converts GridView coordinates to linear coordinates
     *
     * \param x X coordinate
     * \param y Y coordinate
     * \param z Z coordinate
     * \return Linear coordinate of the (x, y, z) pair
     */
    idx operator()(idx x, idx y, idx z) const {
        // EXCEPTION CHECKS

        if (x >= dim_x_ || y >= dim_y_)
            throw exception::OutOfRange("qpp::GridView3D::operator()");
        // END EXCEPTION CHECKS

        return dim_x_ * dim_y_ * z + GridView::operator()(x, y);
    }
};
} /* namespace qpp */

#endif /* CLASSES_CIRCUITS_LAYOUTS_H_ */