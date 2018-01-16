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
* \file types.h
* \brief Type aliases
*/

#ifndef TYPES_H_
#define TYPES_H_

namespace qpp {
/**
* \brief Non-negative integer index
*/
using idx = std::size_t;

/**
* \brief Big integer
*/
using bigint = long long int;

/**
* \brief Complex number in double precision
*/
using cplx = std::complex<double>;

/**
* \brief Complex (double precision) dynamic Eigen column vector
*/
using ket = Eigen::VectorXcd;

/**
* \brief Complex (double precision) dynamic Eigen row vector
*/
using bra = Eigen::RowVectorXcd;

/**
* \brief Complex (double precision) dynamic Eigen matrix
*/
using cmat = Eigen::MatrixXcd;

/**
* \brief Real (double precision) dynamic Eigen matrix
*/
using dmat = Eigen::MatrixXd;

/**
* \brief Dynamic Eigen matrix over the field specified by \a Scalar
*
* Example:
* \code
* // type of mat is Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>
* dyn_mat<float> mat(2, 3);
* \endcode
*/
template <typename Scalar> // Eigen::MatrixX_type (where type = Scalar)
using dyn_mat = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;

/**
* \brief Dynamic Eigen column vector over the field specified by \a Scalar
*
* Example:
* \code
* // type of colvect is Eigen::Matrix<float, Eigen::Dynamic, 1>
* dyn_col_vect<float> colvect(2);
* \endcode
*/
template <typename Scalar> // Eigen::VectorX_type (where type = Scalar)
using dyn_col_vect = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

/**
* \brief Dynamic Eigen row vector over the field specified by \a Scalar
*
* Example:
* \code
* // type of rowvect is Eigen::Matrix<float, 1, Eigen::Dynamic>
* dyn_row_vect<float> rowvect(3);
* \endcode
*/
template <typename Scalar> // Eigen::RowVectorX_type (where type = Scalar)
using dyn_row_vect = Eigen::Matrix<Scalar, 1, Eigen::Dynamic>;

} /* namespace qpp */

#endif /* TYPES_H_ */
