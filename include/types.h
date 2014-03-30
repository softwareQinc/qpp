/* 
 * File:   types.h
 * Author: vlad
 *
 * Created on December 12, 2013, 10:41 PM
 */

#ifndef TYPES_H_
#define	TYPES_H_

#include <Eigen/Dense>
#include <complex>

namespace qpp
{

namespace types
{
// typedefs
typedef std::complex<double> cplx; // complex number double precision

// complex matrix
typedef Eigen::MatrixXcd cmat; // dynamic-size

// double matrix
typedef Eigen::MatrixXd dmat; // dynamic-size

// float matrix
typedef Eigen::MatrixXf fmat; // dynamic-size

// integer matrix
typedef Eigen::MatrixXi imat; // dynamic-size

// Eigen templated matrices

// Matches expression (i.e. A+B) regardless of matrix type
// Use to manipulate the expression as a matrix
// Deduce MatrixType automatically in an expression
template<typename MatrixType>
using EigenExpression=Eigen::MatrixBase<MatrixType>;

// Matches matrix regardless of scalar (i.e. cmat, dmat etc)
// Use it as return type after expression manipulation
// Can not automatically deduce MatrixType (do not use it as function parameter)
template<typename MatrixType>
using TemplatedEigenMatrix=
Eigen::Matrix<typename MatrixType::Scalar, Eigen::Dynamic, Eigen::Dynamic>;

// Matches specific Scalar types (double/cplx etc)
// Can not automatically deduce Scalar from an expression
// Use it for strong type checking
template<typename Scalar>
using ScalarEigenMatrix=Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;

}
}

#endif	/* TYPES_H_ */

