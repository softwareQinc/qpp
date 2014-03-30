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

// General Eigen expression, use this for function parameters that must work
// with expressions also
template<typename Derived> using EigenExpression=Eigen::MatrixBase<Derived>;

// Use TemplatedEigenMatrix<Derived> as return for functions that take
// EigenExpression<Derived> as inputs, so the function work with
// Eigen expressions, as every EigenExpression<Derived> can be implicitly
// converted to TemplatedEigenMatrix<Derived>
template<typename MatrixType> using TemplatedEigenMatrix=
Eigen::Matrix<typename MatrixType::Scalar, Eigen::Dynamic, Eigen::Dynamic>;

}
}

#endif	/* TYPES_H_ */

