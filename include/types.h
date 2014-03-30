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

// Use to match expressions (i.e. A+B)
template<typename MatrixType>
using EigenExpression=Eigen::MatrixBase<MatrixType>;

// Use to match Eigen matrices regardless of scalar (i.e. cmat, dmat etc)
// Use in general to convert from input of type EigenExpression
template<typename MatrixType>
using TemplatedEigenMatrix=
Eigen::Matrix<typename MatrixType::Scalar, Eigen::Dynamic, Eigen::Dynamic>;

// Use to match specific Scalar types (double/cplx etc)
// Can not extract Scalar from an expression
template<typename Scalar>
using ScalarEigenMatrix=Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;

}
}

#endif	/* TYPES_H_ */

