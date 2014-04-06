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
#include "constants.h"

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

// converts Expression (e.g. from an Eigen::MatrixBase<Expression>) to
// Eigen dynamic matrix Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>
template<typename Expression>
using Expression2DynMat=
Eigen::Matrix<typename Expression::Scalar, Eigen::Dynamic, Eigen::Dynamic>;

// Eigen dynamic matrix
template<typename Scalar>
using DynMat=Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;

int myfunc(int a, int b)
{
	//auto i=gt::X;

	return a + b;
}

}
}

#endif	/* TYPES_H_ */

