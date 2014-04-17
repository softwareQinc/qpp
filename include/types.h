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

// complex matrix, dynamic size
typedef Eigen::MatrixXcd cmat;

// double matrix, dynamic size
typedef Eigen::MatrixXd dmat;

// ket, dynamic size
typedef Eigen::Matrix<cplx, Eigen::Dynamic, 1> ket;

// bra, dynamic size
typedef Eigen::Matrix<cplx, 1, Eigen::Dynamic> bra;

// Eigen dynamic matrix
template<typename Scalar>
using DynMat=Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;

} /* namespace types */
} /* namespace qpp */

#endif	/* TYPES_H_ */

