/* 
 * File:   types.h
 * Author: vlad
 *
 * Created on December 12, 2013, 10:41 PM
 */

#ifndef TYPES_H_
#define	TYPES_H_

namespace qpp
{

namespace types
{
// type aliases
using cplx=std::complex<double>;
// complex number double precision

// complex matrix, dynamic size
using cmat=Eigen::MatrixXcd;

// double matrix, dynamic size
using dmat=Eigen::MatrixXd;

// ket, dynamic size
using ket=Eigen::Matrix<cplx, Eigen::Dynamic, 1>;

// bra, dynamic size
using bra=Eigen::Matrix<cplx, 1, Eigen::Dynamic>;

// Eigen dynamic matrix
template<typename Scalar>
using DynMat=Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;

} /* namespace types */
} /* namespace qpp */

#endif	/* TYPES_H_ */

