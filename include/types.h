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

/**
 * @brief Complex number in double precision
 */
using cplx=std::complex<double>;

/**
 * @brief Complex (double precision) dynamic Eigen matrix
 */
using cmat=Eigen::MatrixXcd;

/**
 * @brief Real (double precision) dynamic Eigen matrix
 */
using dmat=Eigen::MatrixXd;

/**
 * @brief Complex (double precision) dynamic Eigen column matrix
 */
using ket=Eigen::Matrix<cplx, Eigen::Dynamic, 1>;

/**
 * @brief Complex (double precision) dynamic Eigen row matrix
 */
using bra=Eigen::Matrix<cplx, 1, Eigen::Dynamic>;

/**
 * @brief Dynamic Eigen matrix over the field specified by @a Scalar
 *
 * Example:
 * @code auto mat = DynMat<float>(2,3); // type of mat is Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> @endcode
 */
template<typename Scalar>
using DynMat=Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;

} /* namespace types */
} /* namespace qpp */

#endif	/* TYPES_H_ */

