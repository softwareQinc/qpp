/* 
 * File:   types.h
 * Author: vlad
 *
 * Created on December 12, 2013, 10:41 PM
 */

#ifndef INCLUDE_TYPES_H_
#define	INCLUDE_TYPES_H_

namespace qpp
{

/**
 * \brief Complex number in double precision
 */
using cplx = std::complex<double>;

/**
 * \brief Dynamic Eigen matrix over the field specified by \a Scalar
 *
 * Example:
 * \code auto mat = DynMat<float>(2,3); // type of mat is Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> \endcode
 */
template<typename Scalar> // Eigen::MatrixX<type> (where type = Scalar)
using DynMat = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;

/**
 * \brief Dynamic Eigen column vector over the field specified by \a Scalar
 *
 * Example:
 * \code auto colvect = DynColVect<float>(2); // type of colvect is Eigen::Matrix<float, Eigen::Dynamic, 1> \endcode
 */
template<typename Scalar> // Eigen::VectorX<type> (where type = Scalar)
using DynColVect = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

/**
 * \brief Dynamic Eigen row vector over the field specified by \a Scalar
 *
 * Example:
 * \code auto rowvect = DynRowVect<float>(3); // type of rowvect is Eigen::Matrix<float, 1, Eigen::Dynamic> \endcode
 */
template<typename Scalar> // Eigen::RowVectorX<type> (where type = Scalar)
using DynRowVect = Eigen::Matrix<Scalar, 1, Eigen::Dynamic>;

/**
 * \brief Complex (double precision) dynamic Eigen column vector
 */
using ket = DynColVect<cplx>;
// Eigen::VectorXcd

/**
 * \brief Complex (double precision) dynamic Eigen row vector
 */
using bra = DynRowVect<cplx>;
// Eigen::RowVectorXcd

/**
 * \brief Complex (double precision) dynamic Eigen matrix
 */
using cmat = DynMat<cplx>;
// Eigen::MatrixXcd;

/**
 * \brief Real (double precision) dynamic Eigen matrix
 */
using dmat = DynMat<double>;
// Eigen::MatrixXd

} /* namespace qpp */

#endif	/* INCLUDE_TYPES_H_ */

