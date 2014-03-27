/* 
 * File:   types.h
 * Author: vlad
 *
 * Created on December 12, 2013, 10:41 PM
 */

#ifndef TYPES_H_
#define	TYPES_H_

#include <complex>
#include <Eigen/Dense>

// wrappers around various useful types

namespace qpp
{

namespace types
{
// typedefs
typedef std::complex<double> cplx; // complex number double precision

// complex matrices
typedef Eigen::MatrixXcd cmat; // dynamic-size

// complex ket vector
typedef Eigen::VectorXcd ket; // dynamic-size ket

// TODO: make functions like kron work with vector also
// TODO: templatize everything to MatrixBase!

// complex bra vector
typedef Eigen::RowVectorXcd bra; // dynamic-size bra

}
}

#endif	/* TYPES_H_ */

