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

namespace qpp
{
namespace types
{
// typedefs
typedef std::complex<double> cplx; // complex number double precision

// complex matrices
typedef Eigen::MatrixXcd cmat; // dynamic-size

}
}

#endif	/* TYPES_H_ */

