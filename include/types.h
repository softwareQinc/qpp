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
typedef Eigen::Matrix2cd cmat2; // 2 x 2
typedef Eigen::Matrix3cd cmat3; // 3 x 3
typedef Eigen::Matrix4cd cmat4; // 4 x 4
typedef Eigen::MatrixXcd cmat; // dynamic-size
//class cmat:public Eigen::MatrixXcd
//{
//public:
//	cmat(size_t m, size_t n):Eigen::MatrixXcd(m,n){};
//};

// complex ket vector
typedef Eigen::Vector2cd ket2; // 2 x 1
typedef Eigen::Vector3cd ket3; // 3 x 1
typedef Eigen::Vector4cd ket4; // 4 x 1
typedef Eigen::VectorXcd ket; // dynamic-size ket

// complex bra vector
typedef Eigen::RowVector2cd bra2; // 1 x 2
typedef Eigen::RowVector3cd bra3; // 1 x 3
typedef Eigen::RowVector4cd bra4; // 1 x 4
typedef Eigen::RowVectorXcd bra; // dynamic-size bra

typedef Eigen::VectorXcd cvect; // dynamic-size complex vector
typedef Eigen::VectorXi ivect; // dynamic-size integer vector

}
}

#endif	/* TYPES_H_ */

