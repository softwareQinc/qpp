/*
 * entropy.h
 *
 *  Created on: Mar 27, 2014
 *      Author: vlad
 */

#ifndef ENTROPY_H_
#define ENTROPY_H_

#include "types.h"
#include "util.h"
#include <cmath>

// entropy functions

namespace qpp
{

// von-Neumann entropy
template<typename Derived>
double entropyS(const Eigen::MatrixBase<Derived> & A)
{
	// get the eigenvalues
	Eigen::MatrixXcd ev=evals(A);
	double result=0;
	for (size_t i=0;i<ev.rows();i++)
		if(std::abs(ev(i))!=0)
			result-=std::abs(ev(i))*std::log2(std::abs(ev(i)));
	return result;
}

// Shannon entropy
//inline double entropyH(std::vector<>){}
}

#endif /* ENTROPY_H_ */
