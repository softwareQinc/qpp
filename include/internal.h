/*
 * internal.h
 *
 *  Created on: Mar 24, 2014
 *      Author: vlad
 */

#ifndef INTERNAL_H_
#define INTERNAL_H_

// internal functions, do not modify

namespace qpp{

// integer index to multi-index, use C-style array for speed
void _n2multiidx(const size_t n, const size_t numdims, const size_t *dims,
		size_t *result);

// multi index to integer index, use C-style array for speed
size_t _multiidx2n(const size_t *midx, const size_t numdims,
		const size_t *dims);

// used inside the #pragma omp parallel for in syspermute
void _syspermute_worker(const size_t numdims, const size_t *cdims,
		const size_t *cperm, const size_t i, const size_t j, size_t &iperm,
		size_t &jperm, const types::cmat &A, types::cmat &result);

// used inside the #pragma omp parallel for in ptranspose
void _ptranspose_worker(const size_t numdims, const size_t numsubsys, const size_t *cdims,
		const size_t *csubsys, const size_t i, const size_t j, size_t &iperm,
		size_t &jperm, const types::cmat &A, types::cmat &result);

}


#endif /* INTERNAL_H_ */
