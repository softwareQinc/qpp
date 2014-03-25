/*
 * internal.h
 *
 *  Created on: Mar 24, 2014
 *      Author: vlad
 */

#ifndef INTERNAL_H_
#define INTERNAL_H_

#include <vector>
#include "types.h"

// internal functions, do not modify
namespace qpp
{
namespace internal
{

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
void _ptranspose_worker(const size_t* midxrow, const size_t numdims,
		const size_t numsubsys, const size_t *cdims, const size_t *csubsys,
		const size_t i, const size_t j, size_t &iperm, size_t &jperm,
		const types::cmat &A, types::cmat &result);

// check square matrix
bool _check_square_mat(const types::cmat &A);

// check that dims is a valid dimension vector
bool _check_dims(const std::vector<size_t> &dims);

// check that all elements in dims equal to dim
bool _check_eq_dims(const std::vector<size_t> &dims, size_t dim);

// check that dims match the dimension of the matrix
bool _check_dims_match_mat(const std::vector<size_t> & dims, const types::cmat &A);

// check that subsys is valid with respect to dims
bool _check_subsys(const std::vector<size_t> & subsys, const std::vector<size_t> & dims);

// check that the permutation is valid with respect to dims
bool _check_perm(const std::vector<size_t> & perm, const std::vector<size_t> & dims);


}
}

#endif /* INTERNAL_H_ */
