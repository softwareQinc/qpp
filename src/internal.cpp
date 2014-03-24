/*
 * internal.cpp
 *
 *  Created on: Mar 24, 2014
 *      Author: vlad
 */

#include "internal.h"

namespace qpp
{
namespace internal
{

// integer index to multi-index, use C-style array for speed
void _n2multiidx(const size_t n, const size_t numdims, const size_t *dims,
		size_t *result)
{
	size_t maxn = 1;
	for (size_t i = 0; i < numdims; i++)
		maxn *= dims[i];
	if (n > maxn - 1)
		throw std::runtime_error(
				"_n2multiidx: Number too large, out of bounds!");

	size_t _n = n;
	for (size_t i = 0; i < numdims; i++)
	{
		result[numdims - i - 1] = _n % static_cast<int>(dims[numdims - i - 1]);
		_n = _n / static_cast<int>(dims[numdims - i - 1]);
	}
}

// multi index to integer index, use C-style array for speed
size_t _multiidx2n(const size_t *midx, const size_t numdims, const size_t *dims)
{
	for (size_t i = 0; i < numdims; i++)
		if (midx[i] >= dims[i])
			throw std::runtime_error(
					"_multiidx2n: Sub-index exceeds corresponding dimension!");

	size_t *part_prod = new size_t[numdims];

	part_prod[numdims - 1] = 1;
	for (size_t j = 1; j < numdims; j++)
		part_prod[numdims - j - 1] = part_prod[numdims - j] * dims[numdims - j];

	size_t result = 0;
	for (size_t i = 0; i < numdims; i++)
		result += midx[i] * part_prod[i];

	delete[] part_prod;

	return result;
}

// used inside the #pragma omp parallel for in syspermute
void _syspermute_worker(const size_t numdims, const size_t *cdims,
		const size_t *cperm, const size_t i, const size_t j, size_t &iperm,
		size_t &jperm, const types::cmat &A, types::cmat &result)
{
	size_t *midxrow = new size_t[numdims];
	size_t *midxcol = new size_t[numdims];
	size_t *midxrowtmp = new size_t[numdims];
	size_t *midxcoltmp = new size_t[numdims];
	size_t *permdims = new size_t[numdims];

	for (size_t i = 0; i < numdims; i++)
		permdims[i] = cdims[cperm[i]]; // permuted dimensions

	// compute the row and col multi-indexes
	_n2multiidx(i, numdims, cdims, midxrow);
	_n2multiidx(j, numdims, cdims, midxcol);

	for (size_t k = 0; k < numdims; k++)
	{
		midxrowtmp[k] = midxrow[cperm[k]]; // permuted multi-indexes
		midxcoltmp[k] = midxcol[cperm[k]]; // permuted multi-indexes
	}

	// move back to integer indexes
	iperm = _multiidx2n(midxrowtmp, numdims, permdims);
	jperm = _multiidx2n(midxcoltmp, numdims, permdims);
	result(iperm, jperm) = A(i, j);

	delete[] midxrow;
	delete[] midxcol;
	delete[] midxrowtmp;
	delete[] midxcoltmp;
	delete[] permdims;
}

// used inside the #pragma omp parallel for in ptranspose
void _ptranspose_worker(const size_t numdims, const size_t numsubsys,
		const size_t *cdims, const size_t *csubsys, const size_t i,
		const size_t j, size_t &iperm, size_t &jperm, const types::cmat &A,
		types::cmat &result)
{

	size_t *midxrow = new size_t[numdims];
	size_t *midxcol = new size_t[numdims];

	// compute the row and col multi-indexes
	_n2multiidx(i, numdims, cdims, midxrow);
	_n2multiidx(j, numdims, cdims, midxcol);

	for (size_t k = 0; k < numsubsys; k++)
		std::swap(midxrow[csubsys[k]], midxcol[csubsys[k]]);

	// move back to integer indexes
	iperm = _multiidx2n(midxrow, numdims, cdims);
	jperm = _multiidx2n(midxcol, numdims, cdims);
	result(iperm, jperm) = A(i, j);

	delete[] midxrow;
	delete[] midxcol;
}
}

}

