/*
 * Quantum++
 *
 * Copyright (c) 2013 - 2014 Vlad Gheorghiu (vgheorgh@gmail.com)
 *
 * This file is part of Quantum++.
 *
 * Quantum++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Quantum++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Quantum++.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef INCLUDE_INTERNAL_FUNCTIONS_H_
#define INCLUDE_INTERNAL_FUNCTIONS_H_

// internal functions, do not use/modify
/**
* \namespace qpp::internal Internal implementation details,
* do not use/modify these functions/classes
*/
namespace qpp
{
namespace internal
{

// integer index to multi-index, use C-style array for speed
// standard lexicographical order, e.g. 00, 01, 10, 11
inline void _n2multiidx(std::size_t n, std::size_t numdims,
        const std::size_t* dims, std::size_t* result)
{
    // no error checks to improve speed
    for (std::size_t i = 0; i < numdims; ++i)
    {
        result[numdims - i - 1] = n % (dims[numdims - i - 1]);
        n /= (dims[numdims - i - 1]);
    }
}

// multi-index to integer index, use C-style array for speed,
// standard lexicographical order, e.g. 00->0, 01->1, 10->2, 11->3
inline std::size_t _multiidx2n(const std::size_t* midx, std::size_t numdims,
        const std::size_t* dims)
{
    // no error checks to improve speed

    // Static allocation for speed!
    // double the size for matrices reshaped as vectors
    std::size_t part_prod[2 * maxn];

    std::size_t result = 0;
    part_prod[numdims - 1] = 1;
    for (std::size_t i = 1; i < numdims; ++i)
    {
        part_prod[numdims - i - 1] =
                part_prod[numdims - i] * dims[numdims - i];
        result += midx[numdims - i - 1] * part_prod[numdims - i - 1];
    }

    return result + midx[numdims - 1];
}

// check square matrix
template<typename Derived>
bool _check_square_mat(const Eigen::MatrixBase<Derived>& A)
{
    const DynMat<typename Derived::Scalar>& rA = A;

    if (rA.rows() != rA.cols())
        return false;

    return true;
}

// check whether input is a vector or not
template<typename Derived>
bool _check_vector(const Eigen::MatrixBase<Derived>& A)
{
    const DynMat<typename Derived::Scalar>& rA = A;

    if (rA.rows() != 1 && rA.cols() != 1)
        return false;

    return true;
}

// check whether input is a row vector or not
template<typename Derived>
bool _check_row_vector(const Eigen::MatrixBase<Derived>& A)
{
    const DynMat<typename Derived::Scalar>& rA = A;

    if (rA.rows() != 1)
        return false;

    return true;
}

// check whether input is a column vector or not
template<typename Derived>
bool _check_col_vector(const Eigen::MatrixBase<Derived>& A)
{
    const DynMat<typename Derived::Scalar>& rA = A;

    if (rA.cols() != 1)
        return false;

    return true;
}

// check non-zero size of object that supports size() function
template<typename T>
bool _check_nonzero_size(const T& x)
{
    if (x.size() == 0)
        return false;

    return true;
}

// check that dims is a valid dimension vector
bool _check_dims(const std::vector<std::size_t>& dims)
{
    if (dims.size() == 0)
        return false;

    if (std::find_if(std::begin(dims), std::end(dims),
            [dims](std::size_t i) -> bool
            {
                if (i == 0) return true;
                else return false;
            }) != std::end(dims))
        return false;

    return true;
}

// check that valid dims match the dimensions
// of valid (non-zero sized) quare matrix
template<typename Derived>
bool _check_dims_match_mat(const std::vector<std::size_t>& dims,
        const Eigen::MatrixBase<Derived>& A)
{
    const DynMat<typename Derived::Scalar>& rA = A;

    std::size_t proddim = 1;
    for (std::size_t i : dims)
        proddim *= i;
    if (proddim != static_cast<std::size_t>(rA.rows()))
        return false;

    return true;
}

// check that valid dims match the dimensions of valid column vector
template<typename Derived>
bool _check_dims_match_cvect(const std::vector<std::size_t>& dims,
        const Eigen::MatrixBase<Derived>& V)
{
    const DynMat<typename Derived::Scalar>& rV = V;

    std::size_t proddim = 1;
    for (std::size_t i : dims)
        proddim *= i;
    if (proddim != static_cast<std::size_t>(rV.rows()))
        return false;

    return true;
}

// check that valid dims match the dimensions of valid row vector
template<typename Derived>
bool _check_dims_match_rvect(const std::vector<std::size_t>& dims,
        const Eigen::MatrixBase<Derived>& V)
{
    const DynMat<typename Derived::Scalar>& rV = V;

    std::size_t proddim = 1;
    for (std::size_t i : dims)
        proddim *= i;
    if (proddim != static_cast<std::size_t>(rV.cols()))
        return false;

    return true;
}

// check that all elements in valid dims equal to dim
bool _check_eq_dims(const std::vector<std::size_t>& dims, std::size_t dim)
{
    for (std::size_t i : dims)
        if (i != dim)
            return false;

    return true;
}

// check that subsys is valid with respect to valid dims
bool _check_subsys_match_dims(const std::vector<std::size_t>& subsys,
        const std::vector<std::size_t>& dims)
{
    //	// check non-zero sized subsystems
    //	if (subsys.size() == 0)
    //		return false;

    // check valid number of subsystems
    if (subsys.size() > dims.size())
        return false;

    // sort the subsystems
    std::vector<std::size_t> subsyssort = subsys;
    std::sort(std::begin(subsyssort), std::end(subsyssort));

    // check duplicates
    if (std::unique(std::begin(subsyssort), std::end(subsyssort))
            != std::end(subsyssort))
        return false;

    // check range of subsystems
    if (std::find_if(std::begin(subsyssort), std::end(subsyssort),
            [dims](std::size_t i) -> bool
            {
                if (i > dims.size() - 1) return true;
                else return false;
            }) != std::end(subsyssort))
        return false;

    return true;
}

// check valid permutation
bool _check_perm(const std::vector<std::size_t>& perm)
{
    if (perm.size() == 0)
        return false;

    std::vector<std::size_t> ordered(perm.size());
    std::iota(std::begin(ordered), std::end(ordered), 0);

    if (std::is_permutation(std::begin(ordered), std::end(ordered),
            std::begin(perm)))
        return true;
    else
        return false;
}

// Kronecker product of 2 matrices, preserve return type
// internal function for the variadic template function wrapper kron()
template<typename Derived1, typename Derived2>
DynMat<typename Derived1::Scalar> _kron2(const Eigen::MatrixBase<Derived1>& A,
        const Eigen::MatrixBase<Derived2>& B)
{
    const DynMat<typename Derived1::Scalar>& rA = A;
    const DynMat<typename Derived2::Scalar>& rB = B;

    // EXCEPTION CHECKS

    // check types
    if (!std::is_same<typename Derived1::Scalar,
            typename Derived2::Scalar>::value)
        throw Exception("qpp::kron()", Exception::Type::TYPE_MISMATCH);

    // check zero-size
    if (!internal::_check_nonzero_size(rA))
        throw Exception("qpp::kron()", Exception::Type::ZERO_SIZE);

    // check zero-size
    if (!internal::_check_nonzero_size(rB))
        throw Exception("qpp::kron()", Exception::Type::ZERO_SIZE);

    std::size_t Acols = static_cast<std::size_t>(rA.cols());
    std::size_t Arows = static_cast<std::size_t>(rA.rows());
    std::size_t Bcols = static_cast<std::size_t>(rB.cols());
    std::size_t Brows = static_cast<std::size_t>(rB.rows());

    DynMat<typename Derived1::Scalar> result;
    result.resize(Arows * Brows, Acols * Bcols);

#pragma omp parallel for collapse(2)
    for (std::size_t j = 0; j < Acols; ++j)
        for (std::size_t i = 0; i < Arows; ++i)
            result.block(i * Brows, j * Bcols, Brows, Bcols) = rA(i, j) * rB;

    return result;

}

// may be useful, extracts variadic template argument pack into a std::vector
template<typename T>
// ends the recursion
void variadic_vector_emplace(std::vector<T>&)
{
}

// may be useful, extracts variadic template argument pack into a std::vector
template<typename T, typename First, typename ... Args>
void variadic_vector_emplace(std::vector<T>& v, First&& first, Args&& ... args)
{
    v.emplace_back(std::forward<First>(first));
    variadic_vector_emplace(v, std::forward<Args>(args)...);
}

} /* namespace internal */
} /* namespace qpp */

#endif /* INCLUDE_INTERNAL_FUNCTIONS_H_ */
