/*
 * Quantum++
 *
 * Copyright (c) 2013 - 2016 Vlad Gheorghiu (vgheorgh@gmail.com)
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

/**
* \file internal/util.h
* \brief Internal utility functions
*/

#ifndef INTERNAL_UTIL_H_
#define INTERNAL_UTIL_H_

namespace qpp
{
/**
* \namespace qpp::internal
* \brief Internal utility functions, do not use/modify
*/
namespace internal
{

// integer index to multi-index, use C-style array for speed
// standard lexicographical order, e.g. 00, 01, 10, 11
inline void _n2multiidx(idx n, idx numdims, const idx* dims, idx* result)
noexcept
{
    // no error checks to improve speed
    for (idx i = 0; i < numdims; ++i)
    {
        result[numdims - i - 1] = n % (dims[numdims - i - 1]);
        n /= (dims[numdims - i - 1]);
    }
}

// multi-index to integer index, use C-style array for speed,
// standard lexicographical order, e.g. 00->0, 01->1, 10->2, 11->3
inline idx _multiidx2n(const idx* midx, idx numdims, const idx* dims)
noexcept
{
    // no error checks to improve speed

    // Static allocation for speed!
    // double the size for matrices reshaped as vectors
    idx part_prod[2 * maxn];

    idx result = 0;
    part_prod[numdims - 1] = 1;
    for (idx i = 1; i < numdims; ++i)
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
    return A.rows() == A.cols();
}

// check whether input is a vector or not
template<typename Derived>
bool _check_vector(const Eigen::MatrixBase<Derived>& A)
{
    return A.rows() == 1 || A.cols() == 1;
}

// check whether input is a row vector or not
template<typename Derived>
bool _check_rvector(const Eigen::MatrixBase<Derived>& A)
{
    return A.rows() == 1;
}

// check whether input is a column vector or not
template<typename Derived>
bool _check_cvector(const Eigen::MatrixBase<Derived>& A)
{
    return A.cols() == 1;
}

// check non-zero size of object that supports size() function
template<typename T>
bool _check_nonzero_size(const T& x) noexcept
{
    return x.size() != 0;
}

// check that all sizes match
template<typename T1, typename T2>
bool _check_matching_sizes(const T1& lhs, const T2& rhs) noexcept
{
    return lhs.size() == rhs.size();
}

// check that dims is a valid dimension vector
inline bool _check_dims(const std::vector<idx>& dims)
{
    if (dims.size() == 0)
        return false;

    return std::find_if(std::begin(dims), std::end(dims),
                        [dims](idx i) -> bool
                        {
                            if (i == 0) return true;
                            else return false;
                        }) == std::end(dims);
}

// check that valid dims match the dimensions
// of valid (non-zero sized) square matrix
template<typename Derived>
bool _check_dims_match_mat(const std::vector<idx>& dims,
                           const Eigen::MatrixBase<Derived>& A)
{
    idx proddim = std::accumulate(std::begin(dims), std::end(dims),
                                  static_cast<idx>(1), std::multiplies<idx>());

    return proddim == static_cast<idx>(A.rows());
}

// check that valid dims match the dimensions of valid column vector
template<typename Derived>
bool _check_dims_match_cvect(const std::vector<idx>& dims,
                             const Eigen::MatrixBase<Derived>& V)
{
    idx proddim = std::accumulate(std::begin(dims), std::end(dims),
                                  static_cast<idx>(1), std::multiplies<idx>());

    return proddim == static_cast<idx>(V.rows());
}

// check that valid dims match the dimensions of valid row vector
template<typename Derived>
bool _check_dims_match_rvect(const std::vector<idx>& dims,
                             const Eigen::MatrixBase<Derived>& V)
{
    idx proddim = std::accumulate(std::begin(dims), std::end(dims),
                                  static_cast<idx>(1), std::multiplies<idx>());;

    return proddim == static_cast<idx>(V.cols());
}

// check that all elements in valid dims equal to dim
inline bool _check_eq_dims(const std::vector<idx>& dims, idx dim) noexcept
{
    for (idx i : dims)
        if (i != dim)
            return false;

    return true;
}

// check that subsys is valid with respect to valid dims
inline bool _check_subsys_match_dims(const std::vector<idx>& subsys,
                                     const std::vector<idx>& dims)
{
    //	// check non-zero sized subsystems
    //	if (subsys.size() == 0)
    //		return false;

    // check valid number of subsystems
    if (subsys.size() > dims.size())
        return false;

    // sort the subsystems
    std::vector<idx> subsyssort = subsys;
    std::sort(std::begin(subsyssort), std::end(subsyssort));

    // check duplicates
    if (std::unique(std::begin(subsyssort), std::end(subsyssort))
        != std::end(subsyssort))
        return false;

    // check range of subsystems
    return std::find_if(std::begin(subsyssort), std::end(subsyssort),
                        [dims](idx i) -> bool
                        {
                            return i > dims.size() - 1;
                        }) == std::end(subsyssort);
}

// check matrix is 2 x 2
template<typename Derived>
bool _check_qubit_matrix(const Eigen::MatrixBase<Derived>& A) noexcept
{
    return A.rows() == 2 && A.cols() == 2;
}

// check column vector is 2 x 1
template<typename Derived>
bool _check_qubit_cvector(const Eigen::MatrixBase<Derived>& V) noexcept
{
    return V.rows() == 2 && V.cols() == 1;
}

// check row vector is 1 x 2
template<typename Derived>
bool _check_qubit_rvector(const Eigen::MatrixBase<Derived>& V) noexcept
{
    return V.rows() == 1 && V.cols() == 2;
}

// check row vector is 1 x 2 or 2 x 1
template<typename Derived>
bool _check_qubit_vector(const Eigen::MatrixBase<Derived>& V) noexcept
{
    return (V.rows() == 1 && V.cols() == 2) ||
           (V.rows() == 2 && V.cols() == 1);
}


// check valid permutation
inline bool _check_perm(const std::vector<idx>& perm)
{
    if (perm.size() == 0)
        return false;

    std::vector<idx> ordered(perm.size());
    std::iota(std::begin(ordered), std::end(ordered), 0);

    return std::is_permutation(std::begin(ordered), std::end(ordered),
                               std::begin(perm));
}

// Kronecker product of 2 matrices, preserve return type
// internal function for the variadic template function wrapper kron()
template<typename Derived1, typename Derived2>
dyn_mat<typename Derived1::Scalar> _kron2(const Eigen::MatrixBase<Derived1>& A,
                                          const Eigen::MatrixBase<Derived2>& B)
{
    const dyn_mat<typename Derived1::Scalar>& rA = A;
    const dyn_mat<typename Derived2::Scalar>& rB = B;

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
    // END EXCEPTION CHECKS

    idx Acols = static_cast<idx>(rA.cols());
    idx Arows = static_cast<idx>(rA.rows());
    idx Bcols = static_cast<idx>(rB.cols());
    idx Brows = static_cast<idx>(rB.rows());

    dyn_mat<typename Derived1::Scalar> result;
    result.resize(Arows * Brows, Acols * Bcols);

#pragma omp parallel for collapse(2)
    for (idx j = 0; j < Acols; ++j) // column major order for speed
        for (idx i = 0; i < Arows; ++i)
            result.block(i * Brows, j * Bcols, Brows, Bcols) = rA(i, j) * rB;

    return result;
}

// Direct sum of 2 matrices, preserve return type
// internal function for the variadic template function wrapper dirsum()
template<typename Derived1, typename Derived2>
dyn_mat<typename Derived1::Scalar> _dirsum2(
        const Eigen::MatrixBase<Derived1>& A,
        const Eigen::MatrixBase<Derived2>& B)
{
    const dyn_mat<typename Derived1::Scalar>& rA = A;
    const dyn_mat<typename Derived2::Scalar>& rB = B;

    // EXCEPTION CHECKS

    // check types
    if (!std::is_same<typename Derived1::Scalar,
            typename Derived2::Scalar>::value)
        throw Exception("qpp::dirsum()", Exception::Type::TYPE_MISMATCH);

    // check zero-size
    if (!internal::_check_nonzero_size(rA))
        throw Exception("qpp::dirsum()", Exception::Type::ZERO_SIZE);

    // check zero-size
    if (!internal::_check_nonzero_size(rB))
        throw Exception("qpp::dirsum()", Exception::Type::ZERO_SIZE);
    // END EXCEPTION CHECKS

    idx Acols = static_cast<idx>(rA.cols());
    idx Arows = static_cast<idx>(rA.rows());
    idx Bcols = static_cast<idx>(rB.cols());
    idx Brows = static_cast<idx>(rB.rows());

    dyn_mat<typename Derived1::Scalar> result =
            dyn_mat<typename Derived1::Scalar>::Zero(Arows + Brows,
                                                     Acols + Bcols);

    result.block(0, 0, Arows, Acols) = rA;
    result.block(Arows, Acols, Brows, Bcols) = rB;

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

#endif /* INTERNAL_UTIL_H_ */
