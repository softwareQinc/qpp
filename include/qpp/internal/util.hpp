/*
 * This file is part of Quantum++.
 *
 * Copyright (c) 2017 - 2026 softwareQ Inc. All rights reserved.
 *
 * MIT License
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

/**
 * @file qpp/internal/util.hpp
 * @brief Internal utility functions
 */

#ifndef QPP_INTERNAL_UTIL_HPP_
#define QPP_INTERNAL_UTIL_HPP_

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <numeric>
#include <type_traits>

#include <Eigen/Dense>

#include "qpp/options.hpp"
#include "qpp/traits.hpp"
#include "qpp/types.hpp"

#include "qpp/classes/exception.hpp"

namespace qpp {
/**
 * @namespace qpp::internal
 * @brief Internal functions, do not use them directly or modify them
 */
namespace internal {
/**
 * @brief Converts an integer index to a multi-index
 *
 * Uses lexicographical order, with the last index varying fastest
 * For dimensions `{2, 2}`, the order is `00, 01, 10, 11`
 *
 * @tparam T Integral type of the input index
 * @tparam U Integral type of the dimensions
 * @tparam V Integral type of the output indices
 *
 * @param n Input integer index
 * @param numdims Number of dimensions
 * @param dims Dimension array of size \a numdims
 * @param result Output multi-index array of size \a numdims
 */
template <typename T, typename U = T, typename V = T>
[[qpp::critical]] inline void n2multiidx(T n, std::size_t numdims,
                                         const U* dims, V* result) noexcept {
    static_assert(std::is_integral_v<T>, "T must be an integral value");
    static_assert(std::is_integral_v<U>, "U must be an integral value");
    static_assert(std::is_integral_v<V>, "V must be an integral value");

    // error checks only in DEBUG version
#ifndef NDEBUG
    if (numdims != 0) {
        std::size_t D = 1;
        for (std::size_t i = 0; i < numdims; ++i) {
            D *= static_cast<std::size_t>(dims[i]);
        }
        assert(static_cast<std::size_t>(n) < D && "n is out of bounds");
    }
#endif

    // main loop: unrolled directionally for cache locality
    for (std::size_t i = numdims; i-- > 0;) {
        const auto dim = static_cast<T>(dims[i]);
        result[i] = static_cast<V>(n % dim);
        n /= dim;
    }
}

/**
 * @brief Converts a multi-index to an integer index
 *
 * Uses lexicographical order, with the last index varying fastest
 * For dimensions `{2, 2}`, the order is `00 -> 0, 01 -> 1, 10 -> 2, 11 -> 3`
 *
 * @tparam V Integral type of the input indices
 * @tparam T Integral type of the output index
 * @tparam U Integral type of the dimensions
 *
 * @param midx Input multi-index array of size \a numdims
 * @param numdims Number of dimensions
 * @param dims Dimension array of size \a numdims
 * @return Integer index
 */
template <typename V, typename T = V, typename U = T>
[[qpp::critical]] T multiidx2n(const V* const midx, std::size_t numdims,
                               const U* const dims) noexcept {
    static_assert(std::is_integral_v<T>, "T must be an integral value");
    static_assert(std::is_integral_v<U>, "U must be an integral value");
    static_assert(std::is_integral_v<V>, "V must be an integral value");

    // error checks only in DEBUG version
    assert(numdims > 0 && numdims < internal::maxn);
#ifndef NDEBUG
    for (std::size_t i = 0; i < numdims; ++i) {
        assert(static_cast<std::size_t>(midx[i]) <
               static_cast<std::size_t>(dims[i]));
    }
#endif

    T result = 0;
    T stride = 1;
    for (std::size_t i = numdims; i-- > 0;) {
        result += static_cast<T>(midx[i]) * stride;
        stride *= static_cast<T>(dims[i]);
    }
    return result;
}

/**
 * @brief Checks whether a matrix is square
 *
 * @tparam Derived Eigen expression type
 * @param A Matrix to check
 * @return True if \a A is square, false otherwise
 */
template <typename Derived>
bool check_square_mat(const Eigen::MatrixBase<Derived>& A) {
    return A.rows() == A.cols();
}

/**
 * @brief Checks whether a matrix is a row vector
 *
 * @tparam Derived Eigen expression type
 * @param A Matrix to check
 * @return True if \a A is a row vector, false otherwise
 */
template <typename Derived>
bool check_rvector(const Eigen::MatrixBase<Derived>& A) {
    return A.rows() == 1;
}

/**
 * @brief Checks whether a matrix is a column vector
 *
 * @tparam Derived Eigen expression type
 * @param A Matrix to check
 * @return True if \a A is a column vector, false otherwise
 */
template <typename Derived>
bool check_cvector(const Eigen::MatrixBase<Derived>& A) {
    return A.cols() == 1;
}

/**
 * @brief Checks whether a matrix is a vector
 *
 * @tparam Derived Eigen expression type
 * @param A Matrix to check
 * @return True if \a A is a row or column vector, false otherwise
 */
template <typename Derived>
bool check_vector(const Eigen::MatrixBase<Derived>& A) {
    return A.rows() == 1 || A.cols() == 1;
}

/**
 * @brief Checks whether an object has non-zero size
 *
 * @tparam T Type with a size() member function
 * @param x Object to check
 * @return True if \a x has non-zero size, false otherwise
 */
template <typename T>
bool check_nonzero_size(const T& x) noexcept {
    return x.size() != 0;
}

/**
 * @brief Checks whether two objects have matching sizes
 *
 * @tparam T1 Type with a size() member function
 * @tparam T2 Type with a size() member function
 * @param lhs First object to check
 * @param rhs Second object to check
 * @return True if \a lhs and \a rhs have the same size, false otherwise
 */
template <typename T1, typename T2>
bool check_matching_sizes(const T1& lhs, const T2& rhs) noexcept {
    return lhs.size() == rhs.size();
}

/**
 * @brief Checks whether a dimension vector is valid
 *
 * @param dims Dimension vector to check
 * @return True if \a dims is non-empty and has no zero entries, false otherwise
 */
inline bool check_dims(const std::vector<idx>& dims) {
    if (dims.empty()) {
        return false;
    }

    return std::find_if(dims.begin(), dims.end(), [dims](idx i) -> bool {
               return (i == 0);
           }) == dims.end();
}

/**
 * @brief Checks whether dimensions match a square matrix
 *
 * @tparam Derived Eigen expression type
 * @param dims Dimension vector to check
 * @param A Square matrix to check
 * @return True if the product of \a dims matches the size of \a A, false
 * otherwise
 */
template <typename Derived>
bool check_dims_match_mat(const std::vector<idx>& dims,
                          const Eigen::MatrixBase<Derived>& A) {
    // error checks only in DEBUG version
    assert(!dims.empty());
    assert(A.rows() == A.cols());

    idx proddim = std::accumulate(dims.begin(), dims.end(), static_cast<idx>(1),
                                  std::multiplies<>());

    return proddim == static_cast<idx>(A.cols());
}

/**
 * @brief Checks whether dimensions match a column vector
 *
 * @tparam Derived Eigen expression type
 * @param dims Dimension vector to check
 * @param A Column vector to check
 * @return True if the product of \a dims matches the size of \a A, false
 * otherwise
 */
template <typename Derived>
bool check_dims_match_cvect(const std::vector<idx>& dims,
                            const Eigen::MatrixBase<Derived>& A) {
    // error checks only in DEBUG version
    assert(!dims.empty());
    assert(A.rows() > 0);
    assert(A.cols() == 1);

    idx proddim = std::accumulate(dims.begin(), dims.end(), static_cast<idx>(1),
                                  std::multiplies<>());

    return proddim == static_cast<idx>(A.rows());
}

/**
 * @brief Checks whether dimensions match a row vector
 *
 * @tparam Derived Eigen expression type
 * @param dims Dimension vector to check
 * @param A Row vector to check
 * @return True if the product of \a dims matches the size of \a A, false
 * otherwise
 */
template <typename Derived>
bool check_dims_match_rvect(const std::vector<idx>& dims,
                            const Eigen::MatrixBase<Derived>& A) {
    // error checks only in DEBUG version
    assert(!dims.empty());
    assert(A.cols() > 0);
    assert(A.rows() == 1);

    idx proddim = std::accumulate(dims.begin(), dims.end(), static_cast<idx>(1),
                                  std::multiplies<>());

    return proddim == static_cast<idx>(A.cols());
}

/**
 * @brief Checks whether dimensions match a vector
 *
 * @tparam Derived Eigen expression type
 * @param dims Dimension vector to check
 * @param A Row or column vector to check
 * @return True if the product of \a dims matches the size of \a A, false
 * otherwise
 */
template <typename Derived>
bool check_dims_match_vect(const std::vector<idx>& dims,
                           const Eigen::MatrixBase<Derived>& A) {
    // error checks only in DEBUG version
    assert(!dims.empty());
    assert(A.rows() > 0);
    assert(A.cols() > 0);
    assert(A.rows() == 1 || A.cols() == 1);

    const idx proddim = std::accumulate(
        dims.begin(), dims.end(), static_cast<idx>(1), std::multiplies<>());

    return proddim == static_cast<idx>(A.size());
}

/**
 * @brief Checks whether all dimensions are equal to a given value
 *
 * @param dims Dimension vector to check
 * @param dim Dimension value to compare against
 * @return True if all entries of \a dims are equal to \a dim, false otherwise
 */
inline bool check_eq_dims(const std::vector<idx>& dims, idx dim) noexcept {
    // error checks only in DEBUG version
    assert(!dims.empty());

    return std::all_of(dims.begin(), dims.end(),
                       [dim](idx i) { return i == dim; });
}

/**
 * @brief Checks whether a vector has no duplicate entries
 *
 * @param v Vector to check
 * @return True if \a v has no duplicate entries, false otherwise
 */
inline bool check_no_duplicates(std::vector<idx> v) {
    std::sort(v.begin(), v.end());
    if (std::unique(v.begin(), v.end()) != v.end()) {
        return false;
    } else {
        return true;
    }
}

/**
 * @brief Checks whether subsystems match dimensions
 *
 * @param subsys Subsystem vector to check
 * @param dims Dimension vector to check against
 * @return True if \a subsys has no duplicates and is in range, false otherwise
 */
inline bool check_subsys_match_dims(const std::vector<idx>& subsys,
                                    const std::vector<idx>& dims) {
    // subsys can be empty

    // check valid number of subsystems
    if (subsys.size() > dims.size()) {
        return false;
    }

    // check no duplicates
    if (!check_no_duplicates(subsys)) {
        return false;
    }

    // check range of subsystems
    return std::find_if(subsys.begin(), subsys.end(), [dims](idx i) -> bool {
               return i + 1 > static_cast<idx>(dims.size());
           }) == subsys.end();
}

/**
 * @brief Checks whether a matrix is 2 x 2
 *
 * @tparam Derived Eigen expression type
 * @param A Matrix to check
 * @return True if \a A is 2 x 2, false otherwise
 */
template <typename Derived>
bool check_qubit_matrix(const Eigen::MatrixBase<Derived>& A) noexcept {
    return A.rows() == 2 && A.cols() == 2;
}

/**
 * @brief Checks whether a column vector is 2 x 1
 *
 * @tparam Derived Eigen expression type
 * @param A Column vector to check
 * @return True if \a A is 2 x 1, false otherwise
 */
template <typename Derived>
bool check_qubit_cvector(const Eigen::MatrixBase<Derived>& A) noexcept {
    return A.rows() == 2 && A.cols() == 1;
}

/**
 * @brief Checks whether a row vector is 1 x 2
 *
 * @tparam Derived Eigen expression type
 * @param A Row vector to check
 * @return True if \a A is 1 x 2, false otherwise
 */
template <typename Derived>
bool check_qubit_rvector(const Eigen::MatrixBase<Derived>& A) noexcept {
    return A.rows() == 1 && A.cols() == 2;
}

/**
 * @brief Checks whether a vector is 1 x 2 or 2 x 1
 *
 * @tparam Derived Eigen expression type
 * @param A Vector to check
 * @return True if \a A is 1 x 2 or 2 x 1, false otherwise
 */
template <typename Derived>
bool check_qubit_vector(const Eigen::MatrixBase<Derived>& A) noexcept {
    return (A.rows() == 1 && A.cols() == 2) || (A.rows() == 2 && A.cols() == 1);
}

/**
 * @brief Checks whether a vector is a valid permutation
 *
 * @param perm Permutation vector to check
 * @return True if \a perm is a valid permutation, false otherwise
 */
inline bool check_perm(const std::vector<idx>& perm) {
    if (perm.empty()) {
        return false;
    }

    std::vector<idx> ordered(perm.size());
    std::iota(ordered.begin(), ordered.end(), 0);

    return std::is_permutation(ordered.begin(), ordered.end(), perm.begin());
}

/**
 * @brief Computes the Kronecker product of two matrices
 *
 * Internal function used by qpp::kron()
 *
 * @tparam Derived1 Eigen expression type of the first matrix
 * @tparam Derived2 Eigen expression type of the second matrix
 * @param A First matrix
 * @param B Second matrix
 * @return Kronecker product of \a A and \a B
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] dyn_mat<typename Derived1::Scalar>
kron2(const Eigen::MatrixBase<Derived1>& A,
      const Eigen::MatrixBase<Derived2>& B) {
    const dyn_mat<typename Derived1::Scalar>& rA = A.derived();
    const dyn_mat<typename Derived2::Scalar>& rB = B.derived();

    // EXCEPTION CHECKS
    // check types
    if (!std::is_same_v<typename Derived1::Scalar, typename Derived2::Scalar>) {
        throw exception::TypeMismatch("qpp::kron()", "A/B");
    }

    // check zero-size
    if (!internal::check_nonzero_size(rA)) {
        throw exception::ZeroSize("qpp::kron()", "A");
    }
    // check zero-size
    if (!internal::check_nonzero_size(rB)) {
        throw exception::ZeroSize("qpp::kron()", "B");
    }
    // END EXCEPTION CHECKS

    idx Acols = static_cast<idx>(rA.cols());
    idx Arows = static_cast<idx>(rA.rows());
    idx Bcols = static_cast<idx>(rB.cols());
    idx Brows = static_cast<idx>(rB.rows());

    dyn_mat<typename Derived1::Scalar> result;
    result.resize(Arows * Brows, Acols * Bcols);

#ifdef QPP_OPENMP
#pragma omp parallel for collapse(2)
#endif // QPP_OPENMP
    // column major order for speed
    for (idx j = 0; j < Acols; ++j) {
        for (idx i = 0; i < Arows; ++i) {
            result.block(i * Brows, j * Bcols, Brows, Bcols) = rA(i, j) * rB;
        }
    }

    return result;
}

/**
 * @brief Computes the direct sum of two matrices
 *
 * Internal function used by qpp::dirsum()
 *
 * @tparam Derived1 Eigen expression type of the first matrix
 * @tparam Derived2 Eigen expression type of the second matrix
 * @param A First matrix
 * @param B Second matrix
 * @return Direct sum of \a A and \a B
 */
template <typename Derived1, typename Derived2>
dyn_mat<typename Derived1::Scalar>
dirsum2(const Eigen::MatrixBase<Derived1>& A,
        const Eigen::MatrixBase<Derived2>& B) {
    const dyn_mat<typename Derived1::Scalar>& rA = A.derived();
    const dyn_mat<typename Derived2::Scalar>& rB = B.derived();

    // EXCEPTION CHECKS
    // check types
    if (!std::is_same_v<typename Derived1::Scalar, typename Derived2::Scalar>) {
        throw exception::TypeMismatch("qpp::dirsum()", "A/B");
    }

    // check zero-size
    if (!internal::check_nonzero_size(rA)) {
        throw exception::ZeroSize("qpp::dirsum()", "A");
    }
    // check zero-size
    if (!internal::check_nonzero_size(rB)) {
        throw exception::ZeroSize("qpp::dirsum()", "B");
    }
    // END EXCEPTION CHECKS

    idx Acols = static_cast<idx>(rA.cols());
    idx Arows = static_cast<idx>(rA.rows());
    idx Bcols = static_cast<idx>(rB.cols());
    idx Brows = static_cast<idx>(rB.rows());

    dyn_mat<typename Derived1::Scalar> result =
        dyn_mat<typename Derived1::Scalar>::Zero(Arows + Brows, Acols + Bcols);

    result.block(0, 0, Arows, Acols) = rA;
    result.block(Arows, Acols, Brows, Bcols) = rB;

    return result;
}

/**
 * @brief Ends variadic argument extraction into a vector
 *
 * @tparam T Vector value type
 * @param v Vector to fill
 */
template <typename T>
void variadic_vector_emplace(std::vector<T>&) {} // ends the recursion

/**
 * @brief Extracts variadic arguments into a vector
 *
 * @tparam T Vector value type
 * @tparam First First argument type
 * @tparam Args Remaining argument types
 * @param v Vector to fill
 * @param first First argument to add
 * @param args Remaining arguments to add
 */
template <typename T, typename First, typename... Args>
void variadic_vector_emplace(std::vector<T>& v, First&& first, Args&&... args) {
    v.emplace_back(std::forward<First>(first));
    variadic_vector_emplace(v, std::forward<Args>(args)...);
}

/**
 * @brief Returns the number of subsystems of equal dimension
 *
 * @param D Total object size
 * @param d Subsystem dimension
 * @return Number of subsystems
 */
inline idx get_num_subsys(idx D, idx d) {
    // error checks only in DEBUG version
    assert(D > 0);
    assert(d > 1);

    auto n = static_cast<idx>(std::llround(std::log2(D) / std::log2(d)));

    return n;
}

/**
 * @brief Returns the dimension of each subsystem
 *
 * Assumes all subsystems have equal dimension
 *
 * @param D Total object size
 * @param n Number of subsystems
 * @return Subsystem dimension
 */
inline idx get_dim_subsys(idx D, idx n) {
    // error checks only in DEBUG version
    assert(n > 0);
    assert(D > 0);

    auto d = (n == 2) ? static_cast<idx>(std::llround(std::sqrt(D)))
                      : static_cast<idx>(std::llround(
                            std::pow(D, 1. / static_cast<realT>(n))));

    return d;
}

/**
 * @brief Computes an integer power and rounds the result
 *
 * Uses std::pow() and rounds the floating-point result back to type \a T
 *
 * @tparam T Integral type
 * @param a Base
 * @param b Exponent
 * @return Rounded value of \a a raised to the power \a b
 */
template <typename T = idx>
inline T ipow_rounded(T a, T b) {
    return static_cast<T>(std::llround(std::pow(a, b)));
}

/**
 * @brief Chops a floating-point or complex value to zero
 *
 * Returns other types unchanged
 *
 * @tparam T Input type
 * @param x Value to check
 * @param chop Chopping threshold
 * @return Zero if \a x is below \a chop, otherwise \a x
 */
template <typename T>
T abs_float_or_cplx_chop(const T& x, realT chop) {
    if constexpr (std::numeric_limits<T>::is_iec559 || is_complex_v<T>) {
        if (std::abs(x) < chop) {
            return 0;
        }
    }

    return x;
}

/**
 * @brief Checks whether a value is negative
 *
 * @tparam T Input type
 * @param t Value to check
 * @return True if \a t is negative, false otherwise
 */
template <typename T>
constexpr bool is_negative(T t) {
    return t < 0;
}

/**
 * @brief Converts a real value to text
 *
 * Uses enough digits to preserve machine precision
 *
 * @tparam T Real type
 * @param d Value to convert
 * @return Text representation of \a d
 */
template <typename T>
std::string real2text(T d) {
    std::stringstream ss;
    ss << std::setprecision(std::numeric_limits<T>::max_digits10);
    ss << d;

    return ss.str();
}

/**
 * @brief Converts text to a real value
 *
 * @tparam T Real type
 * @param str Text representation to convert
 * @return Real value represented by \a str
 */
template <typename T>
T text2real(const std::string& str) {
    return std::strtod(str.c_str(), nullptr);
}

/**
 * @brief Checks whether an index contains given dits on selected subsystems
 *
 * @param i Integer index to check
 * @param dits Dits to compare against
 * @param subsys Subsystems where the dits are checked
 * @param dims Dimension vector
 * @return True if \a i contains \a dits on \a subsys, false otherwise
 */
inline bool idx_contains_dits(idx i, const std::vector<idx>& dits,
                              const std::vector<idx>& subsys,
                              const std::vector<idx>& dims) {
    idx Cstorage[internal::maxn];
    idx n = dims.size();
    idx subsys_size = subsys.size();

    /* get the col multi-indexes of the complement */
    internal::n2multiidx(i, n, dims.data(), Cstorage);
    std::vector<idx> midx_i(Cstorage, Cstorage + n);

    for (idx m = 0; m < subsys_size; m++) {
        if (midx_i[subsys[m]] != dits[m]) {
            return false;
        }
    }

    return true;
}

/**
 * @brief Projects a ket onto given dits on selected subsystems
 *
 * Computes
 * \f$(|\mathrm{dits}\rangle\langle\mathrm{dits}| \otimes I)|\psi\rangle\f$
 *
 * @tparam Derived Scalar type
 * @param psi Ket to project
 * @param dits Dits to project onto
 * @param subsys Subsystems where the dits are fixed
 * @param dims Dimension vector
 * @param D Total dimension
 * @return Projected ket
 */
template <typename Derived>
dyn_col_vect<Derived> project_ket_on_dits(dyn_col_vect<Derived> psi,
                                          const std::vector<idx>& dits,
                                          const std::vector<idx>& subsys,
                                          const std::vector<idx>& dims, idx D) {
#ifdef QPP_OPENMP
#pragma omp parallel for
#endif // QPP_OPENMP
    for (idx i = 0; i < D; ++i) {
        if (!idx_contains_dits(i, dits, subsys, dims)) {
            psi(i) = 0;
        }
    }

    return psi;
}

/**
 * @brief Checks whether all dimensions are equal to a given value
 *
 * @tparam Container Iterable container type
 * @param c Dimension container to check
 * @param d Dimension value to compare against
 * @return True if all entries of \a c are equal to \a d, false otherwise
 */
template <typename Container,
          typename = std::enable_if_t<is_iterable_v<Container>>>
bool all_dims_equal(const Container& c, idx d) {
    return std::all_of(std::begin(c), std::end(c),
                       [d](const auto& elem) { return elem == d; });
}

/**
 * @brief Computes the complement of a subset of indices relative to the full
 * set [0, 1, ..., n-1].
 * This is an O(n) time complexity approach using a boolean marker array.
 * @tparam T Underlying vector type
 * @param target The subset of indices to exclude (subsystems to trace out)
 * @param n The total number of subsystems (size of the full set)
 * @return The complement of \a target
 */
template <typename T>
inline std::vector<T> fast_complement(const std::vector<T>& target, idx n) {
    idx target_size = target.size();
    assert(target_size <= n);

    // Initialize a marker array (equivalent to a bitset) for fast lookup.
    // std::vector<bool> is highly optimized and space-efficient.
    std::vector<bool> is_target(n, false);

    // Mark the indices that are present in the target set. O(M) time, where
    // M is target.size(). This loop is parallelized as there are no write
    // conflicts.
#ifdef QPP_OPENMP
#pragma omp parallel for
#endif // QPP_OPENMP
    for (idx i = 0; i < target_size; ++i) {
        idx t = target[i];
        // Assuming 't' is validated to be between 0 and n-1
        is_target[t] = true;
    }

    std::vector<T> target_bar;
    target_bar.reserve(n - target_size);

    for (idx i = 0; i < n; ++i) {
        if (!is_target[i]) {
            target_bar.push_back(i);
        }
    }

    return target_bar;
}

} /* namespace internal */
} /* namespace qpp */

#endif /* QPP_INTERNAL_UTIL_HPP_ */
