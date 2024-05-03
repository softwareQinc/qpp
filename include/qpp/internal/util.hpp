/*
 * This file is part of Quantum++.
 *
 * Copyright (c) 2017 - 2024 softwareQ Inc. All rights reserved.
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
 * \file internal/util.hpp
 * \brief Internal utility functions
 */

#ifndef QPP_INTERNAL_UTIL_HPP_
#define QPP_INTERNAL_UTIL_HPP_

#include <cassert>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <numeric>
#include <type_traits>

#include <Eigen/Dense>

#include "qpp/options.hpp"
#include "qpp/types.hpp"

#include "qpp/classes/exception.hpp"

namespace qpp {
/**
 * \namespace qpp::internal
 * \brief Internal utility functions, do not use them directly or modify them
 */
namespace internal {
// integer index to multi-index, use C-style array for speed
// standard lexicographical order, e.g., 00, 01, 10, 11
template <typename T, typename U = T, typename V = T>
[[qpp::critical]] void n2multiidx(T n, std::size_t numdims, const U* const dims,
                                  V* result) noexcept {
    static_assert(std::is_integral_v<T>, "T must be an integral value");
    static_assert(std::is_integral_v<U>, "U must be an integral value");
    static_assert(std::is_integral_v<V>, "V must be an integral value");

    // error checks only in DEBUG version
#ifndef NDEBUG
    if (numdims > 0) // numdims equal zero is a no-op
    {
        idx D = 1;
        for (std::size_t i = 0; i < numdims; ++i) {
            D *= dims[i];
        }
        assert(static_cast<idx>(n) < D);
    }
#endif
    // no error checks in release version to improve speed
    for (std::size_t i = 0; i < numdims; ++i) {
        result[numdims - i - 1] = n % (dims[numdims - i - 1]);
        n /= (dims[numdims - i - 1]);
    }
}

// silence g++4.9 bogus warning -Warray-bounds and -Wmaybe-uninitialized
// in qpp::internal::multiidx2n()
#if defined(__GNUC__) && !defined(__clang__) && !defined(__INTEL_COMPILER)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Warray-bounds"
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif
// multi-index to integer index, use C-style array for speed,
// standard lexicographical order, e.g., 00->0, 01->1, 10->2, 11->3
template <typename V, typename T = V, typename U = T>
[[qpp::critical]] T multiidx2n(const V* const midx, std::size_t numdims,
                               const U* const dims) noexcept {
    static_assert(std::is_integral_v<T>, "T must be an integral value");
    static_assert(std::is_integral_v<U>, "U must be an integral value");
    static_assert(std::is_integral_v<V>, "V must be an integral value");

    // error checks only in DEBUG version
    assert(numdims > 0);
    assert(numdims < internal::maxn);
#ifndef NDEBUG
    for (std::size_t i = 0; i < numdims; ++i) {
        assert(static_cast<idx>(midx[i]) < dims[i]);
    }
#endif
    // no error checks in release version to improve speed

    T part_prod = 1;
    T result = 0;
    for (std::size_t i = 1; i < numdims; ++i) {
        part_prod = part_prod * dims[numdims - i];
        result += midx[numdims - i - 1] * part_prod;
    }

    return result + midx[numdims - 1];
}
#if defined(__GNUC__) && !defined(__clang__) && !defined(__INTEL_COMPILER)
#pragma GCC diagnostic pop
#endif

// check square matrix
template <typename Derived>
bool check_square_mat(const Eigen::MatrixBase<Derived>& A) {
    return A.rows() == A.cols();
}

// check whether input is a vector or not
template <typename Derived>
bool check_vector(const Eigen::MatrixBase<Derived>& A) {
    return A.rows() == 1 || A.cols() == 1;
}

// check whether input is a row vector or not
template <typename Derived>
bool check_rvector(const Eigen::MatrixBase<Derived>& A) {
    return A.rows() == 1;
}

// check whether input is a column vector or not
template <typename Derived>
bool check_cvector(const Eigen::MatrixBase<Derived>& A) {
    return A.cols() == 1;
}

// check non-zero size of object that implements size() member function
template <typename T>
bool check_nonzero_size(const T& x) noexcept {
    return x.size() != 0;
}

// check that all sizes match
template <typename T1, typename T2>
bool check_matching_sizes(const T1& lhs, const T2& rhs) noexcept {
    return lhs.size() == rhs.size();
}

// check that dims is a valid dimension vector
inline bool check_dims(const std::vector<idx>& dims) {
    if (dims.empty()) {
        return false;
    }

    return std::find_if(dims.begin(), dims.end(), [dims](idx i) -> bool {
               return (i == 0);
           }) == dims.end();
}

// check that valid dims match the dimensions
// of valid (non-zero sized) square matrix
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

// check that valid dims match the dimensions of valid column vector
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

// check that valid dims match the dimensions of valid row vector
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

// check that all elements in valid dims equal to dim
inline bool check_eq_dims(const std::vector<idx>& dims, idx dim) noexcept {
    // error checks only in DEBUG version
    assert(!dims.empty());

    return std::all_of(dims.begin(), dims.end(),
                       [dim](idx i) { return i == dim; });
}

// check that vector has no duplicates
inline bool check_no_duplicates(std::vector<idx> v) {
    std::sort(v.begin(), v.end());
    if (std::unique(v.begin(), v.end()) != v.end()) {
        return false;
    } else {
        return true;
    }
}

// check that subsys is valid with respect to valid dims
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

// check matrix is 2 x 2
template <typename Derived>
bool check_qubit_matrix(const Eigen::MatrixBase<Derived>& A) noexcept {
    return A.rows() == 2 && A.cols() == 2;
}

// check column vector is 2 x 1
template <typename Derived>
bool check_qubit_cvector(const Eigen::MatrixBase<Derived>& A) noexcept {
    return A.rows() == 2 && A.cols() == 1;
}

// check row vector is 1 x 2
template <typename Derived>
bool check_qubit_rvector(const Eigen::MatrixBase<Derived>& A) noexcept {
    return A.rows() == 1 && A.cols() == 2;
}

// check row vector is 1 x 2 or 2 x 1
template <typename Derived>
bool check_qubit_vector(const Eigen::MatrixBase<Derived>& A) noexcept {
    return (A.rows() == 1 && A.cols() == 2) || (A.rows() == 2 && A.cols() == 1);
}

// check valid permutation
inline bool check_perm(const std::vector<idx>& perm) {
    if (perm.empty()) {
        return false;
    }

    std::vector<idx> ordered(perm.size());
    std::iota(ordered.begin(), ordered.end(), 0);

    return std::is_permutation(ordered.begin(), ordered.end(), perm.begin());
}

// Kronecker product of 2 matrices, preserve return type
// internal function for the variadic template function wrapper qpp::kron()
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
// NOLINTNEXTLINE
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

// Direct sum of 2 matrices, preserve return type
// internal function for the variadic template function wrapper qpp::dirsum()
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

// may be useful, extracts variadic template argument pack into a std::vector
template <typename T>
// ends the recursion
void variadic_vector_emplace(std::vector<T>&) {}

// may be useful, extracts variadic template argument pack into a std::vector
template <typename T, typename First, typename... Args>
void variadic_vector_emplace(std::vector<T>& v, First&& first, Args&&... args) {
    v.emplace_back(std::forward<First>(first));
    variadic_vector_emplace(v, std::forward<Args>(args)...);
}

// returns the number of subsystems (each subsystem assumed of the same
// dimension d) from an object (ket/bra/matrix) of size D
inline idx get_num_subsys(idx D, idx d) {
    // error checks only in DEBUG version
    assert(D > 0);
    assert(d > 1);

    auto n = static_cast<idx>(std::llround(std::log2(D) / std::log2(d)));

    return n;
}

// returns the dimension of a subsystem (each subsystem assumed of the same
// dimension d) from an object (ket/bra/matrix) of size D consisting of n
// subsystems
inline idx get_dim_subsys(idx D, idx n) {
    // error checks only in DEBUG version
    assert(n > 0);
    assert(D > 0);

    auto d = (n == 2) ? static_cast<idx>(std::llround(std::sqrt(D)))
                      : static_cast<idx>(std::llround(
                            std::pow(D, 1. / static_cast<realT>(n))));

    return d;
}

// safe power of idx types, a^b
template <typename T = idx>
inline T safe_pow(T a, T b) {
    return static_cast<T>(std::llround(std::pow(a, b)));
}

// chops a floating-point or complex number to zero, returns it unchanged
// otherwise
template <typename T>
T abs_float_or_cplx_chop(const T& x, realT chop) {
    if constexpr (std::numeric_limits<T>::is_iec559 || is_complex_v<T>) {
        if (std::abs(x) < chop) {
            return 0;
        }
    }
    return x;
}

// self-documented
template <typename T>
constexpr bool is_negative(T t) {
    return t < 0;
}

// converts a real value to a text representation (to machine precision)
template <typename T>
std::string real2text(T d) {
    std::stringstream ss;
    ss << std::setprecision(std::numeric_limits<T>::max_digits10);
    ss << d;

    return ss.str();
}

// converts the text representation of a real value to its corresponding
// numerical value
template <typename T>
T text2real(const std::string& str) {
    return std::strtod(str.c_str(), nullptr);
}

// returns true if the index i as a multi-index contains the multi-index
// dits
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

// computes (|dits><dits| \otimes I)|psi>
template <typename Derived>
dyn_col_vect<Derived> project_ket_on_dits(dyn_col_vect<Derived> psi,
                                          const std::vector<idx>& dits,
                                          const std::vector<idx>& subsys,
                                          const std::vector<idx>& dims, idx D) {
#ifdef QPP_OPENMP
// NOLINTNEXTLINE
#pragma omp parallel for
#endif // QPP_OPENMP
    for (idx i = 0; i < D; ++i) {
        if (!idx_contains_dits(i, dits, subsys, dims)) {
            psi(i) = 0;
        }
    }

    return psi;
}
} /* namespace internal */

} /* namespace qpp */

#endif /* QPP_INTERNAL_UTIL_HPP_ */
