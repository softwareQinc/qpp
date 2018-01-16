/*
 * This file is part of Quantum++.
 *
 * MIT License
 *
 * Copyright (c) 2013 - 2018 Vlad Gheorghiu (vgheorgh@gmail.com)
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
* \file internal/util.h
* \brief Internal utility functions
*/

#ifndef INTERNAL_UTIL_H_
#define INTERNAL_UTIL_H_

namespace qpp {
/**
* \namespace qpp::internal
* \brief Internal utility functions, do not use them directly or modify them
*/
namespace internal {
// integer index to multi-index, use C-style array for speed
// standard lexicographical order, e.g. 00, 01, 10, 11
inline void n2multiidx(idx n, idx numdims, const idx* const dims,
                       idx* result) noexcept {
// error checks only in DEBUG version
#ifndef NDEBUG
    if (numdims > 0) // numdims equal zero is a no-op
    {
        idx D = 1;
        for (idx i = 0; i < numdims; ++i)
            D *= dims[i];
        assert(n < D);
    }
#endif
    // no error checks in release version to improve speed
    for (idx i = 0; i < numdims; ++i) {
        result[numdims - i - 1] = n % (dims[numdims - i - 1]);
        n /= (dims[numdims - i - 1]);
    }
}

// silence g++4.9 bogus warning -Warray-bounds and -Wmaybe-uninitialized
// in qpp::internal::multiidx2n()
#if (__GNUC__ && !__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Warray-bounds"
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif
// multi-index to integer index, use C-style array for speed,
// standard lexicographical order, e.g. 00->0, 01->1, 10->2, 11->3
inline idx multiidx2n(const idx* const midx, idx numdims,
                      const idx* const dims) noexcept {
// error checks only in DEBUG version
#ifndef NDEBUG
    assert(numdims > 0);
#endif

    // no error checks in release version to improve speed

    // Static allocation for speed!
    // double the size for matrices reshaped as vectors
    idx part_prod[2 * maxn];

    idx result = 0;
    part_prod[numdims - 1] = 1;
    for (idx i = 1; i < numdims; ++i) {
        part_prod[numdims - i - 1] = part_prod[numdims - i] * dims[numdims - i];
        result += midx[numdims - i - 1] * part_prod[numdims - i - 1];
    }

    return result + midx[numdims - 1];
}
#if (__GNUC__ && !__clang__)
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

// check non-zero size of object that supports size() function
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
    if (dims.size() == 0)
        return false;

    return std::find_if(std::begin(dims), std::end(dims),
                        [dims](idx i) -> bool {
                            if (i == 0)
                                return true;
                            else
                                return false;
                        }) == std::end(dims);
}

// check that valid dims match the dimensions
// of valid (non-zero sized) square matrix
template <typename Derived>
bool check_dims_match_mat(const std::vector<idx>& dims,
                          const Eigen::MatrixBase<Derived>& A) {
// error checks only in DEBUG version
#ifndef NDEBUG
    assert(dims.size() > 0);
    assert(A.rows() == A.cols());
#endif
    idx proddim = std::accumulate(std::begin(dims), std::end(dims),
                                  static_cast<idx>(1), std::multiplies<idx>());

    return proddim == static_cast<idx>(A.rows());
}

// check that valid dims match the dimensions of valid column vector
template <typename Derived>
bool check_dims_match_cvect(const std::vector<idx>& dims,
                            const Eigen::MatrixBase<Derived>& A) {
// error checks only in DEBUG version
#ifndef NDEBUG
    assert(dims.size() > 0);
    assert(A.rows() > 0);
    assert(A.cols() == 1);
#endif
    idx proddim = std::accumulate(std::begin(dims), std::end(dims),
                                  static_cast<idx>(1), std::multiplies<idx>());

    return proddim == static_cast<idx>(A.rows());
}

// check that valid dims match the dimensions of valid row vector
template <typename Derived>
bool check_dims_match_rvect(const std::vector<idx>& dims,
                            const Eigen::MatrixBase<Derived>& A) {
// error checks only in DEBUG version
#ifndef NDEBUG
    assert(dims.size() > 0);
    assert(A.cols() > 0);
    assert(A.rows() == 1);
#endif
    idx proddim = std::accumulate(std::begin(dims), std::end(dims),
                                  static_cast<idx>(1), std::multiplies<idx>());
    ;

    return proddim == static_cast<idx>(A.cols());
}

// check that all elements in valid dims equal to dim
inline bool check_eq_dims(const std::vector<idx>& dims, idx dim) noexcept {
// error checks only in DEBUG version
#ifndef NDEBUG
    assert(dims.size() > 0);
#endif
    for (idx i : dims)
        if (i != dim)
            return false;

    return true;
}

// check that subsys is valid with respect to valid dims
inline bool check_subsys_match_dims(const std::vector<idx>& subsys,
                                    const std::vector<idx>& dims) {
    // subsys can be empty

    // check valid number of subsystems
    if (subsys.size() > dims.size())
        return false;

    // sort the subsystems
    std::vector<idx> subsyssort = subsys;
    std::sort(std::begin(subsyssort), std::end(subsyssort));

    // check duplicates
    if (std::unique(std::begin(subsyssort), std::end(subsyssort)) !=
        std::end(subsyssort))
        return false;

    // check range of subsystems
    return std::find_if(std::begin(subsyssort), std::end(subsyssort),
                        [dims](idx i) -> bool {
                            return i > dims.size() - 1;
                        }) == std::end(subsyssort);
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
    if (perm.size() == 0)
        return false;

    std::vector<idx> ordered(perm.size());
    std::iota(std::begin(ordered), std::end(ordered), 0);

    return std::is_permutation(std::begin(ordered), std::end(ordered),
                               std::begin(perm));
}

// Kronecker product of 2 matrices, preserve return type
// internal function for the variadic template function wrapper kron()
template <typename Derived1, typename Derived2>
dyn_mat<typename Derived1::Scalar> kron2(const Eigen::MatrixBase<Derived1>& A,
                                         const Eigen::MatrixBase<Derived2>& B) {
    const dyn_mat<typename Derived1::Scalar>& rA = A.derived();
    const dyn_mat<typename Derived2::Scalar>& rB = B.derived();

    // EXCEPTION CHECKS

    // check types
    if (!std::is_same<typename Derived1::Scalar,
                      typename Derived2::Scalar>::value)
        throw exception::TypeMismatch("qpp::kron()");

    // check zero-size
    if (!internal::check_nonzero_size(rA))
        throw exception::ZeroSize("qpp::kron()");

    // check zero-size
    if (!internal::check_nonzero_size(rB))
        throw exception::ZeroSize("qpp::kron()");
    // END EXCEPTION CHECKS

    idx Acols = static_cast<idx>(rA.cols());
    idx Arows = static_cast<idx>(rA.rows());
    idx Bcols = static_cast<idx>(rB.cols());
    idx Brows = static_cast<idx>(rB.rows());

    dyn_mat<typename Derived1::Scalar> result;
    result.resize(Arows * Brows, Acols * Bcols);

#ifdef WITH_OPENMP_
#pragma omp parallel for collapse(2)
#endif // WITH_OPENMP_
    // column major order for speed
    for (idx j = 0; j < Acols; ++j)
        for (idx i = 0; i < Arows; ++i)
            result.block(i * Brows, j * Bcols, Brows, Bcols) = rA(i, j) * rB;

    return result;
}

// Direct sum of 2 matrices, preserve return type
// internal function for the variadic template function wrapper dirsum()
template <typename Derived1, typename Derived2>
dyn_mat<typename Derived1::Scalar>
dirsum2(const Eigen::MatrixBase<Derived1>& A,
        const Eigen::MatrixBase<Derived2>& B) {
    const dyn_mat<typename Derived1::Scalar>& rA = A.derived();
    const dyn_mat<typename Derived2::Scalar>& rB = B.derived();

    // EXCEPTION CHECKS

    // check types
    if (!std::is_same<typename Derived1::Scalar,
                      typename Derived2::Scalar>::value)
        throw exception::TypeMismatch("qpp::dirsum()");

    // check zero-size
    if (!internal::check_nonzero_size(rA))
        throw exception::ZeroSize("qpp::dirsum()");

    // check zero-size
    if (!internal::check_nonzero_size(rB))
        throw exception::ZeroSize("qpp::dirsum()");
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
// dimension d) from an object (ket/bra/density matrix) of size sz
inline idx get_num_subsys(idx sz, idx d) {
// error checks only in DEBUG version
#ifndef NDEBUG
    assert(sz > 0);
    assert(d > 1);
#endif
    return static_cast<idx>(std::llround(std::log2(sz) / std::log2(d)));
}

// returns the dimension of a subsystem (each subsystem assumed of the same
// dimension d) from an object (ket/bra/density matrix) of size sz consisting
// of N subsystems
inline idx get_dim_subsys(idx sz, idx N) {
// error checks only in DEBUG version
#ifndef NDEBUG
    assert(N > 0);
    assert(sz > 0);
#endif
    if (N == 2)
        return static_cast<idx>(std::llround(std::sqrt(sz)));

    return static_cast<idx>(std::llround(std::pow(sz, 1. / N)));
}

// implementation details for pretty formatting
struct Display_Impl_ {
    template <typename T>
    // T must support rows(), cols(), operator()(idx, idx) const
    std::ostream& display_impl_(const T& A, std::ostream& os,
                                double chop = qpp::chop) const {
        std::ostringstream ostr;
        ostr.copyfmt(os); // copy os' state

        std::vector<std::string> vstr;
        std::string strA;

        for (idx i = 0; i < static_cast<idx>(A.rows()); ++i) {
            for (idx j = 0; j < static_cast<idx>(A.cols()); ++j) {
                strA.clear(); // clear the temporary string
                ostr.clear();
                ostr.str(std::string{}); // clear the ostringstream

                // convert to complex
                double re = static_cast<cplx>(A(i, j)).real();
                double im = static_cast<cplx>(A(i, j)).imag();

                if (std::abs(re) < chop && std::abs(im) < chop) {
                    ostr << "0 "; // otherwise segfault on destruction
                    // if using only vstr.push_back("0 ");
                    // bug in MATLAB libmx
                    vstr.push_back(ostr.str());
                } else if (std::abs(re) < chop) {
                    ostr << im;
                    vstr.push_back(ostr.str() + "i");
                } else if (std::abs(im) < chop) {
                    ostr << re;
                    vstr.push_back(ostr.str() + " ");
                } else {
                    ostr << re;
                    strA = ostr.str();

                    strA += (im > 0 ? " + " : " - ");
                    ostr.clear();
                    ostr.str(std::string()); // clear
                    ostr << std::abs(im);
                    strA += ostr.str();
                    strA += "i";
                    vstr.push_back(strA);
                }
            }
        }

        // determine the maximum lenght of the entries in each column
        std::vector<idx> maxlengthcols(A.cols(), 0);

        for (idx i = 0; i < static_cast<idx>(A.rows()); ++i)
            for (idx j = 0; j < static_cast<idx>(A.cols()); ++j)
                if (vstr[i * A.cols() + j].size() > maxlengthcols[j])
                    maxlengthcols[j] = vstr[i * A.cols() + j].size();

        // finally display it!
        for (idx i = 0; i < static_cast<idx>(A.rows()); ++i) {
            os << std::setw(static_cast<int>(maxlengthcols[0])) << std::right
               << vstr[i * A.cols()]; // display first column
            // then the rest
            for (idx j = 1; j < static_cast<idx>(A.cols()); ++j)
                os << std::setw(static_cast<int>(maxlengthcols[j] + 2))
                   << std::right << vstr[i * A.cols() + j];

            if (i < static_cast<idx>(A.rows()) - 1)
                os << std::endl;
        }

        return os;
    }
};

} /* namespace internal */
} /* namespace qpp */

#endif /* INTERNAL_UTIL_H_ */
