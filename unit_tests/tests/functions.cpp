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

#include <vector>
#include "gtest/gtest.h"
#include "qpp.h"

using namespace qpp;

// Unit testing "functions.h"

/******************************************************************************/
/// BEGIN template<typename Derived> cmat qpp::absm(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_absm, AllTests) {}
/******************************************************************************/
/// BEGIN template<typename Container> std::vector<double> qpp::abssq(
///       const Container& c,
///       typename std::enable_if<is_iterable<Container>::value>::type*
///       = nullptr)
TEST(qpp_abssq_container, AllTests) {}
/******************************************************************************/
/// BEGIN template<typename Derived> std::vector<double> qpp::abssq(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_abssq_matrix, AllTests) {}
/******************************************************************************/
/// BEGIN template<typename InputIterator> std::vector<double> qpp::abssq(
///       InputIterator first, InputIterator last)
TEST(qpp_abssq_iterator, AllTests) {}
/******************************************************************************/
/// BEGIN template<typename Derived> dyn_mat<typename Derived::Scalar>
///       qpp::adjoint(const Eigen::MatrixBase<Derived>& A)
TEST(qpp_adjoint, AllTests) {}
/******************************************************************************/
/// BEGIN template<typename Derived1, typename Derived2>
///       dyn_mat<typename Derived1::Scalar> qpp::anticomm(
///       const Eigen::MatrixBase<Derived1>& A,
///       const Eigen::MatrixBase<Derived2>& B)
TEST(qpp_anticomm, AllTests) {}
/******************************************************************************/
/// BEGIN inline cmat qpp::bloch2rho(const std::vector<double>& r)
TEST(qpp_bloch2rho, AllTests) {}
/******************************************************************************/
/// BEGIN template<typename Derived1, typename Derived2>
///       dyn_mat<typename Derived1::Scalar> qpp::comm(
///       const Eigen::MatrixBase<Derived1>& A,
///       const Eigen::MatrixBase<Derived2>& B)
TEST(qpp_comm, AllTests) {}
/******************************************************************************/
/// BEGIN template<typename T> std::vector<T> qpp::complement(
///       std::vector<T> subsys, idx N)
TEST(qpp_complement, AllTests) {}
/******************************************************************************/
/// BEGIN template<typename Derived> dyn_mat<typename Derived::Scalar>
///       qpp::conjugate(const Eigen::MatrixBase<Derived>& A)
TEST(qpp_conjugate, AllTests) {}
/******************************************************************************/
/// BEGIN template<typename Derived> cmat qpp::cosm(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_cosm, AllTests) {}
/******************************************************************************/
/// BEGIN template<typename OutputScalar, typename Derived>
///       dyn_mat <OutputScalar> qpp::cwise(const Eigen::MatrixBase<Derived>& A,
///       OutputScalar (* f)(const typename Derived::Scalar&))
TEST(qpp_cwise, AllTests) {}
/******************************************************************************/
/// BEGIN template<typename Derived> typename Derived::Scalar qpp::det(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_det, AllTests) {}
/******************************************************************************/
/// BEGIN template<typename Derived> dyn_mat<typename Derived::Scalar>
///       qpp::dirsum(const std::initializer_list<Derived>& As)
TEST(qpp_dirsum_initlist, AllTests) {}
/******************************************************************************/
/// BEGIN template<typename T, typename ... Args>
/// dyn_mat<typename T::Scalar> qpp::dirsum(const T& head, const Args& ... tail)
TEST(qpp_dirsum_variadic, AllTests) {}
/******************************************************************************/
/// BEGIN template<typename Derived> dyn_mat<typename Derived::Scalar>
///       qpp::dirsum(const std::vector<Derived>& As)
TEST(qpp_dirsum_vector, AllTests) {}
/******************************************************************************/
/// BEGIN template<typename Derived> dyn_mat<typename Derived::Scalar>
///       qpp::dirsumpow(const Eigen::MatrixBase<Derived>& A, idx n)
TEST(qpp_dirsumpow, AllTests) {}
/******************************************************************************/
/// BEGIN template<typename Derived> std::pair<dyn_col_vect<cplx>, cmat>
///       qpp::eig (const Eigen::MatrixBase<Derived>& A)
TEST(qpp_eig, AllTests) {}
/******************************************************************************/
/// BEGIN template<typename Derived> dyn_col_vect <cplx> qpp::evals(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_evals, AllTests) {}
/******************************************************************************/
/// BEGIN template<typename Derived> cmat qpp::evects(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_evects, AllTests) {}
/******************************************************************************/
/// BEGIN template<typename Derived> cmat qpp::expm(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_expm, AllTests) {}
/******************************************************************************/
/// BEGIN template<typename Derived> cmat qpp::funm(
///       const Eigen::MatrixBase<Derived>& A, cplx (* f)(const cplx&))
TEST(qpp_funm, AllTests) {}
/******************************************************************************/
/// BEGIN template<typename Derived> dyn_mat<typename Derived::Scalar>
///       qpp::grams(const std::initializer_list<Derived>& As)
TEST(qpp_grams_initlist, AllTests) {}
/******************************************************************************/
/// BEGIN template<typename Derived> dyn_mat<typename Derived::Scalar>
///       qpp:::grams(const Eigen::MatrixBase<Derived>& A)
TEST(qpp_grams_matrix, AllTests) {}
/******************************************************************************/
/// BEGIN template<typename Derived> dyn_mat<typename Derived::Scalar>
///       qpp::grams(const std::vector<Derived>& As)
TEST(qpp_grams_vector, AllTests) {}
/******************************************************************************/
/// BEGIN template<typename Derived> std::pair<dyn_col_vect<double>, cmat>
///       qpp::heig(const Eigen::MatrixBase<Derived>& A)
TEST(qpp_heig, AllTests) {}
/******************************************************************************/
/// BEGIN template<typename Derived> dyn_col_vect<double> qpp::hevals(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_hevals, AllTests) {}
/******************************************************************************/
/// BEGIN template<typename Derived> cmat qpp::hevects(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_hevects, AllTests) {}
/******************************************************************************/
/// BEGIN template<typename Derived> dyn_mat<typename Derived::Scalar>
///       qpp::inverse(const Eigen::MatrixBase<Derived>& A)
TEST(qpp_inverse, AllTests) {}
/******************************************************************************/
/// BEGIN template<typename Derived> dyn_mat<typename Derived::Scalar>
///       qpp::kron(const std::initializer_list<Derived>& As)
TEST(qpp_kron_initlist, AllTests) {}
/******************************************************************************/
/// BEGIN template<typename T, typename ... Args> dyn_mat<typename T::Scalar>
///       qpp::kron(const T& head, const Args& ... tail)
TEST(qpp_kron_variadic, AllTests) {}
/******************************************************************************/
/// BEGIN template<typename Derived> dyn_mat<typename Derived::Scalar>
///       qpp::kron(const std::vector<Derived>& As)
TEST(qpp_kron_vector, AllTests) {}
/******************************************************************************/
/// BEGIN template<typename Derived> dyn_mat<typename Derived::Scalar>
///       qpp::kronpow(const Eigen::MatrixBase<Derived>& A, idx n)
TEST(qpp_kronpow, AllTests) {}
/******************************************************************************/
/// BEGIN template<typename Derived> typename Derived::Scalar qpp::logdet(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_logdet, AllTests) {}
/******************************************************************************/
/// BEGIN template<typename Derived> cmat qpp::logm(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_logm, AllTests) {}
/******************************************************************************/
/// BEGIN inline ket qpp::mket(const std::vector<idx>& mask,
///       const std::vector<idx>& dims)
TEST(qpp_mket, AllTests) {}
/******************************************************************************/
/// BEGIN inline ket qpp::mket(const std::vector<idx>& mask, idx d = 2)
TEST(qpp_mket_qubits, AllTests) {}
/******************************************************************************/
/// BEGIN template<char... Bits> ket qpp::operator "" _ket()
TEST(qpp_operator_ket, AllTests) {}
/******************************************************************************/
/// BEGIN template<char... Bits> bra qpp::operator "" _bra()
TEST(qpp_operator_bra, AllTests) {}
/******************************************************************************/
/// BEGIN template<char... Bits> cmat qpp::operator "" _prj()
TEST(qpp_operator_prj, AllTests) {}
/******************************************************************************/
/// BEGIN inline cmat qpp::mprj(const std::vector<idx>& mask,
///       const std::vector<idx>& dims)
TEST(qpp_mprj, AllTests) {}
/******************************************************************************/
/// BEGIN inline cmat qpp::mprj(const std::vector<idx>& mask, idx d = 2)
TEST(qpp_mprj_qubits, AllTests) {}
/******************************************************************************/
/// BEGIN inline idx qpp::multiidx2n(const std::vector<idx>& midx,
///       const std::vector<idx>& dims)
TEST(qpp_multiidx2n, AllTests) {}
/******************************************************************************/
/// BEGIN inline std::vector<idx> qpp::n2multiidx(idx n,
///       const std::vector<idx>& dims)
TEST(qpp_n2multiidx, AllTests) {}
/******************************************************************************/
/// BEGIN template<typename Derived> double qpp::norm(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_norm, AllTests) {}
/******************************************************************************/
/// BEGIN template<typename Derived> dyn_mat<typename Derived::Scalar>
///       qpp::powm(const Eigen::MatrixBase<Derived>& A, idx n)
TEST(qpp_powm, AllTests) {}
/******************************************************************************/
/// BEGIN template<typename Derived> dyn_mat<typename Derived::Scalar>
///       qpp::prj(const Eigen::MatrixBase<Derived>& A)
TEST(qpp_prj, AllTests) {}
/******************************************************************************/
/// BEGIN template<typename Container> typename Container::value_type
///       qpp::prod(const Container& c,
///       typename std::enable_if<is_iterable<Container>::value>::type*
///       = nullptr)
TEST(qpp_prod_container, AllTests) {}
/******************************************************************************/
/// BEGIN template<typename Derived> typename Derived::Scalar qpp::prod(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_prod_matrix, AllTests) {}
/******************************************************************************/
/// BEGIN template<typename InputIterator>
///       typename std::iterator_traits<InputIterator>::value_type
///       qpp::prod(InputIterator first, InputIterator last)
TEST(qpp_prod_iterator, AllTests) {
    std::vector<int> v{1, 2, 3, 4};
    EXPECT_EQ(24, qpp::prod(v.begin(), v.end()));
}
/******************************************************************************/
/// BEGIN template<typename Derived> dyn_mat<typename Derived::Scalar>
///       qpp::reshape(const Eigen::MatrixBase<Derived>& A, idx rows, idx cols)
TEST(qpp_reshape, AllTests) {}
/******************************************************************************/
/// BEGIN template<typename Derived> std::vector<double> qpp::rho2bloch(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_rho2bloch, AllTests) {}
/******************************************************************************/
/// BEGIN template<typename Derived> dyn_col_vect<typename Derived::Scalar>
///       qpp::rho2pure(const Eigen::MatrixBase<Derived>& A)
TEST(qpp_rho2pure, AllTests) {}
/******************************************************************************/
/// BEGIN template<typename Derived> double qpp::schatten(
///       const Eigen::MatrixBase<Derived>& A, double p)
TEST(qpp_schatten, AllTests) {}
/******************************************************************************/
/// BEGIN template<typename Derived> cmat qpp::sinm(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_sinm, AllTests) {}
/******************************************************************************/
/// BEGIN template<typename Derived> cmat qpp::spectralpowm(
///       const Eigen::MatrixBase<Derived>& A, const cplx z)
TEST(qpp_spectralpowm, AllTests) {}
/******************************************************************************/
/// BEGIN template<typename Derived> cmat qpp::sqrtm(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_sqrtm, AllTests) {}
/******************************************************************************/
/// BEGIN template<typename Container> typename Container::value_type
///       qpp::sum(const Container& c,
///       typename std::enable_if<is_iterable<Container>::value>::type*
///       = nullptr)
TEST(qpp_sum_container, AllTests) {}
/// BEGIN template<typename Derived> typename Derived::Scalar qpp::sum(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_sum_matrix, AllTests) {}
/******************************************************************************/
/// BEGIN template<typename InputIterator>
///       typename std::iterator_traits<InputIterator>::value_type
///       qpp::sum(InputIterator first, InputIterator last)
TEST(qpp_sum_iterator, AllTests) {
    std::vector<int> v{0, 1, 2, 3};
    EXPECT_EQ(6, qpp::sum(v.begin(), v.end()));
}
/******************************************************************************/
/// BEGIN template<typename Derived> dyn_col_vect<double> qpp::svals(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_svals, AllTests) {}
/******************************************************************************/
/// BEGIN template<typename Derived>
///       std::tuple<cmat, dyn_col_vect<double>, cmat> qpp::svd(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_svd, AllTests) {}
/******************************************************************************/
/// BEGIN template<typename Derived> cmat qpp::svdU(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_svdU, AllTests) {}
/******************************************************************************/
/// BEGIN template<typename Derived> cmat qpp::svdV(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_svdV, AllTests) {}
/******************************************************************************/
/// BEGIN template<typename Derived> typename Derived::Scalar qpp::trace(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_trace, AllTests) {}
/******************************************************************************/
/// BEGIN template<typename Derived> dyn_mat<typename Derived::Scalar>
///       qpp::transpose(const Eigen::MatrixBase<Derived>& A)
TEST(qpp_transpose, AllTests) {}
/******************************************************************************/
