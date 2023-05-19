#include <vector>
#include "gtest/gtest.h"
#include "qpp.h"

using namespace qpp;

// Unit testing "functions.hpp"

/******************************************************************************/
/// BEGIN template <char... Bits> bra operator"" _bra()
TEST(qpp_literals_operator_bra, AllTests) {}
/******************************************************************************/
/// BEGIN template <char... Bits> ket operator"" _ket()
TEST(qpp_literals_operator_ket, AllTests) {}
/******************************************************************************/
/// BEGIN template <char... Bits> cmat operator"" _prj()
TEST(qpp_literals_operator_prj, AllTests) {}
/******************************************************************************/
/// BEGIN template <typename Derived> cmat absm(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_absm, AllTests) {}
/******************************************************************************/
/// BEGIN template <typename Container> std::vector<realT> abssq(
///       const Container& c,
///       typename std::enable_if<is_iterable<Container>::value>::type*
///       = nullptr)
TEST(qpp_abssq, Container) {}
/******************************************************************************/
/// BEGIN template <typename Derived> std::vector<realT> abssq(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_abssq, Matrix) {}
/******************************************************************************/
/// BEGIN template <typename InputIterator> std::vector<realT> abssq(
///       InputIterator first, InputIterator last)
TEST(qpp_abssq, Iterator) {}
/******************************************************************************/
/// BEGIN template <typename Derived> dyn_mat<typename Derived::Scalar>
///       adjoint(const Eigen::MatrixBase<Derived>& A)
TEST(qpp_adjoint, AllTests) {}
/******************************************************************************/
/// BEGIN template <typename Derived1, typename Derived2>
///       dyn_mat<typename Derived1::Scalar> anticomm(
///       const Eigen::MatrixBase<Derived1>& A,
///       const Eigen::MatrixBase<Derived2>& B)
TEST(qpp_anticomm, AllTests) {}
/******************************************************************************/
/// BEGIN inline cmat bloch2rho(const std::vector<realT>& r)
TEST(qpp_bloch2rho, AllTests) {}
/******************************************************************************/
/// BEGIN template <typename Derived1, typename Derived2>
///       dyn_mat<typename Derived1::Scalar> comm(
///       const Eigen::MatrixBase<Derived1>& A,
///       const Eigen::MatrixBase<Derived2>& B)
TEST(qpp_comm, AllTests) {}
/******************************************************************************/
/// BEGIN inline std::vector<idx> complement(std::vector<idx> subsys, idx n)
TEST(qpp_complement, AllTests) {}
/******************************************************************************/
/// BEGIN template <typename Derived> dyn_mat<typename Derived::Scalar>
///       conjugate(const Eigen::MatrixBase<Derived>& A)
TEST(qpp_conjugate, AllTests) {}
/******************************************************************************/
/// BEGIN template <typename Derived> cmat cosm(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_cosm, AllTests) {}
/******************************************************************************/
/// BEGIN template <typename OutputScalar, typename Derived>
///       dyn_mat <OutputScalar> cwise(const Eigen::MatrixBase<Derived>& A,
///       OutputScalar (*f)(typename Derived::Scalar))
TEST(qpp_cwise, AllTests) {}
/******************************************************************************/
/// BEGIN template <typename Derived> typename Derived::Scalar det(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_det, AllTests) {}
/******************************************************************************/
/// BEGIN template <typename Derived> dyn_mat<typename Derived::Scalar>
///       dirsum(const std::initializer_list<Derived>& As)
TEST(qpp_dirsum, InitList) {}
/******************************************************************************/
/// BEGIN template <typename T, typename ... Args>
/// dyn_mat<typename T::Scalar> dirsum(const T& head, const Args& ... tail)
TEST(qpp_dirsum, Variadic) {}
/******************************************************************************/
/// BEGIN template <typename Derived> dyn_mat<typename Derived::Scalar>
///       dirsum(const std::vector<Derived>& As)
TEST(qpp_dirsum, Vector) {}
/******************************************************************************/
/// BEGIN template <typename Derived> dyn_mat<typename Derived::Scalar>
///       dirsumpow(const Eigen::MatrixBase<Derived>& A, idx n)
TEST(qpp_dirsumpow, AllTests) {}
/******************************************************************************/
/// BEGIN template <typename Derived> std::pair<dyn_col_vect<cplx>, cmat>
///       eig (const Eigen::MatrixBase<Derived>& A)
TEST(qpp_eig, AllTests) {}
/******************************************************************************/
/// BEGIN template <typename Derived> dyn_col_vect <cplx> evals(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_evals, AllTests) {}
/******************************************************************************/
/// BEGIN template <typename Derived> cmat evects(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_evects, AllTests) {}
/******************************************************************************/
/// BEGIN template <typename Derived> cmat expm(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_expm, AllTests) {}
/******************************************************************************/
/// BEGIN template <typename Derived> cmat funm(
///       const Eigen::MatrixBase<Derived>& A, cplx (*f)(const cplx&))
TEST(qpp_funm, AllTests) {}
/******************************************************************************/
/// BEGIN template <typename Derived> dyn_mat<typename Derived::Scalar>
///       grams(const std::initializer_list<Derived>& As)
TEST(qpp_grams, InitList) {}
/******************************************************************************/
/// BEGIN template <typename Derived> dyn_mat<typename Derived::Scalar>
///       :grams(const Eigen::MatrixBase<Derived>& A)
TEST(qpp_grams, Matrix) {}
/******************************************************************************/
/// BEGIN template <typename Derived> dyn_mat<typename Derived::Scalar>
///       grams(const std::vector<Derived>& As)
TEST(qpp_grams, Vector) {}
/******************************************************************************/
/// BEGIN template <typename Derived> std::size_t hash_eigen(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_hash_eigen, AllTests) {}
/// BEGIN template <typename Derived> std::pair<dyn_col_vect<realT>, cmat>
///       heig(const Eigen::MatrixBase<Derived>& A)
TEST(qpp_heig, AllTests) {}
/******************************************************************************/
/// BEGIN template <typename Derived> dyn_col_vect<realT> hevals(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_hevals, AllTests) {}
/******************************************************************************/
/// BEGIN template <typename Derived> cmat hevects(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_hevects, AllTests) {}
/******************************************************************************/
/// BEGIN template <typename Derived> dyn_mat<typename Derived::Scalar>
///       inverse(const Eigen::MatrixBase<Derived>& A)
TEST(qpp_inverse, AllTests) {}
/******************************************************************************/
/// BEGIN template <typename Derived> dyn_mat<typename Derived::Scalar>
///       kron(const std::initializer_list<Derived>& As)
TEST(qpp_kron, InitList) {}
/******************************************************************************/
/// BEGIN template <typename T, typename ... Args> dyn_mat<typename T::Scalar>
///       kron(const T& head, const Args& ... tail)
TEST(qpp_kron, Variadic) {}
/******************************************************************************/
/// BEGIN template <typename Derived> dyn_mat<typename Derived::Scalar>
///       kron(const std::vector<Derived>& As)
TEST(qpp_kron, Vector) {}
/******************************************************************************/
/// BEGIN template <typename Derived> dyn_mat<typename Derived::Scalar>
///       kronpow(const Eigen::MatrixBase<Derived>& A, idx n)
TEST(qpp_kronpow, AllTests) {}
/******************************************************************************/
/// BEGIN template <typename Derived> typename Derived::Scalar logdet(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_logdet, AllTests) {}
/******************************************************************************/
/// BEGIN template <typename Derived> cmat logm(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_logm, AllTests) {}
/******************************************************************************/
/// BEGIN inline ket mket(const std::vector<idx>& mask,
///       const std::vector<idx>& dims)
TEST(qpp_mket, Qudits) {}
/******************************************************************************/
/// BEGIN inline ket mket(const std::vector<idx>& mask, idx d = 2)
TEST(qpp_mket, Qubits) {}
/******************************************************************************/
/// BEGIN inline cmat mprj(const std::vector<idx>& mask,
///       const std::vector<idx>& dims)
TEST(qpp_mprj, Qudits) {}
/******************************************************************************/
/// BEGIN inline cmat mprj(const std::vector<idx>& mask, idx d = 2)
TEST(qpp_mprj, Qubits) {}
/******************************************************************************/
/// BEGIN template <typename V, typename T = V, typename U = T>
///       T multiidx2n(const std::vector<V>& midx, const std::vector<U>& dims)
TEST(qpp_multiidx2n, AllTests) {}
/******************************************************************************/
/// BEGIN <typename T, typename U = T, typename V = T> n2multiidx(T n,
///       const std::vector<U>& dims)
TEST(qpp_n2multiidx, AllTests) {}
/******************************************************************************/
/// BEGIN template <typename Derived> realT norm(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_norm, AllTests) {}
/******************************************************************************/
/// BEGIN template <typename Derived> dyn_mat<typename Derived::Scalar>
///       normalize(const Eigen::MatrixBase<Derived>& A)
TEST(qpp_normalize, AllTests) {}
/******************************************************************************/
/// BEGIN template <typename Derived> dyn_mat<typename Derived::Scalar>
///       powm(const Eigen::MatrixBase<Derived>& A, idx n)
TEST(qpp_powm, AllTests) {}
/******************************************************************************/
/// BEGIN template <typename Derived> dyn_mat<typename Derived::Scalar>
///       prj(const Eigen::MatrixBase<Derived>& A)
TEST(qpp_prj, AllTests) {}
/******************************************************************************/
/// BEGIN template <typename Container> typename Container::value_type
///       prod(const Container& c,
///       typename std::enable_if<is_iterable<Container>::value>::type*
///       = nullptr)
TEST(qpp_prod, Container) {}
/******************************************************************************/
/// BEGIN template <typename T> T prod(const std::initializer_list<T>& Ts)
TEST(qpp_prod, InitList) { EXPECT_EQ(24, prod({1, 2, 3, 4})); }
/******************************************************************************/
/// BEGIN template <typename Derived> typename Derived::Scalar prod(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_prod, Matrix) {}
/******************************************************************************/
/// BEGIN template <typename InputIterator, typename value_type =
///       std::decay_t<decltype(*std::declval<InputIterator>())>>
///       value_type prod(InputIterator first, InputIterator last)
TEST(qpp_prod, Iterator) {
    std::vector<int> v{1, 2, 3, 4};
    EXPECT_EQ(24, prod(v.begin(), v.end()));
}
/******************************************************************************/
/// BEGIN template <typename Derived> dyn_mat<typename Derived::Scalar>
///       reshape(const Eigen::MatrixBase<Derived>& A, idx rows, idx cols)
TEST(qpp_reshape, AllTests) {}
/******************************************************************************/
/// BEGIN template <typename Derived> std::vector<realT> rho2bloch(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_rho2bloch, AllTests) {}
/******************************************************************************/
/// BEGIN template <typename Derived> dyn_col_vect<typename Derived::Scalar>
///       rho2pure(const Eigen::MatrixBase<Derived>& A)
TEST(qpp_rho2pure, AllTests) {}
/******************************************************************************/
/// BEGIN template <typename Derived> realT schatten(
///       const Eigen::MatrixBase<Derived>& A, realT p)
TEST(qpp_schatten, AllTests) {}
/******************************************************************************/
/// BEGIN template <typename Derived> cmat sinm(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_sinm, AllTests) {}
/******************************************************************************/
/// BEGIN template <typename Derived> cmat spectralpowm(
///       const Eigen::MatrixBase<Derived>& A, const cplx z)
TEST(qpp_spectralpowm, AllTests) {}
/******************************************************************************/
/// BEGIN template <typename Derived> cmat sqrtm(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_sqrtm, AllTests) {}
/******************************************************************************/
/// BEGIN template <typename Container> typename Container::value_type
///       sum(const Container& c,
///       typename std::enable_if<is_iterable<Container>::value>::type*
///       = nullptr)
TEST(qpp_sum, Container) {}
/******************************************************************************/
/// BEGIN template <typename T> T sum(const std::initializer_list<T>& Ts)
TEST(qpp_sum, InitList) { EXPECT_EQ(6, sum({0, 1, 2, 3})); }
/// BEGIN template <typename Derived> typename Derived::Scalar sum(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_sum, Matrix) {}
/******************************************************************************/
/// BEGIN template <typename InputIterator, typename value_type =
///       std::decay_t<decltype(*std::declval<InputIterator>())>>
///      value_type sum(InputIterator first, InputIterator last)
TEST(qpp_sum, Iterator) {
    std::vector<int> v{0, 1, 2, 3};
    EXPECT_EQ(6, sum(v.begin(), v.end()));
}
/******************************************************************************/
/// BEGIN template <typename Derived> dyn_col_vect<realT> svals(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_svals, AllTests) {}
/******************************************************************************/
/// BEGIN template <typename Derived>
///       std::tuple<cmat, dyn_col_vect<realT>, cmat> svd(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_svd, AllTests) {}
/******************************************************************************/
/// BEGIN template <typename Derived> cmat svdU(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_svdU, AllTests) {}
/******************************************************************************/
/// BEGIN template <typename Derived> cmat svdV(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_svdV, AllTests) {}
/******************************************************************************/
/// BEGIN template <typename Derived> typename Derived::Scalar trace(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_trace, AllTests) {}
/******************************************************************************/
/// BEGIN template <typename Derived> dyn_mat<typename Derived::Scalar>
///       transpose(const Eigen::MatrixBase<Derived>& A)
TEST(qpp_transpose, AllTests) {}
/******************************************************************************/
/// BEGIN template <typename Derived> std::vector<idx> zket2dits(
///       const Eigen::MatrixBase<Derived>& psi, const std::vector<idx>& dims,
///       realT precision = 1e-12)
TEST(qpp_zket2dits, Qudits) {}
/******************************************************************************/
/// BEGIN template <typename Derived> std::vector<idx> zket2dits(
///       const Eigen::MatrixBase<Derived>& psi, idx d = 2,
///       realT precision = 1e-12)
TEST(qpp_zket2dits, Qubits) {}
/******************************************************************************/
