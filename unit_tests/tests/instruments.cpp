#include "gtest/gtest.h"
#include "qpp.h"

using namespace qpp;

// Unit testing "instruments.h"

/******************************************************************************/
/// BEGIN template<typename Derived> dyn_col_vect<typename Derived::Scalar>
///       qpp::ip(const Eigen::MatrixBase<Derived>& phi,
///       const Eigen::MatrixBase<Derived>& psi,
///       const std::vector<idx>& subsys,
///       const std::vector<idx>& dims)
TEST(qpp_ip, AllTests) {}
/******************************************************************************/
/// BEGIN template<typename Derived> dyn_col_vect<typename Derived::Scalar>
///       qpp::ip(const Eigen::MatrixBase<Derived>& phi,
///       const Eigen::MatrixBase<Derived>& psi,
///       const std::vector<idx>& subsys,
///       idx d = 2)
TEST(qpp_ip_qubits, AllTests) {}
/******************************************************************************/
/// BEGIN template<typename Derived>
///       std::tuple<idx, std::vector<double>, std::vector<cmat>> qpp::measure(
///       const Eigen::MatrixBase<Derived>& A,
///       const cmat& U)
TEST(qpp_measure_full_orthonormal, AllTests) {}
/******************************************************************************/
/// BEGIN template<typename Derived>
///       std::tuple<idx, std::vector<double>, std::vector<cmat>> qpp::measure(
///       const Eigen::MatrixBase<Derived>& A,
///       const cmat& V,
///       const std::vector<idx>& target,
///       const std::vector<idx>& dims)
TEST(qpp_measure_rankone, AllTests) {}
/******************************************************************************/
/// BEGIN template<typename Derived>
///       std::tuple<idx, std::vector<double>, std::vector<cmat>> qpp::measure(
///       const Eigen::MatrixBase<Derived>& A,
///       const cmat& V,
///       const std::vector<idx>& target,
///       idx d = 2)
TEST(qpp_measure_rankone_qubits, AllTests) {}
/******************************************************************************/
/// BEGIN template<typename Derived>
///       std::tuple<idx, std::vector<double>, std::vector<cmat>>
///       qpp::measure(const Eigen::MatrixBase<Derived>& A,
///       const std::initializer_list<cmat>& Ks)
TEST(qpp_measure_full_kraus_initlist, AllTests) {}
/******************************************************************************/
/// BEGIN template<typename Derived>
///       std::tuple<idx, std::vector<double>, std::vector<cmat>>
///       qpp::measure(const Eigen::MatrixBase<Derived>& A,
///       const std::initializer_list<cmat>& Ks,
///       const std::vector<idx>& target,
///       const std::vector<idx>& dims)
TEST(qpp_measure_kraus_initlist, AllTests) {}
/******************************************************************************/
/// BEGIN template<typename Derived>
///       std::tuple<idx, std::vector<double>, std::vector<cmat>>
///       qpp::measure(const Eigen::MatrixBase<Derived>& A,
///       const std::initializer_list<cmat>& Ks,
///       const std::vector<idx>& target,
///       idx d = 2)
TEST(qpp_measure_kraus_initlist_qubits, AllTests) {}
/******************************************************************************/
/// BEGIN template<typename Derived>
///       std::tuple<idx, std::vector<double>, std::vector<cmat>>
///       qpp::measure(const Eigen::MatrixBase<Derived>& A,
///       const std::vector<cmat>& Ks)
TEST(qpp_measure_full_kraus_vector, AllTests) {}
/******************************************************************************/
/// BEGIN template<typename Derived>
///       std::tuple<idx, std::vector<double>, std::vector<cmat>>
///       qpp::measure(const Eigen::MatrixBase<Derived>& A,
///       const std::vector<cmat>& Ks,
///       const std::vector<idx>& target,
///       const std::vector<idx>& dims)
TEST(qpp_measure_kraus_vector, AllTests) {}
/******************************************************************************/
/// BEGIN template<typename Derived>
///       std::tuple<idx, std::vector<double>, std::vector<cmat>>
///       qpp::measure(const Eigen::MatrixBase<Derived>& A,
///       const std::vector<cmat>& Ks,
///       const std::vector<idx>& target,
///       idx d = 2)
TEST(qpp_measure_kraus_vector_qubits, AllTests) {}
/******************************************************************************/
/// BEGIN template<typename Derived> std::tuple<std::vector<idx>, double, cmat>
///       qpp::measure_seq(const Eigen::MatrixBase<Derived>& A,
///       std::vector<idx> target,
///       idx d = 2)
TEST(qpp_measure_seq_qubits, AllTests) {}
/******************************************************************************/
/// BEGIN template<typename Derived> std::tuple<std::vector<idx>, double, cmat>
///       qpp::measure_seq(const Eigen::MatrixBase<Derived>& A,
///       std::vector<idx> target,
///       std::vector<idx> dims)
TEST(qpp_measure_seq, AllTests) {}
/******************************************************************************/
