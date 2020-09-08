#include "gtest/gtest.h"
#include "qpp.h"

using namespace qpp;

// Unit testing "instruments.hpp"

/******************************************************************************/
/// BEGIN template <typename Derived>
///       dyn_mat<typename Derived::Scalar> qpp::discard(
///       const Eigen::MatrixBase<Derived>& A, const std::vector<idx>& target,
///       const std::vector<idx>& dims)
TEST(qpp_discard, Qudits) {}
/******************************************************************************/
/// BEGIN template <typename Derived>
///       dyn_mat<typename Derived::Scalar> qpp::discard(
///       const Eigen::MatrixBase<Derived>& A, const std::vector<idx>& target,
///       idx d = 2)
TEST(qpp_discard, Qubits) {}
/******************************************************************************/
/// BEGIN template <typename Derived> dyn_col_vect<typename Derived::Scalar>
///       qpp::ip(const Eigen::MatrixBase<Derived>& phi,
///       const Eigen::MatrixBase<Derived>& psi, const std::vector<idx>& subsys,
///       const std::vector<idx>& dims)
TEST(qpp_ip, Qudits) {}
/******************************************************************************/
/// BEGIN template <typename Derived> dyn_col_vect<typename Derived::Scalar>
///       qpp::ip(const Eigen::MatrixBase<Derived>& phi,
///       const Eigen::MatrixBase<Derived>& psi, const std::vector<idx>& subsys,
///       idx d = 2)
TEST(qpp_ip, Qubits) {}
/******************************************************************************/
/// BEGIN template <typename Derived>
///       std::tuple<idx, std::vector<double>, std::vector<cmat>> qpp::measure(
///       const Eigen::MatrixBase<Derived>& A, const cmat& U)
TEST(qpp_measure, Orthonormal) {}
/******************************************************************************/
/// BEGIN template <typename Derived>
///       std::tuple<idx, std::vector<double>, std::vector<cmat>> qpp::measure(
///       const Eigen::MatrixBase<Derived>& A, const cmat& V,
///       const std::vector<idx>& target, const std::vector<idx>& dims,
///       bool destructive = true)
TEST(qpp_measure, RankOneQudits) {}
/******************************************************************************/
/// BEGIN template <typename Derived>
///       std::tuple<idx, std::vector<double>, std::vector<cmat>> qpp::measure(
///       const Eigen::MatrixBase<Derived>& A, const cmat& V,
///       const std::vector<idx>& target, idx d = 2, bool destructive = true)
TEST(qpp_measure, RankOneQubits) {}
/******************************************************************************/
/// BEGIN template <typename Derived>
///       std::tuple<idx, std::vector<double>, std::vector<cmat>>
///       qpp::measure(const Eigen::MatrixBase<Derived>& A,
///       const std::initializer_list<cmat>& Ks)
TEST(qpp_measure, KrausInitList) {}
/******************************************************************************/
/// BEGIN template <typename Derived>
///       std::tuple<idx, std::vector<double>, std::vector<cmat>>
///       qpp::measure(const Eigen::MatrixBase<Derived>& A,
///       const std::initializer_list<cmat>& Ks, const std::vector<idx>& target,
///       const std::vector<idx>& dims, bool destructive = true)
TEST(qpp_measure, SubsysKrausInitList) {}
/******************************************************************************/
/// BEGIN template <typename Derived>
///       std::tuple<idx, std::vector<double>, std::vector<cmat>>
///       qpp::measure(const Eigen::MatrixBase<Derived>& A,
///       const std::initializer_list<cmat>& Ks, const std::vector<idx>& target,
///       idx d = 2, bool destructive = true)
TEST(qpp_measure, SubsysKrausInitListQubits) {}
/******************************************************************************/
/// BEGIN template <typename Derived>
///       std::tuple<idx, std::vector<double>, std::vector<cmat>>
///       qpp::measure(const Eigen::MatrixBase<Derived>& A,
///       const std::vector<cmat>& Ks)
TEST(qpp_measure, KrausVector) {}
/******************************************************************************/
/// BEGIN template <typename Derived>
///       std::tuple<idx, std::vector<double>, std::vector<cmat>>
///       qpp::measure(const Eigen::MatrixBase<Derived>& A,
///       const std::vector<cmat>& Ks, const std::vector<idx>& target,
///       const std::vector<idx>& dims, bool destructive = true)
TEST(qpp_measure, SubsysKrausVector) {}
/******************************************************************************/
/// BEGIN template <typename Derived>
///       std::tuple<idx, std::vector<double>, std::vector<cmat>>
///       qpp::measure(const Eigen::MatrixBase<Derived>& A,
///       const std::vector<cmat>& Ks, const std::vector<idx>& target,
///       idx d = 2)
TEST(qpp_measure, SubsysKrausVectorQubits) {}
/******************************************************************************/
/// BEGIN template <typename Derived> std::tuple<std::vector<idx>, double, cmat>
///       qpp::measure_seq(const Eigen::MatrixBase<Derived>& A,
///       const std::vector<idx>& target, idx d = 2, bool destructive = true)
TEST(qpp_measure_seq, Qubits) {}
/******************************************************************************/
/// BEGIN template <typename Derived> std::tuple<std::vector<idx>, double, cmat>
///       qpp::measure_seq(const Eigen::MatrixBase<Derived>& A,
///       std::vector<idx> target, std::vector<idx> dims,
///       bool destructive = true)
TEST(qpp_measure_seq, Qudits) {}
/******************************************************************************/
/// BEGIN template <typename Derived>
///       dyn_mat<typename Derived::Scalar> qpp::reset(
///       const Eigen::MatrixBase<Derived>& A, const std::vector<idx>& target,
///       const std::vector<idx>& dims)
TEST(qpp_reset, Qudits) {}
/******************************************************************************/
/// BEGIN template <typename Derived>
///       dyn_mat<typename Derived::Scalar> qpp::reset(
///       const Eigen::MatrixBase<Derived>& A, const std::vector<idx>& target,
///       idx d = 2)
TEST(qpp_reset, Qubits) {}
/******************************************************************************/
