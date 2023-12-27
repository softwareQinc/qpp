#include "gtest/gtest.h"

#include "qpp/qpp.h"

using namespace qpp;

// Unit testing "instruments.hpp"

/******************************************************************************/
/// BEGIN template <typename Derived>
///       dyn_mat<typename Derived::Scalar> discard(
///       const Eigen::MatrixBase<Derived>& A, const std::vector<idx>& target,
///       const std::vector<idx>& dims)
TEST(qpp_discard, Qudits) {}
/******************************************************************************/
/// BEGIN template <typename Derived>
///       dyn_mat<typename Derived::Scalar> discard(
///       const Eigen::MatrixBase<Derived>& A, const std::vector<idx>& target,
///       idx d = 2)
TEST(qpp_discard, Qubits) {}
/******************************************************************************/
/// BEGIN template <typename Derived> dyn_col_vect<typename Derived::Scalar>
///       ip(const Eigen::MatrixBase<Derived>& phi,
///       const Eigen::MatrixBase<Derived>& psi, const std::vector<idx>& subsys,
///       const std::vector<idx>& dims)
TEST(qpp_ip, Qudits) {}
/******************************************************************************/
/// BEGIN template <typename Derived> dyn_col_vect<typename Derived::Scalar>
///       ip(const Eigen::MatrixBase<Derived>& phi,
///       const Eigen::MatrixBase<Derived>& psi, const std::vector<idx>& subsys,
///       idx d = 2)
TEST(qpp_ip, Qubits) {}
/******************************************************************************/
/// BEGIN template <typename Derived>
///       std::tuple<idx, std::vector<realT>, std::vector<expr_t<Derived>>>
///       measure(const Eigen::MatrixBase<Derived>& A, const cmat& U)
TEST(qpp_measure, Orthonormal) {}
/******************************************************************************/
/// BEGIN template <typename Derived>
///       std::tuple<idx, std::vector<realT>, std::vector<expr_t<Derived>>>
///       measure(const Eigen::MatrixBase<Derived>& A, const cmat& V,
///       const std::vector<idx>& target, const std::vector<idx>& dims,
///       bool destructive = true)
TEST(qpp_measure, RankOneQudits) {}
/******************************************************************************/
/// BEGIN template <typename Derived>
///       std::tuple<idx, std::vector<realT>, std::vector<expr_t<Derived>>>
///       measure(const Eigen::MatrixBase<Derived>& A, const cmat& V,
///       const std::vector<idx>& target, idx d = 2, bool destructive = true)
TEST(qpp_measure, RankOneQubits) {}
/******************************************************************************/
/// BEGIN template <typename Derived>
///       std::tuple<idx, std::vector<realT>, std::vector<expr_t<Derived>>>
///       measure(const Eigen::MatrixBase<Derived>& A,
///       const std::initializer_list<cmat>& Ks)
TEST(qpp_measure, KrausInitList) {}
/******************************************************************************/
/// BEGIN template <typename Derived>
///       std::tuple<idx, std::vector<realT>, std::vector<expr_t<Derived>>>
///       measure(const Eigen::MatrixBase<Derived>& A,
///       const std::initializer_list<cmat>& Ks, const std::vector<idx>& target,
///       const std::vector<idx>& dims, bool destructive = true)
TEST(qpp_measure, SubsysKrausInitList) {}
/******************************************************************************/
/// BEGIN template <typename Derived>
///       std::tuple<idx, std::vector<realT>, std::vector<expr_t<Derived>>>
///       measure(const Eigen::MatrixBase<Derived>& A,
///       const std::initializer_list<cmat>& Ks, const std::vector<idx>& target,
///       idx d = 2, bool destructive = true)
TEST(qpp_measure, SubsysKrausInitListQubits) {}
/******************************************************************************/
/// BEGIN template <typename Derived>
///       std::tuple<idx, std::vector<realT>, std::vector<expr_t<Derived>>>
///       measure(const Eigen::MatrixBase<Derived>& A,
///       const std::vector<cmat>& Ks)
TEST(qpp_measure, KrausVector) {}
/******************************************************************************/
/// BEGIN template <typename Derived>
///       std::tuple<idx, std::vector<realT>, std::vector<expr_t<Derived>>>
///       measure(const Eigen::MatrixBase<Derived>& A,
///       const std::vector<cmat>& Ks, const std::vector<idx>& target,
///       const std::vector<idx>& dims, bool destructive = true)
TEST(qpp_measure, SubsysKrausVector) {}
/******************************************************************************/
/// BEGIN template <typename Derived>
///       std::tuple<idx, std::vector<realT>, std::vector<expr_t<Derived>>>
///       measure(const Eigen::MatrixBase<Derived>& A,
///       const std::vector<cmat>& Ks, const std::vector<idx>& target,
///       idx d = 2)
TEST(qpp_measure, SubsysKrausVectorQubits) {}
/******************************************************************************/
/// BEGIN template <typename Derived>
///       std::tuple<std::vector<idx>, std::vector<realT>, expr_t<Derived>>
///       measure_seq(const Eigen::MatrixBase<Derived>& A,
///       const std::vector<idx>& target, idx d = 2, bool destructive = true)
TEST(qpp_measure_seq, Qubits) {}
/******************************************************************************/
/// BEGIN template <typename Derived>
///       std::tuple<std::vector<idx>, std::vector<realT>, expr_t<Derived>>
///       measure_seq(const Eigen::MatrixBase<Derived>& A,
///       std::vector<idx> target, std::vector<idx> dims,
///       bool destructive = true)
TEST(qpp_measure_seq, Qudits) {
    std::vector<idx> dims = {2, 2, 3, 2, 3, 2};
    std::vector<idx> mask = {0, 1, 2, 1, 2, 0};
    std::vector<idx> target = {4, 2, 3, 0};
    ket psi = mket(mask, dims);
    cmat rho = mprj(mask, dims);

    // destructive
    auto [ms_d, ps_d, psi_final_d] = measure_seq(psi, target, dims, true);
    std::vector<idx> expected_ms_d = {2, 2, 1, 0};
    std::vector<realT> expected_ps_d = {1, 1, 1, 1};
    ket expected_psi_final_d = mket({1, 0}, {2, 2});
    EXPECT_EQ(expected_ms_d, ms_d);
    EXPECT_EQ(expected_ps_d, ps_d);
    EXPECT_NEAR(norm(psi_final_d - expected_psi_final_d), 0, 1e-5);

    // destructive, density matrix
    cmat rho_final_d;
    std::tie(ms_d, ps_d, rho_final_d) = measure_seq(rho, target, dims, true);
    expected_ms_d = {2, 2, 1, 0};
    expected_ps_d = {1, 1, 1, 1};
    cmat expected_rho_final_d = mprj({1, 0}, {2, 2});
    EXPECT_EQ(expected_ms_d, ms_d);
    EXPECT_EQ(expected_ps_d, ps_d);
    EXPECT_NEAR(norm(rho_final_d - expected_rho_final_d), 0, 1e-5);

    // non-destructive, non-repetitive
    auto [ms_nd, ps_nd, psi_final_nd] = measure_seq(psi, target, dims, false);
    std::vector<idx> expected_ms_nd = {2, 2, 1, 0};
    std::vector<realT> expected_ps_nd = {1, 1, 1, 1};
    ket expected_psi_final_nd = psi;
    EXPECT_EQ(expected_ms_nd, ms_nd);
    EXPECT_EQ(expected_ps_nd, ps_nd);
    EXPECT_NEAR(norm(psi_final_nd - expected_psi_final_nd), 0, 1e-5);

    // non-destructive, repetitive, same target (different order)
    target = {3, 0, 4, 2};
    std::tie(ms_nd, ps_nd, psi_final_nd) =
        measure_seq(psi, target, dims, false);
    expected_ms_nd = {1, 0, 2, 2};
    expected_ps_nd = {1, 1, 1, 1};
    expected_psi_final_nd = psi;
    EXPECT_EQ(expected_ms_nd, ms_nd);
    EXPECT_EQ(expected_ps_nd, ps_nd);
    EXPECT_NEAR(norm(psi_final_nd - expected_psi_final_nd), 0, 1e-5);

    // non-destructive, repetitive, same target (different order), density
    // matrix
    target = {3, 0, 4, 2};
    cmat rho_final_nd;
    std::tie(ms_nd, ps_nd, rho_final_nd) =
        measure_seq(rho, target, dims, false);
    expected_ms_nd = {1, 0, 2, 2};
    expected_ps_nd = {1, 1, 1, 1};
    cmat expected_rho_final_nd = rho;
    EXPECT_EQ(expected_ms_nd, ms_nd);
    EXPECT_EQ(expected_ps_nd, ps_nd);
    EXPECT_NEAR(norm(rho_final_nd - expected_rho_final_nd), 0, 1e-5);
}
/******************************************************************************/
/// BEGIN template <typename Derived>
///       dyn_mat<typename Derived::Scalar> reset(
///       const Eigen::MatrixBase<Derived>& A, const std::vector<idx>& target,
///       const std::vector<idx>& dims)
TEST(qpp_reset, Qudits) {}
/******************************************************************************/
/// BEGIN template <typename Derived>
///       dyn_mat<typename Derived::Scalar> reset(
///       const Eigen::MatrixBase<Derived>& A, const std::vector<idx>& target,
///       idx d = 2)
TEST(qpp_reset, Qubits) {}
/******************************************************************************/
/// BEGIN template <typename Derived> std::vector<idx>
///       sample(const Eigen::MatrixBase<Derived>& A,
///       const std::vector<idx>& target, const std::vector<idx>& dims)
TEST(qpp_sample, Qudits) {}
/******************************************************************************/
/// BEGIN template <typename Derived> std::vector<idx>
///       sample(const Eigen::MatrixBase<Derived>& A,
///       const std::vector<idx>& target, idx d = 2)
TEST(qpp_sample, Qubits) {}
/******************************************************************************/
/// BEGIN template <typename Derived> std::map<std::vector<idx>, idx>
///       sample(idx num_samples, const Eigen::MatrixBase<Derived>& A,
///       const std::vector<idx>& target, const std::vector<idx>& dims)
TEST(qpp_sample, MultipleQudits) {}
/******************************************************************************/
/// BEGIN template <typename Derived> std::map<std::vector<idx>, idx>
///       sample(idx num_samples,const Eigen::MatrixBase<Derived>& A,
///       const std::vector<idx>& target, idx d = 2)
TEST(qpp_sample, MultipleQubits) {}
/******************************************************************************/
