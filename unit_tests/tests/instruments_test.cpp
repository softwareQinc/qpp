#include <map>
#include <tuple>
#include <vector>

#include "gtest/gtest.h"

#include "qpp/qpp.hpp"

using namespace qpp;

// Unit testing "qpp/instruments.hpp"

/// BEGIN template <typename Derived>
///       dyn_mat<typename Derived::Scalar> discard(
///       const Eigen::MatrixBase<Derived>& A, const std::vector<idx>& target,
///       const std::vector<idx>& dims)
TEST(qpp_discard, Qudits) {
    const std::vector<idx> dims{3, 2};
    const ket psi = mket({2, 1}, dims);
    EXPECT_LT((discard(psi, {0}, dims) - mket({1})).norm(), 1e-12);
    EXPECT_LT((discard(prj(psi), {0}, dims) - mprj({1})).norm(), 1e-12);
    EXPECT_THROW(discard(ket{}, {0}, dims), exception::ZeroSize);
}

/// BEGIN template <typename Derived>
///       dyn_mat<typename Derived::Scalar> discard(
///       const Eigen::MatrixBase<Derived>& A, const std::vector<idx>& target,
///       idx d = 2)
TEST(qpp_discard, Qubits) {
    const ket psi = mket({1, 0, 1});
    EXPECT_LT((discard(psi, {1}) - mket({1, 1})).norm(), 1e-12);
    EXPECT_LT((discard(prj(psi), {1}) - mprj({1, 1})).norm(), 1e-12);
}

/// BEGIN template <typename Derived> dyn_col_vect<typename Derived::Scalar>
///       ip(const Eigen::MatrixBase<Derived>& phi,
///       const Eigen::MatrixBase<Derived>& psi, const std::vector<idx>& subsys,
///       const std::vector<idx>& dims)
TEST(qpp_ip, Qudits) {
    const ket phi = mket({2}, 3);
    const ket psi = mket({2, 1}, std::vector<idx>{3, 2});
    EXPECT_LT((ip(phi, psi, {0}, {3, 2}) - mket({1})).norm(), 1e-12);
    EXPECT_NEAR(0, norm(ip(mket({0}, 3), psi, {0}, {3, 2})), 1e-12);
}

/// BEGIN template <typename Derived> dyn_col_vect<typename Derived::Scalar>
///       ip(const Eigen::MatrixBase<Derived>& phi,
///       const Eigen::MatrixBase<Derived>& psi, const std::vector<idx>& subsys,
///       idx d = 2)
TEST(qpp_ip, Qubits) {
    EXPECT_LT((ip(st.z1, mket({1, 0}), {0}) - st.z0).norm(), 1e-12);
    EXPECT_NEAR(0, norm(ip(st.z0, mket({1, 0}), {0})), 1e-12);
}

/// BEGIN template <typename Derived>
///       std::tuple<idx, std::vector<realT>, std::vector<expr_t<Derived>>>
///       measure(const Eigen::MatrixBase<Derived>& A, const cmat& U)
TEST(qpp_measure, Orthonormal) {
    auto [m, p, states] = measure(st.z1, gt.Id(2));
    EXPECT_EQ(1u, m);
    EXPECT_EQ((std::vector<realT>{0, 1}), p);
    EXPECT_LT((states[m] - st.z1).norm(), 1e-12);
    auto [mr, pr, post] = measure(mprj({0}), gt.Id(2));
    EXPECT_EQ(0u, mr);
    EXPECT_EQ((std::vector<realT>{1, 0}), pr);
    EXPECT_LT((post[mr] - mprj({0})).norm(), 1e-12);
}

/// BEGIN template <typename Derived>
///       std::tuple<idx, std::vector<realT>, std::vector<expr_t<Derived>>>
///       measure(const Eigen::MatrixBase<Derived>& A, const cmat& V,
///       const std::vector<idx>& target, const std::vector<idx>& dims,
///       bool destructive = true)
TEST(qpp_measure, RankOneQudits) {
    const ket psi = mket({2, 1}, std::vector<idx>{3, 2});
    auto [m, p, states] = measure(psi, gt.Id(3), {0}, {3, 2});
    EXPECT_EQ(2u, m);
    EXPECT_EQ((std::vector<realT>{0, 0, 1}), p);
    EXPECT_LT((states[m] - mket({1})).norm(), 1e-12);
    auto [mn, pn, nondestructive] = measure(psi, gt.Id(3), {0}, {3, 2}, false);
    EXPECT_EQ(2u, mn);
    EXPECT_LT((nondestructive[mn] - psi).norm(), 1e-12);
}

/// BEGIN template <typename Derived>
///       std::tuple<idx, std::vector<realT>, std::vector<expr_t<Derived>>>
///       measure(const Eigen::MatrixBase<Derived>& A, const cmat& V,
///       const std::vector<idx>& target, idx d = 2, bool destructive = true)
TEST(qpp_measure, RankOneQubits) {
    const ket psi = mket({1, 0});
    auto [m, p, states] = measure(psi, gt.Id(2), {0});
    EXPECT_EQ(1u, m);
    EXPECT_EQ((std::vector<realT>{0, 1}), p);
    EXPECT_LT((states[m] - st.z0).norm(), 1e-12);
    auto [mn, pn, nondestructive] = measure(psi, gt.Id(2), {0}, 2, false);
    EXPECT_EQ(1u, mn);
    EXPECT_LT((nondestructive[mn] - psi).norm(), 1e-12);
}

/// BEGIN template <typename Derived>
///       std::tuple<idx, std::vector<realT>, std::vector<expr_t<Derived>>>
///       measure(const Eigen::MatrixBase<Derived>& A,
///       const std::initializer_list<cmat>& Ks)
TEST(qpp_measure, KrausInitList) {
    auto [m, p, states] = measure(st.z1, {mprj({0}), mprj({1})});
    EXPECT_EQ(1u, m);
    EXPECT_EQ((std::vector<realT>{0, 1}), p);
    EXPECT_LT((states[m] - st.z1).norm(), 1e-12);
}

/// BEGIN template <typename Derived>
///       std::tuple<idx, std::vector<realT>, std::vector<expr_t<Derived>>>
///       measure(const Eigen::MatrixBase<Derived>& A,
///       const std::initializer_list<cmat>& Ks, const std::vector<idx>& target,
///       const std::vector<idx>& dims, bool destructive = true)
TEST(qpp_measure, SubsysKrausInitList) {
    const ket psi = mket({2, 1}, std::vector<idx>{3, 2});
    auto [m, p, states] = measure(
        psi, {mprj({0}, 3), mprj({1}, 3), mprj({2}, 3)}, {0}, {3, 2}, false);
    EXPECT_EQ(2u, m);
    EXPECT_EQ((std::vector<realT>{0, 0, 1}), p);
    EXPECT_LT((states[m] - psi).norm(), 1e-12);
}

/// BEGIN template <typename Derived>
///       std::tuple<idx, std::vector<realT>, std::vector<expr_t<Derived>>>
///       measure(const Eigen::MatrixBase<Derived>& A,
///       const std::initializer_list<cmat>& Ks, const std::vector<idx>& target,
///       idx d = 2, bool destructive = true)
TEST(qpp_measure, SubsysKrausInitListQubits) {
    const ket psi = mket({1, 0});
    auto [m, p, states] = measure(psi, {mprj({0}), mprj({1})}, {0}, 2, false);
    EXPECT_EQ(1u, m);
    EXPECT_EQ((std::vector<realT>{0, 1}), p);
    EXPECT_LT((states[m] - psi).norm(), 1e-12);
}

/// BEGIN template <typename Derived>
///       std::tuple<idx, std::vector<realT>, std::vector<expr_t<Derived>>>
///       measure(const Eigen::MatrixBase<Derived>& A,
///       const std::vector<cmat>& Ks)
TEST(qpp_measure, KrausVector) {
    const std::vector<cmat> Ks{mprj({0}), mprj({1})};
    auto [m, p, states] = measure(mprj({1}), Ks);
    EXPECT_EQ(1u, m);
    EXPECT_EQ((std::vector<realT>{0, 1}), p);
    EXPECT_LT((states[m] - mprj({1})).norm(), 1e-12);
}

/// BEGIN template <typename Derived>
///       std::tuple<idx, std::vector<realT>, std::vector<expr_t<Derived>>>
///       measure(const Eigen::MatrixBase<Derived>& A,
///       const std::vector<cmat>& Ks, const std::vector<idx>& target,
///       const std::vector<idx>& dims, bool destructive = true)
TEST(qpp_measure, SubsysKrausVector) {
    const std::vector<cmat> Ks{mprj({0}, 3), mprj({1}, 3), mprj({2}, 3)};
    auto [m, p, states] = measure(mprj({2, 1}, {3, 2}), Ks, {0}, {3, 2});
    EXPECT_EQ(2u, m);
    EXPECT_EQ((std::vector<realT>{0, 0, 1}), p);
    EXPECT_LT((states[m] - mprj({1})).norm(), 1e-12);
}

/// BEGIN template <typename Derived>
///       std::tuple<idx, std::vector<realT>, std::vector<expr_t<Derived>>>
///       measure(const Eigen::MatrixBase<Derived>& A,
///       const std::vector<cmat>& Ks, const std::vector<idx>& target,
///       idx d = 2)
TEST(qpp_measure, SubsysKrausVectorQubits) {
    const std::vector<cmat> Ks{mprj({0}), mprj({1})};
    const ket psi = mket({1, 0});
    auto [m, p, states] = measure(psi, Ks, {0}, 2, false);
    EXPECT_EQ(1u, m);
    EXPECT_EQ((std::vector<realT>{0, 1}), p);
    EXPECT_LT((states[m] - psi).norm(), 1e-12);
}

/// BEGIN template <typename Derived>
///       std::tuple<std::vector<idx>, std::vector<realT>, expr_t<Derived>>
///       measure_seq(const Eigen::MatrixBase<Derived>& A,
///       const std::vector<idx>& target, idx d = 2, bool destructive = true)
TEST(qpp_measure_seq, Qubits) {
    const ket psi = mket({1, 0, 1});
    auto [m, p, state] = measure_seq(psi, {2, 0}, 2, true);
    EXPECT_EQ((std::vector<idx>{1, 1}), m);
    EXPECT_EQ((std::vector<realT>{1, 1}), p);
    EXPECT_LT((state - st.z0).norm(), 1e-12);
    auto [mn, pn, unchanged] = measure_seq(psi, {2, 0}, 2, false);
    EXPECT_EQ((std::vector<idx>{1, 1}), mn);
    EXPECT_LT((unchanged - psi).norm(), 1e-12);
}

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

/// BEGIN template <typename Derived>
///       dyn_mat<typename Derived::Scalar> reset(
///       const Eigen::MatrixBase<Derived>& A, const std::vector<idx>& target,
///       const std::vector<idx>& dims)
TEST(qpp_reset, Qudits) {
    const ket psi = mket({2, 1}, std::vector<idx>{3, 2});
    EXPECT_LT((reset(psi, {0}, {3, 2}) - mket({0, 1}, {3, 2})).norm(), 1e-12);
    EXPECT_LT((reset(prj(psi), {0}, {3, 2}) - mprj({0, 1}, {3, 2})).norm(),
              1e-12);
}

/// BEGIN template <typename Derived>
///       dyn_mat<typename Derived::Scalar> reset(
///       const Eigen::MatrixBase<Derived>& A, const std::vector<idx>& target,
///       idx d = 2)
TEST(qpp_reset, Qubits) {
    const ket psi = mket({1, 1});
    EXPECT_LT((reset(psi, {1}) - mket({1, 0})).norm(), 1e-12);
    EXPECT_LT((reset(prj(psi), {1}) - mprj({1, 0})).norm(), 1e-12);
}

/// BEGIN template <typename Derived> std::tuple<std::vector<idx>, realT>
///       sample(const Eigen::MatrixBase<Derived>& A,
///       const std::vector<idx>& target, const std::vector<idx>& dims)
TEST(qpp_sample, Qudits) {
    const ket psi = mket({2, 1}, std::vector<idx>{3, 2});
    auto [m, p] = sample(psi, {1, 0}, {3, 2});
    EXPECT_EQ((std::vector<idx>{1, 2}), m);
    EXPECT_NEAR(1, p, 1e-12);
    auto [mr, pr] = sample(prj(psi), {0}, {3, 2});
    EXPECT_EQ((std::vector<idx>{2}), mr);
    EXPECT_NEAR(1, pr, 1e-12);
}

/// BEGIN template <typename Derived> std::tuple<std::vector<idx>, realT>
///       sample(const Eigen::MatrixBase<Derived>& A,
///       const std::vector<idx>& target, idx d = 2)
TEST(qpp_sample, Qubits) {
    const ket psi = mket({1, 0, 1});
    auto [m, p] = sample(psi, {2, 0});
    EXPECT_EQ((std::vector<idx>{1, 1}), m);
    EXPECT_NEAR(1, p, 1e-12);
}

/// BEGIN template <typename Derived> std::map<std::vector<idx>, idx>
///       sample(idx num_samples, const Eigen::MatrixBase<Derived>& A,
///       const std::vector<idx>& target, const std::vector<idx>& dims)
TEST(qpp_sample, MultipleRunsQudits) {
    const ket psi = mket({2, 1}, std::vector<idx>{3, 2});
    const auto counts = sample(100, psi, {0, 1}, {3, 2});
    ASSERT_EQ(1u, counts.size());
    EXPECT_EQ(100u, counts.at(std::vector<idx>{2, 1}));
}

/// BEGIN template <typename Derived> std::map<std::vector<idx>, idx>
///       sample(idx num_samples,const Eigen::MatrixBase<Derived>& A,
///       const std::vector<idx>& target, idx d = 2)
TEST(qpp_sample, MultipleRunsQubits) {
    const auto counts = sample(100, mprj({1, 0}), {0, 1});
    ASSERT_EQ(1u, counts.size());
    EXPECT_EQ(100u, counts.at(std::vector<idx>{1, 0}));
}
