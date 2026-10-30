#include <cmath>
#include <vector>

#include "gtest/gtest.h"

#include "qpp/qpp.hpp"

using namespace qpp;

// Unit testing "qpp/operations.hpp"

/// BEGIN template <typename Derived1, typename Derived2> expr_t<Derived1>
///       apply(const Eigen::MatrixBase<Derived1>& state,
///       const Eigen::MatrixBase<Derived2>& A, const std::vector<idx>& target,
///       const std::vector<idx>& dims)
TEST(qpp_apply, Qudits) {
    // pure states
    idx d = 3;
    // 1 qudit

    // 2 qudits
    ket psi = mket({0, 0}, d) + mket({1, 0}, d) / std::sqrt(2);
    ket result = apply(psi, gt.Xd(d), {0}, d);
    ket expected = mket({1, 0}, d) + mket({2, 0}, d) / std::sqrt(2);
    EXPECT_NEAR(0, norm(result - expected), 1e-5);

    // 4 qudits

    // mixed states
    // 2 qudits
    cmat rho = 0.6 * prj(mket({0, 0}, d) + mket({1, 0}, d) / std::sqrt(2)) +
               0.4 * prj(mket({0, 1}, d) + mket({1, 1}, d) / std::sqrt(2));
    cmat result_rho = apply(rho, gt.Xd(d), {0}, d);
    cmat expected_rho =
        0.6 * prj(mket({1, 0}, d) + mket({2, 0}, d) / std::sqrt(2)) +
        0.4 * prj(mket({1, 1}, d) + mket({2, 1}, d) / std::sqrt(2));
    EXPECT_NEAR(0, norm(result_rho - expected_rho), 1e-5);
}

/// BEGIN template <typename Derived1, typename Derived2> expr_t<Derived1>
///       apply(const Eigen::MatrixBase<Derived1>& state,
///       const Eigen::MatrixBase<Derived2>& A, const std::vector<idx>& target,
///       idx d = 2)
TEST(qpp_apply, Qubits) {
    // pure states
    // 1 qubit
    ket psi = 1_ket;
    // X, Y, Z and H
    ket resultX = apply(psi, gt.X, {0}, std::vector<idx>({2}));
    EXPECT_EQ(0_ket, resultX);
    ket resultY = apply(psi, gt.Y, {0}, std::vector<idx>({2}));
    EXPECT_EQ(-1_i * 0_ket, resultY);
    ket resultZ = apply(psi, gt.Z, {0}, std::vector<idx>({2}));
    EXPECT_EQ(-1_ket, resultZ);
    ket resultH = apply(psi, gt.H, {0}, std::vector<idx>({2}));
    EXPECT_NEAR(0, norm(resultH - (0_ket - 1_ket) / std::sqrt(2)), 1e-5);

    // 2 qubits
    psi = 0.8 * 00_ket + 0.6 * 11_ket;
    resultX = apply(psi, gt.X, {1}, {2, 2});
    EXPECT_EQ(0.8 * 01_ket + 0.6 * 10_ket, resultX);
    resultY = apply(psi, gt.Y, {1}, {2, 2});
    EXPECT_EQ(0.8_i * 01_ket - 0.6_i * 10_ket, resultY);
    resultZ = apply(psi, gt.Z, {1}, {2, 2});
    EXPECT_EQ(0.8 * 00_ket - 0.6 * 11_ket, resultZ);
    ket resultCNOT = apply(psi, gt.CNOT, {0, 1}, {2, 2});
    EXPECT_EQ(0.8 * 00_ket + 0.6 * 10_ket, resultCNOT);
    resultCNOT = apply(psi, gt.CNOT, {1, 0}, {2, 2});
    EXPECT_EQ(0.8 * 00_ket + 0.6 * 01_ket, resultCNOT);
    ket resultZZ = apply(psi, kron(gt.Z, gt.Z), {0, 1}, {2, 2});
    EXPECT_EQ(0.8 * 00_ket + 0.6 * 11_ket, resultZZ);

    // 4 qubits
    psi = 0.8 * 0000_ket + 0.6 * 1111_ket;
    ket resultXZ = apply(psi, kron(gt.X, gt.Z), {1, 2}, {2, 2, 2, 2});
    EXPECT_EQ(0.8 * 0100_ket - 0.6 * 1011_ket, resultXZ);
    resultXZ = apply(psi, kron(gt.X, gt.Z), {2, 1}, {2, 2, 2, 2});
    EXPECT_EQ(0.8 * 0010_ket - 0.6 * 1101_ket, resultXZ);
    ket resultTOF = apply(psi, gt.TOF, {1, 2, 0}, {2, 2, 2, 2});
    EXPECT_EQ(0.8 * 0000_ket + 0.6 * 0111_ket, resultTOF);
}

/// BEGIN template <typename Derived> cmat apply(
///       const Eigen::MatrixBase<Derived>& A, const std::vector<cmat>& Ks)
TEST(qpp_apply, Channel) {
    const cmat rho = prj(st.x0);
    EXPECT_LT((qpp::apply(rho, std::vector<cmat>{gt.Id(2)}) - rho).norm(),
              1e-12);
    EXPECT_LT(
        (qpp::apply(prj(st.z0), std::vector<cmat>{gt.X}) - prj(st.z1)).norm(),
        1e-12);
}

/// BEGIN template <typename Derived> cmat apply(
///       const Eigen::MatrixBase<Derived>& A, const std::vector<cmat>& Ks,
///       const std::vector<idx>& target, const std::vector<idx>& dims)
TEST(qpp_apply, SubsysChannel) {
    const cmat rho = mprj({2, 1}, std::vector<idx>{3, 2});
    EXPECT_LT((apply(rho, std::vector<cmat>{gt.Xd(3)}, {0}, {3, 2}) -
               mprj({0, 1}, std::vector<idx>{3, 2}))
                  .norm(),
              1e-12);
}

/// BEGIN template <typename Derived> cmat apply(
///       const Eigen::MatrixBase<Derived>& A, const std::vector<cmat>& Ks,
///       const std::vector<idx>& target, idx d = 2)
TEST(qpp_apply, SubsysChannelQubits) {
    const cmat rho = mprj({1, 0});
    EXPECT_LT((apply(rho, std::vector<cmat>{gt.X}, {1}) - mprj({1, 1})).norm(),
              1e-12);
}

/// BEGIN template <typename Derived1, typename Derived2> expr_t<Derived1>
///       applyCTRL(const Eigen::MatrixBase<Derived1>& state,
///       const Eigen::MatrixBase<Derived2>& A, const std::vector<idx>& ctrl,
///       const std::vector<idx>& target, const std::vector<idx>& dims,
///       std::optional<std::vector<idx>> shift = std::nullopt)
TEST(qpp_applyCTRL, Qudits) {
    idx d = 3;
    idx n = 4;
    std::vector<idx> dims(n, d); // 4 qutrits
    idx D = prod(dims);          // total dimension

    std::vector<idx> ctrl{0, 2};   // control
    std::vector<idx> target{1, 3}; // target

    idx Dtarget = 1; // dimension of the target subsystems
    for (auto i : target) {
        Dtarget *= dims[i]; // compute it here
    }

    // ket
    ket psi = mket({2, 1, 2, 1}, d);
    psi = applyCTRL(psi, kron(gt.Xd(3), gt.Xd(3)), ctrl, target, dims);
    ket expected = mket({2, 0, 2, 0}, d);
    EXPECT_NEAR(0, norm(psi - expected), 1e-5);

    // ket, with shift
    std::vector<idx> shift{1, 2};
    psi = mket({2, 1, 1, 1}, d); // will behave as if |0101>
    psi = applyCTRL(psi, kron(gt.Xd(3), gt.Xd(3)), ctrl, target, dims, shift);
    expected = mket({2, 1, 1, 1}, d);
    EXPECT_NEAR(0, norm(psi - expected), 1e-5);

    // density matrix
    psi = mket({2, 0, 2, 0}, d);
    cmat rho =
        applyCTRL(prj(psi), kron(gt.Xd(3), gt.Xd(3)), ctrl, target, dims);
    cmat expected_rho = mprj({2, 2, 2, 2}, d);
    EXPECT_NEAR(0, norm(rho - expected_rho), 1e-5);

    // density matrix, with shift
    shift = {1, 2};
    psi = mket({1, 0, 0, 0}, d); // will behave as if |2020>
    rho = applyCTRL(prj(psi), kron(gt.Xd(3), gt.Xd(3)), ctrl, target, dims,
                    shift);
    expected_rho = mprj({1, 2, 0, 2}, d);
    EXPECT_NEAR(0, norm(rho - expected_rho), 1e-5);

    // some random n qudit pure state
    psi = randket(D);

    rho = psi * adjoint(psi); // the corresponding density matrix
    cmat U = randU(Dtarget);  // some random unitary on the target

    // applyCTRL on pure state
    ket A = applyCTRL(psi, U, ctrl, target, dims);

    // applyCTRL on density matrix
    cmat B = applyCTRL(rho, U, ctrl, target, dims);

    // result when using CTRL-U|psi><psi|CTRL-U^\dagger
    cmat result_psi = A * adjoint(A);
    // result when using CTRL-U(rho)CTRL-U^\dagger
    const cmat& result_rho = B;

    realT res = norm(result_psi - result_rho);
    EXPECT_NEAR(0, res, 1e-5);
}

/// BEGIN template <typename Derived1, typename Derived2> expr_t<Derived1>
///       applyCTRL(const Eigen::MatrixBase<Derived1>& state,
///       const Eigen::MatrixBase<Derived2>& A, const std::vector<idx>& ctrl,
///       const std::vector<idx>& target, idx d = 2,
///       std::optional<std::vector<idx>> shift = std::nullopt)
TEST(qpp_applyCTRL, Qubits) {
    EXPECT_LT((applyCTRL(mket({1, 0}), gt.X, {0}, {1}) - mket({1, 1})).norm(),
              1e-12);
    EXPECT_LT((applyCTRL(mket({0, 0}), gt.X, {0}, {1}) - mket({0, 0})).norm(),
              1e-12);
    EXPECT_LT((applyCTRL(mprj({1, 0}), gt.X, {0}, {1}) - mprj({1, 1})).norm(),
              1e-12);
}

/// BEGIN template <typename Derived1, typename Derived2> expr_t<Derived1>
///       applyCTRL_fan(const Eigen::MatrixBase<Derived1>& state,
///       const Eigen::MatrixBase<Derived2>& A, const std::vector<idx>& ctrl,
///       const std::vector<idx>& target, const std::vector<idx>& dims,
///       std::optional<std::vector<idx>> shift = std::nullopt)
TEST(qpp_applyCTRL_fan, Qudits) {
    idx d = 3;
    idx n = 4;
    std::vector<idx> dims(n, d); // 4 qutrits
    idx D = prod(dims);          // total dimension

    std::vector<idx> ctrl{0, 2};   // control
    std::vector<idx> target{1, 3}; // target

    // ket
    ket psi = mket({2, 1, 2, 1}, d);
    psi = applyCTRL_fan(psi, gt.Xd(3), ctrl, target, dims);
    ket expected = mket({2, 0, 2, 0}, d);
    EXPECT_NEAR(0, norm(psi - expected), 1e-5);

    // ket, with shift
    std::vector<idx> shift{1, 2};
    psi = mket({2, 1, 1, 1}, d); // will behave as if |0101>
    psi = applyCTRL_fan(psi, gt.Xd(3), ctrl, target, dims, shift);
    expected = mket({2, 1, 1, 1}, d);
    EXPECT_NEAR(0, norm(psi - expected), 1e-5);

    // density matrix
    psi = mket({2, 0, 2, 0}, d);
    cmat rho = applyCTRL_fan(prj(psi), gt.Xd(3), ctrl, target, dims);
    cmat expected_rho = mprj({2, 2, 2, 2}, d);
    EXPECT_NEAR(0, norm(rho - expected_rho), 1e-5);

    // density matrix, with shift
    shift = {1, 2};
    psi = mket({1, 0, 0, 0}, d); // will behave as if |2020>
    rho = applyCTRL_fan(prj(psi), gt.Xd(3), ctrl, target, dims, shift);
    expected_rho = mprj({1, 2, 0, 2}, d);
    EXPECT_NEAR(0, norm(rho - expected_rho), 1e-5);

    // some random n qudit pure state
    psi = randket(D);

    rho = psi * adjoint(psi); // the corresponding density matrix
    cmat U = randU(d);        // some random unitary

    // applyCTRL on pure state
    ket A = applyCTRL_fan(psi, U, ctrl, target, dims);

    // applyCTRL on density matrix
    cmat B = applyCTRL_fan(rho, U, ctrl, target, dims);

    // result when using CTRL-U|psi><psi|CTRL-U^\dagger
    cmat result_psi = A * adjoint(A);
    // result when using CTRL-U(rho)CTRL-U^\dagger
    const cmat& result_rho = B;

    realT res = norm(result_psi - result_rho);
    EXPECT_NEAR(0, res, 1e-5);
}

/// BEGIN template <typename Derived1, typename Derived2> expr_t<Derived1>
///       applyCTRL_fan(const Eigen::MatrixBase<Derived1>& state,
///       const Eigen::MatrixBase<Derived2>& A, const std::vector<idx>& ctrl,
///       const std::vector<idx>& target, idx d = 2,
///       std::optional<std::vector<idx>> shift = std::nullopt)
TEST(qpp_applyCTRL_fan, Qubits) {
    const idx n = 6; // total qubits
    const idx D = internal::ipow_rounded<idx>(2, n);

    // worst-ish case: half controls, half targets (disjoint)
    const idx m = n / 2;
    std::vector<idx> ctrl(m);
    std::vector<idx> target(m);
    for (idx i = 0; i < m; ++i) {
        ctrl[i] = i;       // 0..m-1
        target[i] = m + i; // m..2m-1
    }

    // -------------------------------------------------------------------------
    // A) ket: deterministic, no shift
    // -------------------------------------------------------------------------
    {
        // controls all |1>, targets all |0> -> X flips all targets -> targets
        // |1>
        std::vector<idx> midx(n, 0);
        for (idx i = 0; i < m; ++i) {
            midx[ctrl[i]] = 1;
        }
        ket psi = mket(midx, 2);

        // signature expects shift, so pass the "no shift" vector
        const std::vector<idx> shift(m, 0);
        ket out = qpp::internal::kernels::qubit::apply_ctrl_fan_psi(
            psi, gt.X, ctrl, target, shift, n);

        std::vector<idx> exp_midx = midx;
        for (idx i = 0; i < m; ++i) {
            exp_midx[target[i]] = 1;
        }
        ket expected = mket(exp_midx, 2);

        EXPECT_NEAR(0.0, norm(out - expected), 1e-12);
    }

    // -------------------------------------------------------------------------
    // B) ket: deterministic, with shift
    // -------------------------------------------------------------------------
    {
        // negate first and last control
        std::vector<idx> shift(m, 0);
        shift[0] = 1;
        shift[m - 1] = 1;

        // ctrl pattern required: (0,1,0) on ctrl[0..2]
        std::vector<idx> midx(n, 0);
        midx[ctrl[0]] = 0;
        midx[ctrl[1]] = 1;
        midx[ctrl[2]] = 0;
        ket psi = mket(midx, 2);

        ket out = qpp::internal::kernels::qubit::apply_ctrl_fan_psi(
            psi, gt.X, ctrl, target, shift, n);

        std::vector<idx> exp_midx = midx;
        for (idx i = 0; i < m; ++i) {
            exp_midx[target[i]] = 1;
        }
        ket expected = mket(exp_midx, 2);

        EXPECT_NEAR(0.0, norm(out - expected), 1e-12);
    }

    // -------------------------------------------------------------------------
    // C) ket: control pattern doesn't match -> unchanged
    // -------------------------------------------------------------------------
    {
        std::vector<idx> shift(m, 0);
        shift[0] = 1;
        shift[m - 1] = 1;

        // start from the matching pattern, then break ctrl[1]
        std::vector<idx> midx(n, 0);
        midx[ctrl[0]] = 0;
        midx[ctrl[1]] = 1;
        midx[ctrl[2]] = 0;
        midx[ctrl[1]] = 0; // break (non-negated) control condition

        ket psi = mket(midx, 2);
        ket out = qpp::internal::kernels::qubit::apply_ctrl_fan_psi(
            psi, gt.X, ctrl, target, shift, n);

        EXPECT_NEAR(0.0, norm(out - psi), 1e-12);
    }

    // -------------------------------------------------------------------------
    // D) density matrix: deterministic, no shift
    // -------------------------------------------------------------------------
    {
        std::vector<idx> midx(n, 0);
        for (idx i = 0; i < m; ++i) {
            midx[ctrl[i]] = 1;
        }
        ket psi = mket(midx, 2);

        const std::vector<idx> shift(m, 0);
        cmat rho = qpp::internal::kernels::qubit::apply_ctrl_fan_rho(
            prj(psi), gt.X, ctrl, target, shift, n);

        std::vector<idx> exp_midx = midx;
        for (idx i = 0; i < m; ++i) {
            exp_midx[target[i]] = 1;
        }
        cmat expected_rho = prj(mket(exp_midx, 2));

        EXPECT_NEAR(0.0, norm(rho - expected_rho), 1e-12);
    }

    // -------------------------------------------------------------------------
    // E) density matrix: deterministic, with shift
    // -------------------------------------------------------------------------
    {
        std::vector<idx> shift(m, 0);
        shift[0] = 1;
        shift[m - 1] = 1;

        std::vector<idx> midx(n, 0);
        midx[ctrl[0]] = 0;
        midx[ctrl[1]] = 1;
        midx[ctrl[2]] = 0;
        ket psi = mket(midx, 2);

        cmat rho = qpp::internal::kernels::qubit::apply_ctrl_fan_rho(
            prj(psi), gt.X, ctrl, target, shift, n);

        std::vector<idx> exp_midx = midx;
        for (idx i = 0; i < m; ++i) {
            exp_midx[target[i]] = 1;
        }
        cmat expected_rho = prj(mket(exp_midx, 2));

        EXPECT_NEAR(0.0, norm(rho - expected_rho), 1e-12);
    }

    // -------------------------------------------------------------------------
    // F) random consistency: ket vs density (same operation)
    // -------------------------------------------------------------------------
    {
        ket psi = randket(D);
        cmat rho = psi * adjoint(psi);

        const std::vector<idx> shift(m, 0);
        const cmat U = randU(2);

        ket A = qpp::internal::kernels::qubit::apply_ctrl_fan_psi(
            psi, U, ctrl, target, shift, n);
        cmat B = qpp::internal::kernels::qubit::apply_ctrl_fan_rho(
            rho, U, ctrl, target, shift, n);

        const cmat result_psi = A * adjoint(A);
        const realT res = norm(result_psi - B);

        EXPECT_NEAR(0.0, norm(result_psi - adjoint(result_psi)), 1e-10);
        EXPECT_NEAR(0.0, norm(B - adjoint(B)), 1e-10);
        EXPECT_NEAR(0.0, res, 1e-10);
    }
}

/// BEGIN template <typename Derived> expr_t<Derived>
///       applyQFT(const Eigen::MatrixBase<Derived>& A,
///       const std::vector<idx>& target, idx d = 2, bool swap = true)
TEST(qpp_applyQFT, AllTests) {
    const ket psi = mket({0, 0});
    const ket expected = (mket({0, 0}) + mket({0, 1})) / std::sqrt(2.);
    EXPECT_LT((applyQFT(psi, {1}) - expected).norm(), 1e-12);
    EXPECT_LT((applyQFT(prj(psi), {1}) - prj(expected)).norm(), 1e-12);
    EXPECT_LT((applyTFQ(applyQFT(psi, {1}), {1}) - psi).norm(), 1e-10);
}

/// BEGIN template <typename Derived> expr_t<Derived>
///       applyTFQ(const Eigen::MatrixBase<Derived>& A,
///       const std::vector<idx>& target, idx d = 2, bool swap = true)
TEST(qpp_applyTQF, AllTests) {
    const ket psi = mket({0, 1});
    const ket zero = mket({0, 0});
    const ket expected = (mket({0, 0}) + mket({0, 1})) / std::sqrt(2.);
    EXPECT_LT((applyTFQ(zero, {1}) - expected).norm(), 1e-12);
    EXPECT_LT((applyQFT(applyTFQ(psi, {0, 1}), {0, 1}) - psi).norm(), 1e-10);
    EXPECT_LT((applyTFQ(applyQFT(prj(psi), {0, 1}), {0, 1}) - prj(psi)).norm(),
              1e-10);
}

/// BEGIN inline std::vector<cmat> choi2kraus(const cmat& A, idx Din, idx Dout)
TEST(qpp_choi2kraus, AllTests) {
    const std::vector<cmat> Ks{gt.X};
    const cmat choi = kraus2choi(Ks);
    EXPECT_LT((kraus2super(choi2kraus(choi, 2, 2)) - kraus2super(Ks)).norm(),
              1e-10);
}

/// BEGIN inline cmat choi2super(const cmat& A, idx Din, idx Dout)
TEST(qpp_choi2super, AllTests) {
    const cmat choi = kraus2choi(std::vector<cmat>{gt.X});
    EXPECT_LT(
        (choi2super(choi, 2, 2) - kraus2super(std::vector<cmat>{gt.X})).norm(),
        1e-12);
}

/// BEGIN inline cmat kraus2choi(const std::vector<cmat>& Ks)
TEST(qpp_kraus2choi, AllTests) {
    const cmat choi = kraus2choi(std::vector<cmat>{gt.Id(2)});
    const ket omega = mket({0, 0}) + mket({1, 1});
    EXPECT_LT((choi - omega * adjoint(omega)).norm(), 1e-12);
}

/// BEGIN inline cmat kraus2super(const std::vector<cmat>& Ks)
TEST(qpp_kraus2super, AllTests) {
    EXPECT_LT((kraus2super(std::vector<cmat>{gt.Id(2)}) - cmat::Identity(4, 4))
                  .norm(),
              1e-12);
}

/// BEGIN template <typename Derived> dyn_mat<typename Derived::Scalar>
///       ptrace(const Eigen::MatrixBase<Derived>& A,
///       const std::vector<idx>& target, const std::vector<idx>& dims)
TEST(qpp_ptrace, Qudits) {
    const ket psi = mket({2, 1}, std::vector<idx>{3, 2});
    EXPECT_LT((ptrace(psi, {0}, {3, 2}) - mprj({1})).norm(), 1e-12);
    EXPECT_LT((ptrace(prj(psi), {1}, {3, 2}) - mprj({2}, 3)).norm(), 1e-12);
}

/// BEGIN template <typename Derived> dyn_mat<typename Derived::Scalar>
///       ptrace(const Eigen::MatrixBase<Derived>& A,
///       const std::vector<idx>& target, idx d = 2)
TEST(qpp_ptrace, Qubits) {
    EXPECT_LT((ptrace(st.b00, {0}) - cmat::Identity(2, 2) / 2.).norm(), 1e-12);
    EXPECT_LT((ptrace(mprj({1, 0}), {1}) - mprj({1})).norm(), 1e-12);
}

/// BEGIN template <typename Derived> dyn_mat<typename Derived::Scalar>
///       ptrace1(const Eigen::MatrixBase<Derived>& A,
///       const std::vector<idx>& dims)
TEST(qpp_ptrace1, Qudits) {
    const ket psi = mket({2, 1}, std::vector<idx>{3, 2});
    EXPECT_LT((ptrace1(psi, {3, 2}) - mprj({1})).norm(), 1e-12);
}

/// BEGIN template <typename Derived> dyn_mat<typename Derived::Scalar>
///       ptrace1(const Eigen::MatrixBase<Derived>& A, idx d = 2)
TEST(qpp_ptrace1, Qubits) {
    EXPECT_LT((ptrace1(mprj({1, 0})) - mprj({0})).norm(), 1e-12);
}

/// BEGIN template <typename Derived> dyn_mat<typename Derived::Scalar>
///       ptrace2(const Eigen::MatrixBase<Derived>& A,
///       const std::vector<idx>& dims)
TEST(qpp_ptrace2, Qudits) {
    const ket psi = mket({2, 1}, std::vector<idx>{3, 2});
    EXPECT_LT((ptrace2(prj(psi), {3, 2}) - mprj({2}, 3)).norm(), 1e-12);
}

/// BEGIN template <typename Derived> dyn_mat<typename Derived::Scalar>
///       ptrace2(const Eigen::MatrixBase<Derived>& A, idx d = 2)
TEST(qpp_ptrace2, Qubits) {
    EXPECT_LT((ptrace2(mprj({1, 0})) - mprj({1})).norm(), 1e-12);
}

/// BEGIN template <typename Derived> dyn_mat<typename Derived::Scalar>
///       ptranspose(const Eigen::MatrixBase<Derived>& A,
///       const std::vector<idx>& target, const std::vector<idx>& dims)
TEST(qpp_ptranspose, Qudits) {
    const cmat rho = prj(st.mes(3));
    const cmat pt = ptranspose(rho, {0}, {3, 3});
    EXPECT_LT((ptranspose(pt, {0}, {3, 3}) - rho).norm(), 1e-12);
    EXPECT_NEAR(1, std::real(trace(pt)), 1e-12);
}

/// BEGIN template <typename Derived> dyn_mat<typename Derived::Scalar>
///       ptranspose(const Eigen::MatrixBase<Derived>& A,
///       const std::vector<idx>& target, idx d = 2)
TEST(qpp_ptranspose, Qubits) {
    const cmat rho = prj(st.b00);
    const cmat pt = ptranspose(rho, {1});
    EXPECT_LT((ptranspose(pt, {1}) - rho).norm(), 1e-12);
    EXPECT_GT((pt - rho).norm(), 0.1);
}

/// BEGIN template <typename Derived> expr_t<Derived>
///       QFT(const Eigen::MatrixBase<Derived>& A, idx d = 2, bool swap = true)
TEST(qpp_QFT, AllTests) {
    const ket psi = mket({0, 0});
    const ket expected = ket::Ones(4) / 2.;
    EXPECT_LT((QFT(psi) - expected).norm(), 1e-12);
    EXPECT_LT((QFT(prj(psi)) - prj(expected)).norm(), 1e-12);
    EXPECT_LT((TFQ(QFT(psi)) - psi).norm(), 1e-10);
}

/// BEGIN template <typename Derived> dyn_col_vect<typename Derived::Scalar>
///       qRAM(const Eigen::MatrixBase<Derived>& psi, const qram& data,
///       idx DqRAM)
TEST(qpp_qRAM, AllTests) {
    const ket psi = st.x0;
    const qram data{1, 0};
    const ket expected = (mket({0, 1}) + mket({1, 0})) / std::sqrt(2.);
    EXPECT_LT((qRAM(psi, data, 2) - expected).norm(), 1e-12);
    EXPECT_THROW(qRAM(psi, data, 1), exception::OutOfRange);
}

/// BEGIN inline cmat super2choi(const cmat& A)
TEST(qpp_super2choi, AllTests) {
    const cmat super = kraus2super(std::vector<cmat>{gt.X});
    EXPECT_LT((super2choi(super) - kraus2choi(std::vector<cmat>{gt.X})).norm(),
              1e-12);
}

/// BEGIN inline cmat super2kraus(const cmat& A)
TEST(qpp_super2kraus, AllTests) {
    const cmat super = kraus2super(std::vector<cmat>{gt.X});
    EXPECT_LT((kraus2super(super2kraus(super)) - super).norm(), 1e-10);
}

/// BEGIN template <typename Derived> dyn_mat<typename Derived::Scalar>
///       syspermute(const Eigen::MatrixBase<Derived>& A,
///       const std::vector<idx>& perm, const std::vector<idx>& dims)
TEST(qpp_syspermute, Qudits) {
    const ket psi = mket({2, 1}, std::vector<idx>{3, 2});
    EXPECT_LT(
        (syspermute(psi, {1, 0}, {3, 2}) - mket({1, 2}, std::vector<idx>{2, 3}))
            .norm(),
        1e-12);
    EXPECT_LT((syspermute(prj(psi), {1, 0}, {3, 2}) -
               mprj({1, 2}, std::vector<idx>{2, 3}))
                  .norm(),
              1e-12);
}

/// BEGIN template <typename Derived> dyn_mat<typename Derived::Scalar>
///       syspermute(const Eigen::MatrixBase<Derived>& A,
///       const std::vector<idx>& perm, idx d = 2)
TEST(qpp_syspermute, Qubits) {
    EXPECT_LT((syspermute(mket({1, 0, 1}), {2, 0, 1}) - mket({1, 1, 0})).norm(),
              1e-12);
}

/// BEGIN template <typename Derived> expr_t<Derived>
///       TFQ(const Eigen::MatrixBase<Derived>& A, idx d = 2, bool swap = true)
TEST(qpp_TFQ, AllTests) {
    const ket psi = mket({1, 0});
    const ket zero = mket({0, 0});
    const ket expected = ket::Ones(4) / 2.;
    EXPECT_LT((TFQ(zero) - expected).norm(), 1e-12);
    EXPECT_LT((QFT(TFQ(psi)) - psi).norm(), 1e-10);
    EXPECT_LT((TFQ(QFT(prj(psi))) - prj(psi)).norm(), 1e-10);
}
