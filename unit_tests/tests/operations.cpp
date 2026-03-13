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

    idx d = 3;
    // 1 qudit

    // 2 qudits
    psi = mket({0, 0}, d) + mket({1, 0}, d) / std::sqrt(2);
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
TEST(qpp_apply, Qubits) {}

/// BEGIN template <typename Derived> cmat apply(
///       const Eigen::MatrixBase<Derived>& A, const std::vector<cmat>& Ks)
TEST(qpp_apply, Channel) {}

/// BEGIN template <typename Derived> cmat apply(
///       const Eigen::MatrixBase<Derived>& A, const std::vector<cmat>& Ks,
///       const std::vector<idx>& target, const std::vector<idx>& dims)
TEST(qpp_apply, SubsysChannel) {}

/// BEGIN template <typename Derived> cmat apply(
///       const Eigen::MatrixBase<Derived>& A, const std::vector<cmat>& Ks,
///       const std::vector<idx>& target, idx d = 2)
TEST(qpp_apply, SubsysChannelQubits) {}

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
TEST(qpp_applyCTRL, Qubits) {}

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
    const idx D = internal::safe_pow<idx>(2, n);

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
TEST(qpp_applyQFT, AllTests) {}

/// BEGIN template <typename Derived> expr_t<Derived>
///       applyTFQ(const Eigen::MatrixBase<Derived>& A,
///       const std::vector<idx>& target, idx d = 2, bool swap = true)
TEST(qpp_applyTQF, AllTests) {}

/// BEGIN inline std::vector<cmat> choi2kraus(const cmat& A, idx Din, idx Dout)
TEST(qpp_choi2kraus, AllTests) {}

/// BEGIN inline cmat choi2super(const cmat& A, idx Din, idx Dout)
TEST(qpp_choi2super, AllTests) {}

/// BEGIN inline cmat kraus2choi(const std::vector<cmat>& Ks)
TEST(qpp_kraus2choi, AllTests) {}

/// BEGIN inline cmat kraus2super(const std::vector<cmat>& Ks)
TEST(qpp_kraus2super, AllTests) {}

/// BEGIN template <typename Derived> dyn_mat<typename Derived::Scalar>
///       ptrace(const Eigen::MatrixBase<Derived>& A,
///       const std::vector<idx>& target, const std::vector<idx>& dims)
TEST(qpp_ptrace, Qudits) {}

/// BEGIN template <typename Derived> dyn_mat<typename Derived::Scalar>
///       ptrace(const Eigen::MatrixBase<Derived>& A,
///       const std::vector<idx>& target, idx d = 2)
TEST(qpp_ptrace, Qubits) {}

/// BEGIN template <typename Derived> dyn_mat<typename Derived::Scalar>
///       ptrace1(const Eigen::MatrixBase<Derived>& A,
///       const std::vector<idx>& dims)
TEST(qpp_ptrace1, Qudits) {}

/// BEGIN template <typename Derived> dyn_mat<typename Derived::Scalar>
///       ptrace1(const Eigen::MatrixBase<Derived>& A, idx d = 2)
TEST(qpp_ptrace1, Qubits) {}

/// BEGIN template <typename Derived> dyn_mat<typename Derived::Scalar>
///       ptrace2(const Eigen::MatrixBase<Derived>& A,
///       const std::vector<idx>& dims)
TEST(qpp_ptrace2, Qudits) {}

/// BEGIN template <typename Derived> dyn_mat<typename Derived::Scalar>
///       ptrace2(const Eigen::MatrixBase<Derived>& A, idx d = 2)
TEST(qpp_ptrace2, Qubits) {}

/// BEGIN template <typename Derived> dyn_mat<typename Derived::Scalar>
///       ptranspose(const Eigen::MatrixBase<Derived>& A,
///       const std::vector<idx>& target, const std::vector<idx>& dims)
TEST(qpp_ptranspose, Qudits) {}

/// BEGIN template <typename Derived> dyn_mat<typename Derived::Scalar>
///       ptranspose(const Eigen::MatrixBase<Derived>& A,
///       const std::vector<idx>& target, idx d = 2)
TEST(qpp_ptranspose, Qubits) {}

/// BEGIN template <typename Derived> expr_t<Derived>
///       QFT(const Eigen::MatrixBase<Derived>& A, idx d = 2, bool swap = true)
TEST(qpp_QFT, AllTests) {}

/// BEGIN template <typename Derived> dyn_col_vect<typename Derived::Scalar>
///       qRAM(const Eigen::MatrixBase<Derived>& psi, const qram& data,
///       idx DqRAM)
TEST(qpp_qRAM, AllTests) {}

/// BEGIN inline cmat super2choi(const cmat& A)
TEST(qpp_super2choi, AllTests) {}

/// BEGIN inline cmat super2kraus(const cmat& A)
TEST(qpp_super2kraus, AllTests) {}

/// BEGIN template <typename Derived> dyn_mat<typename Derived::Scalar>
///       syspermute(const Eigen::MatrixBase<Derived>& A,
///       const std::vector<idx>& perm, const std::vector<idx>& dims)
TEST(qpp_syspermute, Qudits) {}

/// BEGIN template <typename Derived> dyn_mat<typename Derived::Scalar>
///       syspermute(const Eigen::MatrixBase<Derived>& A,
///       const std::vector<idx>& perm, idx d = 2)
TEST(qpp_syspermute, Qubits) {}

/// BEGIN template <typename Derived> expr_t<Derived>
///       TFQ(const Eigen::MatrixBase<Derived>& A, idx d = 2, bool swap = true)
TEST(qpp_TFQ, AllTests) {}
