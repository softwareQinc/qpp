#include <cmath>
#include <complex>

#include "gtest/gtest.h"

#include "qpp.h"

using namespace qpp;

// Unit testing "classes/gates.hpp"

/******************************************************************************/
TEST(qpp_Gates_CTRL, Qudits) {
    idx d = 3; // qutrits

    // CNOT control-target on 2 qutrits
    cmat CTRL1 = gt.CTRL(gt.Xd(d), {0}, {1}, 2, d);
    EXPECT_NEAR(0, norm(CTRL1 * mket({0, 0}, {d, d}) - mket({0, 0}, {d, d})),
                1e-5);
    EXPECT_NEAR(0, norm(CTRL1 * mket({0, 1}, {d, d}) - mket({0, 1}, {d, d})),
                1e-5);
    EXPECT_NEAR(0, norm(CTRL1 * mket({0, 2}, {d, d}) - mket({0, 2}, {d, d})),
                1e-5);
    EXPECT_NEAR(0, norm(CTRL1 * mket({1, 0}, {d, d}) - mket({1, 1}, {d, d})),
                1e-5);
    EXPECT_NEAR(0, norm(CTRL1 * mket({1, 1}, {d, d}) - mket({1, 2}, {d, d})),
                1e-5);
    EXPECT_NEAR(0, norm(CTRL1 * mket({1, 2}, {d, d}) - mket({1, 0}, {d, d})),
                1e-5);
    EXPECT_NEAR(0, norm(CTRL1 * mket({2, 0}, {d, d}) - mket({2, 2}, {d, d})),
                1e-5);
    EXPECT_NEAR(0, norm(CTRL1 * mket({2, 1}, {d, d}) - mket({2, 0}, {d, d})),
                1e-5);
    EXPECT_NEAR(0, norm(CTRL1 * mket({2, 2}, {d, d}) - mket({2, 1}, {d, d})),
                1e-5);

    // CNOT target-control on 2 qutrits
    cmat CTRL2 = gt.CTRL(gt.Xd(d), {1}, {0}, 2, d);
    EXPECT_NEAR(0, norm(CTRL2 * mket({0, 0}, {d, d}) - mket({0, 0}, {d, d})),
                1e-5);
    EXPECT_NEAR(0, norm(CTRL2 * mket({0, 1}, {d, d}) - mket({1, 1}, {d, d})),
                1e-5);
    EXPECT_NEAR(0, norm(CTRL2 * mket({0, 2}, {d, d}) - mket({2, 2}, {d, d})),
                1e-5);
    EXPECT_NEAR(0, norm(CTRL2 * mket({1, 0}, {d, d}) - mket({1, 0}, {d, d})),
                1e-5);
    EXPECT_NEAR(0, norm(CTRL2 * mket({1, 1}, {d, d}) - mket({2, 1}, {d, d})),
                1e-5);
    EXPECT_NEAR(0, norm(CTRL2 * mket({1, 2}, {d, d}) - mket({0, 2}, {d, d})),
                1e-5);
    EXPECT_NEAR(0, norm(CTRL2 * mket({2, 0}, {d, d}) - mket({2, 0}, {d, d})),
                1e-5);
    EXPECT_NEAR(0, norm(CTRL2 * mket({2, 1}, {d, d}) - mket({0, 1}, {d, d})),
                1e-5);
    EXPECT_NEAR(0, norm(CTRL2 * mket({2, 2}, {d, d}) - mket({1, 2}, {d, d})),
                1e-5);

    // multiple Control-X-X, partial testing
    cmat CTRL3 = gt.CTRL(kron(gt.Xd(d), gt.Xd(3)), {1, 4}, {2, 3}, 6, d);
    ket psi1 = mket({0, 1, 2, 2, 1, 1}, {d, d, d, d, d, d});
    ket res1 = mket({0, 1, 0, 0, 1, 1}, {d, d, d, d, d, d});
    EXPECT_NEAR(0, norm(CTRL3 * psi1 - res1), 1e-5);

    ket psi2 = mket({0, 1, 2, 2, 2, 1}, {d, d, d, d, d, d});
    ket res2 = mket({0, 1, 2, 2, 2, 1}, {d, d, d, d, d, d});
    EXPECT_NEAR(0, norm(CTRL3 * psi2 - res2), 1e-5);

    ket psi3 = mket({1, 2, 1, 0, 2, 2}, {d, d, d, d, d, d});
    ket res3 = mket({1, 2, 0, 2, 2, 2}, {d, d, d, d, d, d});
    EXPECT_NEAR(0, norm(CTRL3 * psi3 - res3), 1e-5);
}
/******************************************************************************/
/// BEGIN template <typename Derived> dyn_mat<typename Derived::Scalar>
///       Gates::CTRL(const Eigen::MatrixBase<Derived>& A,
///       const std::vector<idx>& ctrl, const std::vector<idx>& target, idx N,
///       idx d = 2) const
/******************************************************************************/
TEST(qpp_Gates_CTRL, Qubits) {
    // CNOT control-target on 2 qubits
    cmat CTRL1 = gt.CTRL(gt.X, {0}, {1}, 2);
    EXPECT_EQ(CTRL1, gt.CNOT);

    // CNOT target-control on 2 qubits
    cmat CTRL2 = gt.CTRL(gt.X, {1}, {0}, 2);
    EXPECT_EQ(CTRL2, gt.CNOTba);

    // TOFFOLI
    cmat CTRL3 = gt.CTRL(gt.X, {0, 1}, {2}, 3);
    EXPECT_EQ(CTRL3, gt.TOF);
    CTRL3 = gt.CTRL(gt.X, {0, 1}, {2}, 3, 2); // test non-default args
    EXPECT_EQ(CTRL3, gt.TOF);

    // random gate as multiple control on 2 qubits
    // cmat U = randU();
    cmat U = gt.X;
    cmat CTRL4 = gt.CTRL(U, {0, 2}, {1}, 3);
    ket psi1 = mket({0, 0, 1});
    ket res1 = mket({0, 0, 1});
    EXPECT_NEAR(0, norm(CTRL4 * psi1 - res1), 1e-5);

    ket psi2 = mket({1, 1, 1});
    ket res2 = kron(st.z1, U * st.z1, st.z1);
    EXPECT_NEAR(0, norm(CTRL4 * psi2 - res2), 1e-5);
}
/******************************************************************************/
/// BEGIN template <typename Derived>
///       dyn_mat<typename Derived::Scalar> Gates::expandout(
///       const Eigen::MatrixBase<Derived>& A, idx pos,
///       const std::initializer_list<idx>& dims) const
TEST(qpp_Gates_expandout, InitList) {
    // single qutrit (degenerate case) random gate expansion
    cmat U = randU(3);
    EXPECT_EQ(gt.expandout(U, 0, {3}), U);
}
/******************************************************************************/
/// BEGIN template <typename Derived>
///       dyn_mat<typename Derived::Scalar> Gates::expandout(
///       const Eigen::MatrixBase<Derived>& A, idx pos,
///       const std::vector<idx>& dims) const
TEST(qpp_Gates_expandout, Qudits) {
    // single qubit (degenerate case) random gate expansion
    cmat U = randU();
    EXPECT_EQ(gt.expandout(U, 0, std::vector<idx>{2}), U);

    // 4 qutrits, identity on qutrit 3 expansion
    EXPECT_EQ(gt.expandout(gt.Id(3), 2, {3, 3, 3, 3}), gt.Id(81));

    // 3 qubits, X on qubit 2 expansion
    EXPECT_EQ(gt.expandout(gt.X, 1, {2, 2, 2}), kron(gt.Id2, gt.X, gt.Id2));
}
/******************************************************************************/
/// BEGIN template <typename Derived> dyn_mat<typename Derived::Scalar>
///       expandout(const Eigen::MatrixBase<Derived>& A, idx pos, idx N,
///       idx d = 2) const
TEST(qpp_Gates_expandout, Qubits) {
    // single qubit (degenerate case) random gate expansion
    cmat U = randU();
    EXPECT_EQ(gt.expandout(U, 0, 1), U);

    // 4 qutrits, identity on qutrit 3 expansion
    EXPECT_EQ(gt.expandout(gt.Id(3), 2, 4, 3), gt.Id(81));

    // 3 qubits, X on qubit 2 expansion
    EXPECT_EQ(gt.expandout(gt.X, 1, 3), kron(gt.Id2, gt.X, gt.Id2));
}
/******************************************************************************/
/// BEGIN cmat Gates::Fd(idx D = 2) const
TEST(qpp_Gates_Fd, AllTests) {
    EXPECT_NEAR(0, norm(gt.Fd(1) - gt.Id(1)), 1e-5);

    EXPECT_NEAR(0, norm(gt.Fd(2) - gt.H), 1e-5);

    cmat F3(3, 3);
    cplx o3 = omega(3);
    F3 << 1, 1, 1, 1, o3, o3 * o3, 1, o3 * o3, o3;
    F3 /= std::sqrt(3);
    EXPECT_NEAR(0, norm(gt.Fd(3) - F3), 1e-5);

    cmat F4(4, 4);
    cplx o4 = omega(4);
    F4 << 1, 1, 1, 1, 1, o4, o4 * o4, o4 * o4 * o4, 1, o4 * o4, 1, o4 * o4, 1,
        o4 * o4 * o4, o4 * o4, o4;
    F4 /= std::sqrt(4);
    EXPECT_NEAR(0, norm(gt.Fd(4) - F4), 1e-5);
}
/******************************************************************************/
/// BEGIN template <typename Derived> dyn_mat<typename Derived::Scalar>
///       Gates::GATE(const Eigen::MatrixBase<Derived>& A,
///       const std::vector<idx>& target, idx n,
///       const std::vector<idx>& dims) const
TEST(qpp_Gates_GATE, Qudits) {
    idx d = 3;
    idx n = 4; // 4 qutrits
    std::vector<idx> dims(n, d);

    ket psi = mket({1, 1, 1, 1}, d);
    cmat U = kron(gt.Xd(d), gt.Xd(d) * gt.Xd(d));
    cmat res = gt.GATE(U, {1, 3}, dims);
    EXPECT_NEAR(0, norm(res * psi - mket({1, 2, 1, 0}, d)), 1e-5);
    res = gt.GATE(U, {3, 1}, dims);
    EXPECT_NEAR(0, norm(res * psi - mket({1, 0, 1, 2}, d)), 1e-5);

    // random matrix on 2 qutrits
    U = randU(d * d);
    res = gt.GATE(U, {1, 3}, dims);
    EXPECT_NEAR(0, norm(res * psi - apply(psi, U, {1, 3}, dims)), 1e-5);
    res = gt.GATE(U, {3, 1}, dims);
    EXPECT_NE(0, norm(res * psi - apply(psi, U, {1, 3}, dims)));
}
/******************************************************************************/
/// BEGIN template <typename Derived> dyn_mat<typename Derived::Scalar>
///       Gates::GATE(const Eigen::MatrixBase<Derived>& A,
///       const std::vector<idx>& target, idx n, idx d = 2) const
TEST(qpp_Gates_GATE, Qubits) {
    idx n = 4; // 4 qubits

    ket psi = mket({1, 1, 1, 1});
    cmat U = gt.CNOT;
    cmat res = gt.GATE(U, {1, 3}, 4);
    EXPECT_NEAR(0, norm(res * psi - mket({1, 1, 1, 0})), 1e-5);
    res = gt.GATE(U, {3, 1}, 4);
    EXPECT_NEAR(0, norm(res * psi - mket({1, 0, 1, 1})), 1e-5);

    // random matrix on 2 qubits
    U = randU(4);
    res = gt.GATE(U, {1, 3}, 4);
    EXPECT_NEAR(0, norm(res * psi - apply(psi, U, {1, 3})), 1e-5);
    res = gt.GATE(U, {3, 1}, 4);
    EXPECT_NE(0, norm(res * psi - apply(psi, U, {1, 3})));
}
/******************************************************************************/
/// BEGIN std::optional<std::string> Gates::get_name(const cmat& U) const
TEST(qpp_Gates_get_name, AllTests) {}
/******************************************************************************/
/// BEGIN template <typename Derived = cmat>
///       Gates::Id(idx D = 2) const
TEST(qpp_Gates_Id, AllTests) {
    EXPECT_EQ(gt.Id(1), cmat::Identity(1, 1));
    EXPECT_EQ(gt.Id(2), cmat::Identity(2, 2));
    EXPECT_EQ(gt.Id(3), cmat::Identity(3, 3));
    EXPECT_EQ(gt.Id(100), cmat::Identity(100, 100));
}
/******************************************************************************/
/// BEGIN cmat Gates::MODMUL(idx a, idx N) const
TEST(qpp_Gates_MODMUL, AllTests) {}
/******************************************************************************/
/// BEGIN cmat Gates::Rn(realT theta, const std::vector<realT>& n) const
TEST(qpp_Gates_Rn, AllTests) {
    // |z0> stays invariant (up to a phase) if rotated by any angle
    // around the Z axis
    EXPECT_NEAR(0, norm(st.pz0 - prj(gt.Rn(2.345, {0, 0, 1}) * st.z0)), 1e-5);

    // |z0> gets a (-1) phase if rotated by 2pi around the X axis
    EXPECT_NEAR(0, norm(st.z0 + gt.Rn(2 * pi, {1, 0, 0}) * st.z0), 1e-5);

    // |z0> gets a (-1) phase if rotated by 2pi around the Y axis
    EXPECT_NEAR(0, norm(st.z0 + gt.Rn(2 * pi, {0, 1, 0}) * st.z0), 1e-5);

    // rotate |x0> by pi/2 around the Z axis, must obtain |y0> (up to a phase)
    EXPECT_NEAR(0, norm(st.py0 - prj(gt.Rn(pi / 2, {0, 0, 1}) * st.x0)), 1e-5);

    // rotate |y0> by pi/2 around the X axis, must obtain |z0> (up to a phase)
    EXPECT_NEAR(0, norm(st.pz0 - prj(gt.Rn(pi / 2, {1, 0, 0}) * st.y0)), 1e-5);

    // rotate |z0> by pi/2 around the Y axis, must obtain |x0> (up to a phase)
    EXPECT_NEAR(0, norm(st.px0 - prj(gt.Rn(pi / 2, {0, 1, 0}) * st.z0)), 1e-5);

    // rotate |y0> by pi around the Z axis, must obtain |y1> (up to a phase)
    EXPECT_NEAR(0, norm(st.py1 - prj(gt.Rn(pi, {0, 0, 1}) * st.y0)), 1e-5);
}
/******************************************************************************/
/// BEGIN cmat Gates::RX(realT theta) const
TEST(qpp_Gates_RX, AllTests) {}
/******************************************************************************/
/// BEGIN cmat Gates::RY(realT theta) const
TEST(qpp_Gates_RY, AllTests) {}
/******************************************************************************/
/// BEGIN cmat Gates::RZ(realT theta) const
TEST(qpp_Gates_RZ, AllTests) {}
/******************************************************************************/
/// BEGIN cmat SWAPd(idx D = 2) const
TEST(qpp_Gates_SWAPd, AllTests) {
    for (idx D = 1; D < 6; ++D) {
        ket psi = randket(D);
        ket phi = randket(D);
        ket a = kron(psi, phi);
        ket b = gt.SWAPd(D) * kron(phi, psi);
        EXPECT_NEAR(0, norm(a - b), 1e-5);
    }
}
/******************************************************************************/
/// BEGIN cmat Gates::Xd(idx D = 2) const
TEST(qpp_Gates_Xd, AllTests) {
    for (idx D = 1; D < 10; ++D) {
        cmat Xd = gt.Xd(D);
        for (idx i = 0; i < D; ++i) {
            ket psi = mket({i}, D);
            ket res = mket({static_cast<idx>((i + 1) % D)}, D);
            EXPECT_NEAR(0, norm(res - Xd * psi), 1e-5);
        }
    }
}
/******************************************************************************/
/// BEGIN cmat Gates::Zd(idx D = 2) const
TEST(qpp_Gates_Zd, AllTests) {
    for (idx D = 1; D < 10; ++D) {
        cmat Zd = gt.Zd(D);
        cplx oD = omega(D);
        for (idx i = 0; i < D; ++i) {
            ket psi = mket({i}, D);
            ket res = static_cast<cplx>(std::pow(oD, i)) * psi;
            EXPECT_NEAR(0, norm(res - Zd * psi), 1e-5);
        }
    }
}
/******************************************************************************/
