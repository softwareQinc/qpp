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

#include <cmath>
#include <complex>
#include "gtest/gtest.h"
#include "qpp.h"

using namespace qpp;

// Unit testing "classes/gates.h"

/******************************************************************************/
/// BEGIN template<typename Derived> dyn_mat<typename Derived::Scalar>
///        qpp::Gates::CTRL(const Eigen::MatrixBase<Derived>& A,
///        const std::vector<idx>& ctrl,
///        const std::vector<idx>& subsys,
///        idx N,
///        idx d = 2) const
TEST(qpp_Gates_CTRL, Qubits) {
    // CNOT control-target on 2 qubits
    cmat CTRL1 = gt.CTRL(gt.X, {0}, {1}, 2);
    EXPECT_EQ(CTRL1, gt.CNOT);

    // CNOT target-control on 2 qubits
    cmat CTRL2 = gt.CTRL(gt.X, {1}, {0}, 2);
    EXPECT_EQ(CTRL2, gt.CNOTba);

    // TOFOLI
    cmat CTRL3 = gt.CTRL(gt.X, {0, 1}, {2}, 3);
    EXPECT_EQ(CTRL3, gt.TOF);
    CTRL3 = gt.CTRL(gt.X, {0, 1}, {2}, 3, 2); // test non-default args
    EXPECT_EQ(CTRL3, gt.TOF);

    // random gate as multiple control on 2 qubits
    cmat U = randU();
    cmat CTRL4 = gt.CTRL(U, {0, 2}, {1}, 3);
    ket psi1 = mket({0, 0, 1});
    ket res1 = mket({0, 0, 1});
    EXPECT_NEAR(0, norm(CTRL4 * psi1 - res1), 1e-7);

    ket psi2 = mket({1, 1, 1});
    ket res2 = kron(st.z1, U * st.z1, st.z1);
    EXPECT_NEAR(0, norm(CTRL4 * psi2 - res2), 1e-7);
}

TEST(qpp_Gates_CTRL, Qudits) {
    idx D = 3; // qutrits

    // CNOT control-target on 2 qutrits
    cmat CTRL1 = gt.CTRL(gt.Xd(3), {0}, {1}, 2, D);
    EXPECT_NEAR(0, norm(CTRL1 * mket({0, 0}, {D, D}) - mket({0, 0}, {D, D})),
                1e-7);
    EXPECT_NEAR(0, norm(CTRL1 * mket({0, 1}, {D, D}) - mket({0, 1}, {D, D})),
                1e-7);
    EXPECT_NEAR(0, norm(CTRL1 * mket({0, 2}, {D, D}) - mket({0, 2}, {D, D})),
                1e-7);
    EXPECT_NEAR(0, norm(CTRL1 * mket({1, 0}, {D, D}) - mket({1, 1}, {D, D})),
                1e-7);
    EXPECT_NEAR(0, norm(CTRL1 * mket({1, 1}, {D, D}) - mket({1, 2}, {D, D})),
                1e-7);
    EXPECT_NEAR(0, norm(CTRL1 * mket({1, 2}, {D, D}) - mket({1, 0}, {D, D})),
                1e-7);
    EXPECT_NEAR(0, norm(CTRL1 * mket({2, 0}, {D, D}) - mket({2, 2}, {D, D})),
                1e-7);
    EXPECT_NEAR(0, norm(CTRL1 * mket({2, 1}, {D, D}) - mket({2, 0}, {D, D})),
                1e-7);
    EXPECT_NEAR(0, norm(CTRL1 * mket({2, 2}, {D, D}) - mket({2, 1}, {D, D})),
                1e-7);

    // CNOT target-control on 2 qutrits
    cmat CTRL2 = gt.CTRL(gt.Xd(3), {1}, {0}, 2, D);
    EXPECT_NEAR(0, norm(CTRL2 * mket({0, 0}, {D, D}) - mket({0, 0}, {D, D})),
                1e-7);
    EXPECT_NEAR(0, norm(CTRL2 * mket({0, 1}, {D, D}) - mket({1, 1}, {D, D})),
                1e-7);
    EXPECT_NEAR(0, norm(CTRL2 * mket({0, 2}, {D, D}) - mket({2, 2}, {D, D})),
                1e-7);
    EXPECT_NEAR(0, norm(CTRL2 * mket({1, 0}, {D, D}) - mket({1, 0}, {D, D})),
                1e-7);
    EXPECT_NEAR(0, norm(CTRL2 * mket({1, 1}, {D, D}) - mket({2, 1}, {D, D})),
                1e-7);
    EXPECT_NEAR(0, norm(CTRL2 * mket({1, 2}, {D, D}) - mket({0, 2}, {D, D})),
                1e-7);
    EXPECT_NEAR(0, norm(CTRL2 * mket({2, 0}, {D, D}) - mket({2, 0}, {D, D})),
                1e-7);
    EXPECT_NEAR(0, norm(CTRL2 * mket({2, 1}, {D, D}) - mket({0, 1}, {D, D})),
                1e-7);
    EXPECT_NEAR(0, norm(CTRL2 * mket({2, 2}, {D, D}) - mket({1, 2}, {D, D})),
                1e-7);

    // multiple Control-X-X, partial testing
    cmat CTRL3 = gt.CTRL(kron(gt.Xd(3), gt.Xd(3)), {1, 4}, {2, 3}, 6, 3);
    ket psi1 = mket({0, 1, 2, 2, 1, 1}, {D, D, D, D, D, D});
    ket res1 = mket({0, 1, 0, 0, 1, 1}, {D, D, D, D, D, D});
    EXPECT_NEAR(0, norm(CTRL3 * psi1 - res1), 1e-7);

    ket psi2 = mket({0, 1, 2, 2, 2, 1}, {D, D, D, D, D, D});
    ket res2 = mket({0, 1, 2, 2, 2, 1}, {D, D, D, D, D, D});
    EXPECT_NEAR(0, norm(CTRL3 * psi2 - res2), 1e-7);

    ket psi3 = mket({1, 2, 1, 0, 2, 2}, {D, D, D, D, D, D});
    ket res3 = mket({1, 2, 0, 2, 2, 2}, {D, D, D, D, D, D});
    EXPECT_NEAR(0, norm(CTRL3 * psi3 - res3), 1e-7);
}
/******************************************************************************/
/// BEGIN template<typename Derived>
///       dyn_mat<typename Derived::Scalar> qpp::Gates::expandout(
///       const Eigen::MatrixBase<Derived>& A,
///       idx pos,
///       const std::initializer_list<idx>& dims) const
TEST(qpp_Gates_expandout_init_list, AllTests) {
    // single qutrit (degenerate case) random gate expansion
    cmat U = randU(3);
    EXPECT_EQ(gt.expandout(U, 0, {3}), U);
}

/******************************************************************************/
/// BEGIN template<typename Derived>
///       dyn_mat<typename Derived::Scalar> qpp::Gates::expandout(
///       const Eigen::MatrixBase<Derived>& A,
///       idx pos,
///       const std::vector<idx>& dims) const
TEST(qpp_Gates_expandout_vector, AllTests) {
    // single qubit (degenerate case) random gate expansion
    cmat U = randU();
    EXPECT_EQ(gt.expandout(U, 0, std::vector<idx>{2}), U);

    // 4 qutrits, identity on qutrit 3 expansion
    EXPECT_EQ(gt.expandout(gt.Id(3), 2, {3, 3, 3, 3}), gt.Id(81));

    // 3 qubits, X on qubit 2 expansion
    EXPECT_EQ(gt.expandout(gt.X, 1, {2, 2, 2}), kron(gt.Id2, gt.X, gt.Id2));
}
/******************************************************************************/
/// BEGIN template<typename Derived> dyn_mat<typename Derived::Scalar>
///       qpp::expandout(const Eigen::MatrixBase<Derived>& A,
///       idx pos,
///       idx N,
///       idx d = 2) const
TEST(qpp_Gates_expandout_qubits, AllTests) {
    // single qubit (degenerate case) random gate expansion
    cmat U = randU();
    EXPECT_EQ(gt.expandout(U, 0, 1), U);

    // 4 qutrits, identity on qutrit 3 expansion
    EXPECT_EQ(gt.expandout(gt.Id(3), 2, 4, 3), gt.Id(81));

    // 3 qubits, X on qubit 2 expansion
    EXPECT_EQ(gt.expandout(gt.X, 1, 3), kron(gt.Id2, gt.X, gt.Id2));
}
/******************************************************************************/
/// BEGIN cmat qpp::Gates::Fd(idx D = 2) const
TEST(qpp_Gates_Fd, AllTests) {
    EXPECT_NEAR(0, norm(gt.Fd(1) - gt.Id(1)), 1e-7);

    EXPECT_NEAR(0, norm(gt.Fd(2) - gt.H), 1e-7);

    cmat F3(3, 3);
    cplx o3 = omega(3);
    F3 << 1, 1, 1, 1, o3, o3 * o3, 1, o3 * o3, o3;
    F3 /= std::sqrt(3);
    EXPECT_NEAR(0, norm(gt.Fd(3) - F3), 1e-7);

    cmat F4(4, 4);
    cplx o4 = omega(4);
    F4 << 1, 1, 1, 1, 1, o4, o4 * o4, o4 * o4 * o4, 1, o4 * o4, 1, o4 * o4, 1,
        o4 * o4 * o4, o4 * o4, o4;
    F4 /= std::sqrt(4);
    EXPECT_NEAR(0, norm(gt.Fd(4) - F4), 1e-7);
}
/******************************************************************************/
/// BEGIN  template<typename Derived = Eigen::MatrixXcd>
///        qpp::Gates::Id(idx D = 2) const
TEST(qpp_Gates_Id, AllTests) {
    EXPECT_EQ(gt.Id(1), Eigen::MatrixXcd::Identity(1, 1));
    EXPECT_EQ(gt.Id(2), Eigen::MatrixXcd::Identity(2, 2));
    EXPECT_EQ(gt.Id(3), Eigen::MatrixXcd::Identity(3, 3));
    EXPECT_EQ(gt.Id(100), Eigen::MatrixXcd::Identity(100, 100));
}
/******************************************************************************/
/// BEGIN  cmat qpp::Gates::Rn(double theta, const std::vector<double>& n) const
TEST(qpp_Gates_Rn, AllTests) {
    // |z0> stays invariant (up to a phase) if rotated by any angle
    // around the Z axis
    EXPECT_NEAR(0, norm(st.pz0 - prj(gt.Rn(2.345, {0, 0, 1}) * st.z0)), 1e-7);

    // |z0> gets a (-1) phase if rotated by 2pi around the X axis
    EXPECT_NEAR(0, norm(st.z0 + gt.Rn(2 * pi, {1, 0, 0}) * st.z0), 1e-7);

    // |z0> gets a (-1) phase if rotated by 2pi around the Y axis
    EXPECT_NEAR(0, norm(st.z0 + gt.Rn(2 * pi, {0, 1, 0}) * st.z0), 1e-7);

    // rotate |x0> by pi/2 around the Z axis, must obtain |y0> (up to a phase)
    EXPECT_NEAR(0, norm(st.py0 - prj(gt.Rn(pi / 2, {0, 0, 1}) * st.x0)), 1e-7);

    // rotate |y0> by pi/2 around the X axis, must obtain |z0> (up to a phase)
    EXPECT_NEAR(0, norm(st.pz0 - prj(gt.Rn(pi / 2, {1, 0, 0}) * st.y0)), 1e-7);

    // rotate |z0> by pi/2 around the Y axis, must obtain |x0> (up to a phase)
    EXPECT_NEAR(0, norm(st.px0 - prj(gt.Rn(pi / 2, {0, 1, 0}) * st.z0)), 1e-7);

    // rotate |y0> by pi around the Z axis, must obtain |y1> (up to a phase)
    EXPECT_NEAR(0, norm(st.py1 - prj(gt.Rn(pi, {0, 0, 1}) * st.y0)), 1e-7);
}
/******************************************************************************/
/// BEGIN cmat qpp::Gates::Xd(idx D = 2) const
TEST(qpp_Gates_Xd, AllTests) {
    for (idx D = 1; D < 10; ++D) {
        cmat Xd = gt.Xd(D);
        for (idx i = 0; i < D; ++i) {
            ket psi = mket({i}, D);
            ket res = mket({(i + 1) % D}, D);
            EXPECT_NEAR(0, norm(res - Xd * psi), 1e-7);
        }
    }
}
/******************************************************************************/
/// BEGIN cmat qpp::Gates::Zd(idx D = 2) const
TEST(qpp_Gates_Zd, AllTests) {
    for (idx D = 1; D < 10; ++D) {
        cmat Zd = gt.Zd(D);
        cplx oD = omega(D);
        for (idx i = 0; i < D; ++i) {
            ket psi = mket({i}, D);
            ket res = std::pow(oD, i) * psi;
            EXPECT_NEAR(0, norm(res - Zd * psi), 1e-7);
        }
    }
}
/******************************************************************************/
