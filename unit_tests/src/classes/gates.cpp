/*
 * Quantum++
 *
 * Copyright (c) 2013 - 2016 Vlad Gheorghiu (vgheorgh@gmail.com)
 *
 * This file is part of Quantum++.
 *
 * Quantum++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Quantum++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Quantum++.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <gtest/gtest.h>
#include <qpp.h>
using namespace qpp;

// Write your unit tests here. Some examples are provided below.

/******************************************************************************/
/// BEGIN template<typename Derived> dyn_mat<typename Derived::Scalar>
///        qpp::Gates::CTRL(const Eigen::MatrixBase<Derived>& A,
///        const std::vector<idx>& ctrl,
///        const std::vector<idx>& subsys,
///        idx n, idx d = 2) const
TEST(qpp_Gates_CTRL, Qubits)
{
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
    cmat U = rand<cmat>(2, 2);
    cmat CTRL4 = gt.CTRL(U, {0, 2}, {1}, 3);
    ket psi1 = mket({0, 0, 1});
    ket res1 = mket({0, 0, 1});
    EXPECT_NEAR(0, norm(CTRL4 * psi1 - res1), 1e-10);

    ket psi2 = mket({1, 1, 1});
    ket res2 = kron(st.z1, U * st.z1, st.z1);
    EXPECT_NEAR(0, norm(CTRL4 * psi2 - res2), 1e-10);
}

TEST(qpp_Gates_CTRL, Qudits)
{
    idx D = 3; // qutrits

    // CNOT control-target on 2 qutrits
    cmat CTRL1 = gt.CTRL(gt.Xd(3), {0}, {1}, 2, D);
    EXPECT_NEAR(0, qpp::norm(CTRL1 * qpp::mket({0, 0}, {D, D}) -
                                     qpp::mket({0, 0}, {D,D})), 1e-10);
    EXPECT_NEAR(0, qpp::norm(CTRL1 * qpp::mket({0, 1}, {D, D}) -
                                     qpp::mket({0, 1}, {D, D})), 1e-10);
    EXPECT_NEAR(0, qpp::norm(CTRL1 * qpp::mket({0, 2}, {D, D}) -
                                     qpp::mket({0, 2}, {D, D})), 1e-10);
    EXPECT_NEAR(0, qpp::norm(CTRL1 * qpp::mket({1, 0}, {D, D}) -
                                     qpp::mket({1, 1}, {D, D})), 1e-10);
    EXPECT_NEAR(0, qpp::norm(CTRL1 * qpp::mket({1, 1}, {D, D}) -
                                     qpp::mket({1, 2}, {D, D})), 1e-10);
    EXPECT_NEAR(0, qpp::norm(CTRL1 * qpp::mket({1, 2}, {D, D}) -
                                     qpp::mket({1, 0}, {D, D})), 1e-10);
    EXPECT_NEAR(0, qpp::norm(CTRL1 * qpp::mket({2, 0}, {D, D}) -
                                     qpp::mket({2, 2}, {D, D})), 1e-10);
    EXPECT_NEAR(0, qpp::norm(CTRL1 * qpp::mket({2, 1}, {D, D}) -
                                     qpp::mket({2, 0}, {D, D})), 1e-10);
    EXPECT_NEAR(0, qpp::norm(CTRL1 * qpp::mket({2, 2}, {D, D}) -
                                     qpp::mket({2, 1}, {D, D})), 1e-10);

    // CNOT target-control on 2 qutrits
    cmat CTRL2 = gt.CTRL(gt.Xd(3), {1}, {0}, 2, D);
    EXPECT_NEAR(0, qpp::norm(CTRL2 * qpp::mket({0, 0}, {D, D}) -
                                     qpp::mket({0, 0}, {D,D})), 1e-10);
    EXPECT_NEAR(0, qpp::norm(CTRL2 * qpp::mket({0, 1}, {D, D}) -
                                     qpp::mket({1, 1}, {D, D})), 1e-10);
    EXPECT_NEAR(0, qpp::norm(CTRL2 * qpp::mket({0, 2}, {D, D}) -
                                     qpp::mket({2, 2}, {D, D})), 1e-10);
    EXPECT_NEAR(0, qpp::norm(CTRL2 * qpp::mket({1, 0}, {D, D}) -
                                     qpp::mket({1, 0}, {D, D})), 1e-10);
    EXPECT_NEAR(0, qpp::norm(CTRL2 * qpp::mket({1, 1}, {D, D}) -
                                     qpp::mket({2, 1}, {D, D})), 1e-10);
    EXPECT_NEAR(0, qpp::norm(CTRL2 * qpp::mket({1, 2}, {D, D}) -
                                     qpp::mket({0, 2}, {D, D})), 1e-10);
    EXPECT_NEAR(0, qpp::norm(CTRL2 * qpp::mket({2, 0}, {D, D}) -
                                     qpp::mket({2, 0}, {D, D})), 1e-10);
    EXPECT_NEAR(0, qpp::norm(CTRL2 * qpp::mket({2, 1}, {D, D}) -
                                     qpp::mket({0, 1}, {D, D})), 1e-10);
    EXPECT_NEAR(0, qpp::norm(CTRL2 * qpp::mket({2, 2}, {D, D}) -
                                     qpp::mket({1, 2}, {D, D})), 1e-10);

    // multiple Control-X-X, partial testing
    cmat CTRL3 = gt.CTRL(qpp::kron(gt.Xd(3), gt.Xd(3)), {1, 4}, {2, 3}, 6, 3);
    ket psi1 = mket({0, 1, 2, 2, 1, 1}, {D, D, D, D ,D ,D});
    ket res1 = mket({0, 1, 0, 0, 1, 1}, {D, D, D, D ,D ,D});
    EXPECT_NEAR(0, qpp::norm(CTRL3 * psi1 - res1), 1e-10);

    ket psi2 = mket({0, 1, 2, 2, 2, 1}, {D, D, D, D ,D ,D});
    ket res2 = mket({0, 1, 2, 2, 2, 1}, {D, D, D, D ,D ,D});
    EXPECT_NEAR(0, qpp::norm(CTRL3 * psi2 - res2), 1e-10);

    ket psi3 = mket({1, 2, 1, 0, 2, 2}, {D, D, D, D ,D ,D});
    ket res3 = mket({1, 2, 0, 2, 2, 2}, {D, D, D, D ,D ,D});
    EXPECT_NEAR(0, qpp::norm(CTRL3 * psi3 - res3), 1e-10);
}
/******************************************************************************/
/// BEGIN template<typename Derived>
///       dyn_mat<typename Derived::Scalar> qpp::Gates::expandout(
///       const Eigen::MatrixBase<Derived>& A, idx pos,
///       const std::vector<idx>& dims) const
TEST(qpp_Gates_expandout, AllTests)
{

}
/******************************************************************************/
/// BEGIN cmat qpp::Gates::Fd(idx D) const
TEST(qpp_Gates_Fd, AllTests)
{

}
/******************************************************************************/
/// BEGIN  template<typename Derived = Eigen::MatrixXcd>
///        qpp::Gates::Id(idx D) const
TEST(qpp_Gates_Id, AllTests)
{

}
/******************************************************************************/
/// BEGIN  cmat qpp::Gates::Rn(double theta, const std::vector<double>& n) const
TEST(qpp_Gates_Rn, AllTests)
{

}
/******************************************************************************/
/// BEGIN cmat qpp::Gates::Xd(idx D) const
TEST(qpp_Gates_Xd, AllTests)
{

}
/******************************************************************************/
/// BEGIN cmat qpp::Gates::Zd(idx D) const
TEST(qpp_Gates_Zd, AllTests)
{

}
/******************************************************************************/
