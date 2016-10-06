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
TEST(qpp_Gates_CTRL, AllTests)
{

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
