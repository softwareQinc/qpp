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

#include "gtest/gtest.h"
#include "qpp.h"

using namespace qpp;

// Unit testing "entropies.h"

/******************************************************************************/
/// BEGIN template<typename Derived> double qpp::entropy(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_entropy_matrix, AllTests)
{

}
/******************************************************************************/
/// BEGIN inline double qpp::entropy(const std::vector<double>& prob)
TEST(qpp_entropy_vector, AllTests)
{

}
/******************************************************************************/
/// BEGIN template<typename Derived> double qpp::qmutualinfo(
///       const Eigen::MatrixBase<Derived>& A,
///       const std::vector<idx>& subsysA,
///       const std::vector<idx>& subsysB,
///       const std::vector<idx>& dims)
TEST(qpp_qmutualinfo, AllTests)
{

}
/******************************************************************************/
/// BEGIN template<typename Derived> double qpp::qmutualinfo(
///       const Eigen::MatrixBase<Derived>& A,
///       const std::vector<idx>& subsysA,
///       const std::vector<idx>& subsysB,
///       idx d = 2)
TEST(qpp_mutualinfo_qubits, AllTests)
{

}
/******************************************************************************/
/// BEGIN template<typename Derived> double qpp::renyi(
///       const Eigen::MatrixBase<Derived>& A, double alpha)
TEST(qpp_renyi_matrix, AllTests)
{

}
/******************************************************************************/
/// BEGIN qpp::inline double renyi(const std::vector<double>& prob,
///       double alpha)
TEST(qpp_renyi_vector, AllTests)
{

}
/******************************************************************************/
/// BEGIN template<typename Derived> double qpp::tsallis(
///       const Eigen::MatrixBase<Derived>& A, double q)
TEST(qpp_tsallis_matrix, AllTests)
{

}
/******************************************************************************/
/// BEGIN inline double qpp::tsallis(const std::vector<double>& prob, double q)
TEST(qpp_tsallis_vector, AllTests)
{

}
/******************************************************************************/
