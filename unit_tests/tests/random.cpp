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

// Unit testing <random.h>

/******************************************************************************/
/// BEGIN inline bigint qpp::rand(bigint a, bigint b)
TEST(qpp_rand_integer, AllTests)
{

}
/******************************************************************************/
/// BEGIN inline double qpp::rand(double a, double b)
TEST(qpp_rand_double, AllTests)
{

}
/******************************************************************************/
/// BEGIN template<> inline cmat qpp::rand(
///       idx rows, idx cols, double a, double b)
TEST(qpp_rand_cmat, AllTests)
{

}
/******************************************************************************/
/// BEGIN template<> inline dmat qpp::rand(
///       idx rows, idx cols, double a, double b)
TEST(qpp_rand_dmat, AllTests)
{

}
/******************************************************************************/
/// BEGIN inline cmat qpp::randH(idx D)
TEST(qpp_randH, AllTests)
{

}
/******************************************************************************/
/// BEGIN inline idx qpp::randidx(
///       idx a = std::numeric_limits<idx>::min(),
///       idx b = std::numeric_limits<idx>::max())
TEST(qpp_randidx, AllTests)
{

}
/******************************************************************************/
/// BEGIN inline ket qpp::randket(idx D)
TEST(qpp_randket, AllTests)
{

}
/******************************************************************************/
/// BEGIN inline std::vector<cmat> qpp::randkraus(idx N, idx D)
TEST(qpp_randkraus, AllTests)
{

}
/******************************************************************************/
/// BEGIN inline double qpp::randn(double mean = 0, double sigma = 1)
TEST(qpp_randn_double, AllTests)
{

}
/******************************************************************************/
/// BEGIN template<> inline cmat qpp::randn(
///       idx rows, idx cols, double mean, double sigma)
TEST(qpp_randn_cmat, AllTests)
{

}
/******************************************************************************/
/// BEGIN template<> inline dmat qpp::randn(
///       idx rows, idx cols, double mean, double sigma)
TEST(qpp_randn_dmat, AllTests)
{

}
/******************************************************************************/
/// BEGIN inline std::vector<idx> qpp::randperm(idx n)
TEST(qpp_randperm, AllTests)
{

}
/******************************************************************************/
/// BEGIN inline cmat qpp::randrho(idx D)
TEST(qpp_randrho, AllTests)
{

}
/******************************************************************************/
/// BEGIN inline cmat qpp::randU(idx D)
TEST(qpp_randU, AllTests)
{

}
/******************************************************************************/
/// BEGIN inline cmat qpp::randV(idx Din, idx Dout)
TEST(qpp_randV, AllTests)
{

}
/******************************************************************************/
