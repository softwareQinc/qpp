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

// Unit testing <statistics.h>

/******************************************************************************/
/// BEGIN template<typename Container> double qpp::avg(
///       const std::vector<double>& prob,
///       const Container& X,
///       typename std::enable_if<is_iterable<Container>::value>::type*
///       = nullptr)
TEST(qpp_avg, AllTests)
{

}
/******************************************************************************/
/// BEGIN template<typename Container> double qpp::cor(
///       const dmat& probXY,
///       const Container& X,
///       const Container& Y,
///       typename std::enable_if<is_iterable<Container>::value>::type*
///       = nullptr)
TEST(qpp_cor, AllTests)
{

}
/******************************************************************************/
/// BEGIN template<typename Container> double qpp::cov(
///       const dmat& probXY,
///       const Container& X,
///       const Container& Y,
///       typename std::enable_if<is_iterable<Container>::value>::type*
///       = nullptr)
TEST(qpp_cov, AllTests)
{

}
/******************************************************************************/
/// BEGIN inline std::vector<double> qpp::marginalX(const dmat& probXY)
TEST(qpp_marginalX, AllTests)
{

}
/******************************************************************************/
/// BEGIN inline std::vector<double> qpp::marginalY(const dmat& probXY)
TEST(qpp_marginalY, AllTests)
{

}
/******************************************************************************/
/// BEGIN template<typename Container> double qpp::sigma(
///       const std::vector<double>& prob,
///       const Container& X,
///       typename std::enable_if<is_iterable<Container>::value>::type*
///       = nullptr)
TEST(qpp_sigma, AllTests)
{

}
/******************************************************************************/
/// BEGIN inline std::vector<double> qpp::uniform(idx N)
TEST(qpp_uniform, AllTests)
{

}
/******************************************************************************/
/// BEGIN template<typename Container> double qpp::var(
///       const std::vector<double>& prob,
///       const Container& X,
///       typename std::enable_if<is_iterable<Container>::value>::type*
///       = nullptr)
TEST(qpp_var, AllTests)
{

}
/******************************************************************************/
