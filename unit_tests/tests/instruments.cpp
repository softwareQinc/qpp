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

#include "gtest/gtest.h"
#include "qpp.h"

using namespace qpp;

// Unit testing "instruments.h"

/******************************************************************************/
/// BEGIN template<typename Derived> dyn_col_vect<typename Derived::Scalar>
///       qpp::ip(const Eigen::MatrixBase<Derived>& phi,
///       const Eigen::MatrixBase<Derived>& psi,
///       const std::vector<idx>& subsys,
///       const std::vector<idx>& dims)
TEST(qpp_ip, AllTests) {}
/******************************************************************************/
/// BEGIN template<typename Derived> dyn_col_vect<typename Derived::Scalar>
///       qpp::ip(const Eigen::MatrixBase<Derived>& phi,
///       const Eigen::MatrixBase<Derived>& psi,
///       const std::vector<idx>& subsys,
///       idx d = 2)
TEST(qpp_ip_qubits, AllTests) {}
/******************************************************************************/
/// BEGIN template<typename Derived>
///       std::tuple<idx, std::vector<double>, std::vector<cmat>> qpp::measure(
///       const Eigen::MatrixBase<Derived>& A,
///       const cmat& U)
TEST(qpp_measure_full_orthonormal, AllTests) {}
/******************************************************************************/
/// BEGIN template<typename Derived>
///       std::tuple<idx, std::vector<double>, std::vector<cmat>> qpp::measure(
///       const Eigen::MatrixBase<Derived>& A,
///       const cmat& V,
///       const std::vector<idx>& subsys,
///       const std::vector<idx>& dims)
TEST(qpp_measure_rankone, AllTests) {}
/******************************************************************************/
/// BEGIN template<typename Derived>
///       std::tuple<idx, std::vector<double>, std::vector<cmat>> qpp::measure(
///       const Eigen::MatrixBase<Derived>& A,
///       const cmat& V,
///       const std::vector<idx>& subsys,
///       idx d = 2)
TEST(qpp_measure_rankone_qubits, AllTests) {}
/******************************************************************************/
/// BEGIN template<typename Derived>
///       std::tuple<idx, std::vector<double>, std::vector<cmat>>
///       qpp::measure(const Eigen::MatrixBase<Derived>& A,
///       const std::initializer_list<cmat>& Ks)
TEST(qpp_measure_full_kraus_initlist, AllTests) {}
/******************************************************************************/
/// BEGIN template<typename Derived>
///       std::tuple<idx, std::vector<double>, std::vector<cmat>>
///       qpp::measure(const Eigen::MatrixBase<Derived>& A,
///       const std::initializer_list<cmat>& Ks,
///       const std::vector<idx>& subsys,
///       const std::vector<idx>& dims)
TEST(qpp_measure_kraus_initlist, AllTests) {}
/******************************************************************************/
/// BEGIN template<typename Derived>
///       std::tuple<idx, std::vector<double>, std::vector<cmat>>
///       qpp::measure(const Eigen::MatrixBase<Derived>& A,
///       const std::initializer_list<cmat>& Ks,
///       const std::vector<idx>& subsys,
///       idx d = 2)
TEST(qpp_measure_kraus_initlist_qubits, AllTests) {}
/******************************************************************************/
/// BEGIN template<typename Derived>
///       std::tuple<idx, std::vector<double>, std::vector<cmat>>
///       qpp::measure(const Eigen::MatrixBase<Derived>& A,
///       const std::vector<cmat>& Ks)
TEST(qpp_measure_full_kraus_vector, AllTests) {}
/******************************************************************************/
/// BEGIN template<typename Derived>
///       std::tuple<idx, std::vector<double>, std::vector<cmat>>
///       qpp::measure(const Eigen::MatrixBase<Derived>& A,
///       const std::vector<cmat>& Ks,
///       const std::vector<idx>& subsys,
///       const std::vector<idx>& dims)
TEST(qpp_measure_kraus_vector, AllTests) {}
/******************************************************************************/
/// BEGIN template<typename Derived>
///       std::tuple<idx, std::vector<double>, std::vector<cmat>>
///       qpp::measure(const Eigen::MatrixBase<Derived>& A,
///       const std::vector<cmat>& Ks,
///       const std::vector<idx>& subsys,
///       idx d = 2)
TEST(qpp_measure_kraus_vector_qubits, AllTests) {}
/******************************************************************************/
/// BEGIN template<typename Derived> std::tuple<std::vector<idx>, double, cmat>
///       qpp::measure_seq(const Eigen::MatrixBase<Derived>& A,
///       std::vector<idx> subsys,
///       idx d = 2)
TEST(qpp_measure_seq_qubits, AllTests) {}
/******************************************************************************/
/// BEGIN template<typename Derived> std::tuple<std::vector<idx>, double, cmat>
///       qpp::measure_seq(const Eigen::MatrixBase<Derived>& A,
///       std::vector<idx> subsys,
///       std::vector<idx> dims)
TEST(qpp_measure_seq, AllTests) {}
/******************************************************************************/
