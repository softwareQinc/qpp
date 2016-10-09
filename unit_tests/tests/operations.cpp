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

// Unit testing <operations.h>

/******************************************************************************/
/// BEGIN template<typename Derived1, typename Derived2>
///       dyn_mat<typename Derived1::Scalar> qpp::apply(
///       const Eigen::MatrixBase<Derived1>& state,
///       const Eigen::MatrixBase<Derived2>& A,
///       const std::vector<idx>& subsys,
///       const std::vector<idx>& dims)
TEST(qpp_apply, AllTests)
{

}
/******************************************************************************/
/// BEGIN template<typename Derived1, typename Derived2>
///       dyn_mat<typename Derived1::Scalar> qpp::apply(
///       const Eigen::MatrixBase<Derived1>& state,
///       const Eigen::MatrixBase<Derived2>& A,
///       const std::vector<idx>& subsys,
///       idx d = 2)
TEST(qpp_apply_qubits, AllTests)
{

}
/******************************************************************************/
/// BEGIN template<typename Derived> cmat qpp::apply(
///       const Eigen::MatrixBase<Derived>& rho, const std::vector<cmat>& Ks)
TEST(qpp_apply_full_kraus, AllTests)
{

}
/******************************************************************************/
/// BEGIN template<typename Derived> cmat qpp::apply(
///       const Eigen::MatrixBase<Derived>& rho,
///       const std::vector<cmat>& Ks,
///       const std::vector<idx>& subsys,
///       const std::vector<idx>& dims)
TEST(qpp_apply_kraus, AllTests)
{

}
/******************************************************************************/
/// BEGIN template<typename Derived> cmat qpp::apply(
///       const Eigen::MatrixBase<Derived>& rho,
///       const std::vector<cmat>& Ks,
///       const std::vector<idx>& subsys,
///       idx d = 2)
TEST(qpp_apply_kraus_qubits, AllTests)
{

}
/******************************************************************************/
/// BEGIN template<typename Derived1, typename Derived2>
///       dyn_mat<typename Derived1::Scalar> qpp::applyCTRL(
///       const Eigen::MatrixBase<Derived1>& state,
///       const Eigen::MatrixBase<Derived2>& A,
///       const std::vector<idx>& ctrl,
///       const std::vector<idx>& subsys,
///       const std::vector<idx>& dims)
TEST(qpp_applyCTRL, NonEmptyControl)
{
    std::vector<idx> dims{2, 2, 2, 2};  // 3 qubits
    idx D = prod(dims);                 // total dimension

    std::vector<idx> ctrl{2, 0};           // where we apply the control
    std::vector<idx> target{1, 3};   // target

    idx Dtarget = 1;                    // dimension of the target subsystems
    for(idx i = 0; i < target.size(); ++i)
    Dtarget *= dims[target[i]];         // compute it here

    // some random n qudit pure state
    ket psi = randket(D);

    cmat rho = psi * adjoint(psi); // the corresponding density matrix
    cmat U = randU(Dtarget);       // some random unitary on the target

    // applyCTRL on pure state
    ket A = applyCTRL(psi, U, ctrl, target, dims);

    // applyCTRL on density matrix
    cmat B = applyCTRL(rho, U, ctrl, target, dims);

    // result when using CTRL-U|psi><psi|CTRL-U^\dagger
    cmat result_psi = A * adjoint(A);
    // result when using CTRL-U(rho)CTRL-U^\dagger
    cmat result_rho = B;

    double res = norm(result_psi - result_rho);
    EXPECT_NEAR (0, res, 1e-10);
}

TEST(qpp_applyCTRL, EmptyControl)
{
    std::vector<idx> dims{2, 2, 2, 2};  // 3 qubits
    idx D = prod(dims);                 // total dimension

    std::vector<idx> ctrl{};           // where we apply the control
    std::vector<idx> target{1, 0, 3};   // target

    idx Dtarget = 1;                    // dimension of the target subsystems
    for(idx i = 0; i < target.size(); ++i)
    Dtarget *= dims[target[i]];         // compute it here

    // some random n qudit pure state
    ket psi = randket(D);

    cmat rho = psi * adjoint(psi); // the corresponding density matrix
    cmat U = randU(Dtarget);       // some random unitary on the target

    // applyCTRL on pure state
    ket A = applyCTRL(psi, U, ctrl, target, dims);

    // applyCTRL on density matrix
    cmat B = applyCTRL(rho, U, ctrl, target, dims);

    // result when using CTRL-U|psi><psi|CTRL-U^\dagger
    cmat result_psi = A * adjoint(A);
    // result when using CTRL-U(rho)CTRL-U^\dagger
    cmat result_rho = B;

    double res = norm(result_psi - result_rho);
    EXPECT_NEAR (0, res, 1e-10);
}
/******************************************************************************/
/// BEGIN template<typename Derived1, typename Derived2>
///       dyn_mat<typename Derived1::Scalar> qpp::applyCTRL(
///       const Eigen::MatrixBase<Derived1>& state,
///       const Eigen::MatrixBase<Derived2>& A,
///       const std::vector<idx>& ctrl,
///       const std::vector<idx>& subsys,
///       idx d = 2)
TEST(qpp_applyCTRL_qubits, AllTests)
{

}
/******************************************************************************/
/// BEGIN inline std::vector<cmat> qpp::choi2kraus(const cmat& A)
TEST(qpp_choi2kraus, AllTests)
{

}
/******************************************************************************/
/// BEGIN inline cmat qpp::choi2super(const cmat& A)
TEST(qpp_choi2super, AllTests)
{

}
/******************************************************************************/
/// BEGIN inline cmat qpp::kraus2choi(const std::vector<cmat>& Ks)
TEST(qpp_kraus2choi, AllTests)
{

}
/******************************************************************************/
/// BEGIN inline cmat qpp::kraus2super(const std::vector<cmat>& Ks)
TEST(qpp_kraus2super, AllTests)
{

}
/******************************************************************************/
/// BEGIN template<typename Derived> dyn_mat<typename Derived::Scalar>
///       qpp::ptrace(const Eigen::MatrixBase<Derived>& A,
///       const std::vector<idx>& subsys,
///       const std::vector<idx>& dims)
TEST(qpp_ptrace, AllTests)
{

}
/******************************************************************************/
/// BEGIN template<typename Derived> dyn_mat<typename Derived::Scalar>
///       qpp::ptrace(const Eigen::MatrixBase<Derived>& A,
///       const std::vector<idx>& subsys,
///       idx d = 2)
TEST(qpp_ptrace_qubits, AllTests)
{

}
/******************************************************************************/
/// BEGIN template<typename Derived> dyn_mat<typename Derived::Scalar>
///       qpp::ptrace1(const Eigen::MatrixBase<Derived>& A,
///       const std::vector<idx>& dims)
TEST(qpp_ptrace1, AllTests)
{

}
/******************************************************************************/
/// BEGIN template<typename Derived> dyn_mat<typename Derived::Scalar>
///       qpp::ptrace2(const Eigen::MatrixBase<Derived>& A,
///       const std::vector<idx>& dims)
TEST(qpp_ptrace2, AllTests)
{

}
/******************************************************************************/
/// BEGIN template<typename Derived> dyn_mat<typename Derived::Scalar>
///       qpp::ptranspose(const Eigen::MatrixBase<Derived>& A,
///       const std::vector<idx>& subsys,
///       const std::vector<idx>& dims)
TEST(qpp_ptranspose, AllTests)
{

}
/******************************************************************************/
/// BEGIN template<typename Derived> dyn_mat<typename Derived::Scalar>
///       qpp::ptranspose(const Eigen::MatrixBase<Derived>& A,
///       const std::vector<idx>& subsys,
///       idx d = 2)
TEST(qpp_ptranspose_qubits, AllTests)
{

}
/******************************************************************************/
/// BEGIN inline cmat qpp::super2choi(const cmat& A)
TEST(qpp_super2choi, AllTests)
{

}
/******************************************************************************/
/// BEGIN template<typename Derived> dyn_mat<typename Derived::Scalar>
///       qpp::syspermute(const Eigen::MatrixBase<Derived>& A,
///       const std::vector<idx>& perm,
///       const std::vector<idx>& dims)
TEST(qpp_syspermute, AllTests)
{

}
/******************************************************************************/
/// BEGIN template<typename Derived> dyn_mat<typename Derived::Scalar>
///       qpp::syspermute(const Eigen::MatrixBase<Derived>& A,
///       const std::vector<idx>& perm,
///       idx d = 2)
TEST(qpp_syspermute_qubits, AllTests)
{

}
/******************************************************************************/
