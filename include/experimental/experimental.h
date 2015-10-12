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

/**
* \file experimental/experimental.h
* \brief Experimental/test functions/classes
*/

#ifndef EXPERIMENTAL_EXPERIMENTAL_H_
#define EXPERIMENTAL_EXPERIMENTAL_H_

namespace qpp
{
/**
* \namespace qpp::experimental
* \brief Experimental/test functions/classes, do not use or modify
*/
namespace experimental
{

/**
* \brief Partial trace
* \see qpp::ptrace1(), qpp::ptrace2()
*
*  Partial trace of the multi-partite state vector or density matrix
*  over a list of subsystems
*
* \param A Eigen expression
* \param subsys Subsystem indexes
* \param dims Dimensions of the multi-partite system
* \return Partial trace \f$Tr_{subsys}(\cdot)\f$ over the subsytems \a subsys
* in a multi-partite system, as a dynamic matrix
* over the same scalar field as \a A
*/
template<typename Derived>
dyn_mat<typename Derived::Scalar> ptrace(const Eigen::MatrixBase<Derived>& A,
                                         const std::vector<idx>& subsys,
                                         const std::vector<idx>& dims)
{
    const dyn_mat<typename Derived::Scalar>& rA = A;

    // EXCEPTION CHECKS

    // check zero-size
    if (!internal::_check_nonzero_size(rA))
        throw Exception("qpp::ptrace()", Exception::Type::ZERO_SIZE);

    // check that dims is a valid dimension vector
    if (!internal::_check_dims(dims))
        throw Exception("qpp::ptrace()", Exception::Type::DIMS_INVALID);

    // check that subsys are valid
    if (!internal::_check_subsys_match_dims(subsys, dims))
        throw Exception("qpp::ptrace()",
                        Exception::Type::SUBSYS_MISMATCH_DIMS);
    // END EXCEPTION CHECKS

    idx D = static_cast<idx>(rA.rows());
    idx n = dims.size();
    idx nsubsys = subsys.size();
    idx nsubsysbar = n - nsubsys;
    idx dimsubsys = 1;
    for (idx i = 0; i < nsubsys; ++i)
        dimsubsys *= dims[subsys[i]];
    idx dimsubsysbar = D / dimsubsys;

    dyn_mat<typename Derived::Scalar> result =
            dyn_mat<typename Derived::Scalar>(dimsubsysbar, dimsubsysbar);

    //************ ket ************//
    if (internal::_check_cvector(rA)) // we have a ket
    {
        // check that dims match the dimension of A
        if (!internal::_check_dims_match_cvect(dims, rA))
            throw Exception("qpp::ptrace()",
                            Exception::Type::DIMS_MISMATCH_CVECTOR);

        if (subsys.size() == dims.size())
        {
            result(0, 0) = (adjoint(rA) * rA).value();
            return result;
        }

        if (subsys.size() == 0)
            return rA * adjoint(rA);

        return result;
    }
        //************ density matrix ************//
    else if (internal::_check_square_mat(rA)) // we have a density operator
    {
        // check that dims match the dimension of A
        if (!internal::_check_dims_match_mat(dims, rA))
            throw Exception("qpp::ptrace()",
                            Exception::Type::DIMS_MISMATCH_MATRIX);

        if (subsys.size() == dims.size())
        {
            result(0, 0) = rA.trace();
            return result;
        }

        if (subsys.size() == 0)
            return rA;

        return result;
    }
        //************ Exception: not ket nor density matrix ************//
    else
        throw Exception("qpp::ptrace()",
                        Exception::Type::MATRIX_NOT_SQUARE_OR_CVECTOR);
}

/**
* \brief Partial trace
* \see qpp::ptrace1(), qpp::ptrace2()
*
*  Partial trace of the multi-partite state vector or density matrix
*  over a list of subsystems
*
* \param A Eigen expression
* \param subsys Subsystem indexes
* \param d Subsystem dimensions
* \return Partial trace \f$Tr_{subsys}(\cdot)\f$ over the subsytems \a subsys
* in a multi-partite system, as a dynamic matrix
* over the same scalar field as \a A
*/
template<typename Derived>
dyn_mat<typename Derived::Scalar> ptrace(const Eigen::MatrixBase<Derived>& A,
                                         const std::vector<idx>& subsys,
                                         idx d = 2)
{
    const dyn_mat<typename Derived::Scalar>& rA = A;

    // EXCEPTION CHECKS

    // check zero size
    if (!internal::_check_nonzero_size(rA))
        throw Exception("qpp::ptrace()", Exception::Type::ZERO_SIZE);

    // check valid dims
    if (d == 0)
        throw Exception("qpp::ptrace()", Exception::Type::DIMS_INVALID);
    // END EXCEPTION CHECKS

    idx n =
            static_cast<idx>(std::llround(std::log2(rA.rows()) /
                                          std::log2(d)));
    std::vector<idx> dims(n, d); // local dimensions vector

    return ptrace(rA, subsys, dims);
}

} /* namespace experimental */
} /* namespace qpp */

#endif /* EXPERIMENTAL_EXPERIMENTAL_H_ */
