/*
 * Quantum++
 *
 * Copyright (c) 2013 - 2015 Vlad Gheorghiu (vgheorgh@gmail.com)
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
* \file experimental/test.h
* \brief Experimental/test functions/classes
*/

#ifndef EXPERIMENTAL_TEST_H_
#define EXPERIMENTAL_TEST_H_

namespace qpp
{
/**
* \namespace qpp::experimental
* \brief Experimental/test functions/classes, do not use or modify
*/
namespace experimental
{

// sequential partial measurements
/**
* \brief Sequentially measures the part \a subsys
* of the multi-partite state vector or density matrix \a A
* in the computational basis
* \see qpp::measure()
*
* \param A Eigen expression
* \param subsys Subsystem indexes that are measured
* \param dims Dimensions of the multi-partite system
* \return Tuple of: 1. Vector of outcome results of the
* measurement, 2. Outcome probability, and 3. Post-measurement normalized state
*/
template<typename Derived>
std::tuple<std::vector<idx>, double, cmat>
measure_seq(
        const Eigen::MatrixBase<Derived>& A,
        std::vector<idx> subsys,
        std::vector<idx> dims)
{
    const cmat& rA = A;

    // EXCEPTION CHECKS

    // check zero-size
    if (!internal::_check_nonzero_size(rA))
        throw Exception("qpp::measure_seq()", Exception::Type::ZERO_SIZE);

    // check that dimension is valid
    if (!internal::_check_dims(dims))
        throw Exception("qpp::measure_seq()", Exception::Type::DIMS_INVALID);

    // check that dims match rho matrix
    if (!internal::_check_dims_match_mat(dims, rA))
        throw Exception("qpp::measure_seq()",
                Exception::Type::DIMS_MISMATCH_MATRIX);

    // check subsys is valid w.r.t. dims
    if (!internal::_check_subsys_match_dims(subsys, dims))
        throw Exception("qpp::measure_seq()",
                Exception::Type::SUBSYS_MISMATCH_DIMS);

    // END EXCEPTION CHECKS

    std::vector<idx> result;
    double prob = 1;
    cmat outstate = rA;

    // sort subsys in decreasing order
    std::sort(std::begin(subsys), std::end(subsys), std::greater<idx>());

    //************ density matrix or column vector ************//
    if (internal::_check_square_mat(rA) || internal::_check_col_vector(rA))
    {
        while (subsys.size() > 0)
        {
            auto tmp = measure(outstate, gt.Id(dims.at(subsys.at(0))),
                    {subsys.at(0)}, dims);
            result.push_back(std::get<0>(tmp));
            prob *= std::get<1>(tmp).at(std::get<0>(tmp));
            outstate = std::get<2>(tmp).at(std::get<0>(tmp));

            // remove the subsystem
            dims.erase(std::next(dims.begin(), subsys.at(0)));
            subsys.erase(subsys.begin());
        }
    }
    else
        throw Exception("qpp::measure_seq()",
                Exception::Type::MATRIX_NOT_SQUARE_OR_CVECTOR);

    return std::make_tuple(result, prob, outstate);
}

// sequential partial measurements
/**
* \brief Sequentially measures the part \a subsys
* of the multi-partite state vector or density matrix \a A
* in the computational basis
* \see qpp::measure()
*
* \param A Eigen expression
* \param subsys Subsystem indexes that are measured
* \param d Subsystem dimensions
* \return Tuple of: 1. Vector of measurement outcomes, 2.
* Outcome probability, and 3. Post-measurement normalized state
*/
template<typename Derived>
std::tuple<std::vector<idx>, double, cmat>
measure_seq(
        const Eigen::MatrixBase<Derived>& A,
        std::vector<idx> subsys, idx d = 2)
{
    const cmat& rA = A;

    // check zero size
    if (!internal::_check_nonzero_size(rA))
        throw Exception("qpp::measure_seq()", Exception::Type::ZERO_SIZE);

    idx n =
            static_cast<idx>(std::llround(std::log2(rA.rows()) /
                    std::log2(d)));
    std::vector<idx> dims(n, d); // local dimensions vector

    return measure_seq(rA, subsys, dims);
}

} /* namespace experimental */
} /* namespace qpp */

#endif /* EXPERIMENTAL_TEST_H_ */
