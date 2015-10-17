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
* \class qpp::experimental::MatrixView
* \brief Matrix view class
*/
template<typename Derived>
class MatrixView
{
private:
    const Eigen::MatrixBase<Derived>& _data;
public:
    MatrixView(const Eigen::MatrixBase<Derived>& exp): _data(exp) {}
    // disable temporaries, as they don't bind to _data via function arguments
    MatrixView(Eigen::MatrixBase<Derived>&& exp) = delete;
    typename Derived::Scalar operator()(std::size_t i, std::size_t j) const
    {
        return _data(i, j);
    }
    explicit operator
    Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>()
    const
    {
        return _data;
    }
    virtual ~MatrixView() = default;
};

template<typename Derived>
MatrixView<Derived> make_MatrixView(Eigen::MatrixBase<Derived>& A)
{
    return MatrixView<Derived>(A);
}

} /* namespace experimental */
} /* namespace qpp */

#endif /* EXPERIMENTAL_EXPERIMENTAL_H_ */
