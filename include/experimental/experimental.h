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
    const std::vector<idx> _subsys;
    const std::vector<idx> _dims;
public:
    // TODO: exception checking
    MatrixView(const Eigen::MatrixBase<Derived>& exp,
               const std::vector<idx> subsys,
               const std::vector<idx> dims) :
            _data(exp), _subsys(subsys), _dims(dims)
    { }

    MatrixView(const Eigen::MatrixBase<Derived>& exp,
               const std::vector<idx> subsys,
               idx d = 2) :
            MatrixView(exp, subsys, std::vector<idx>(subsys.size(), d))
    { }

    // disable temporaries, as they don't bind to _data via function arguments
    MatrixView(Eigen::MatrixBase<Derived>&& exp) = delete;

    typename Derived::Scalar operator()(std::size_t i, std::size_t j) const
    {
        static idx rowmidx0[maxn], colmidx0[maxn];
        static idx rowmidx1[maxn], colmidx1[maxn];
        internal::_n2multiidx(i, _dims.size(), _dims.data(), rowmidx0);
        internal::_n2multiidx(j, _dims.size(), _dims.data(), colmidx0);
        // TODO: check here
        for(idx i = 0 ; i < _dims.size(); ++i)
        {
            rowmidx1[i] = rowmidx0[_subsys[i]];
            colmidx1[i] = colmidx0[_subsys[i]];
        }
        i = internal::_multiidx2n(rowmidx1, _dims.size(), _dims.data());
        j = internal::_multiidx2n(colmidx1, _dims.size(), _dims.data());

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
MatrixView<Derived> make_MatrixView(Eigen::MatrixBase<Derived>& A,
                                    const std::vector<idx> subsys,
                                    const std::vector<idx> dims)
{
    return MatrixView<Derived>(A, subsys, dims);
}

template<typename Derived>
MatrixView<Derived> make_MatrixView(Eigen::MatrixBase<Derived>& A,
                                    const std::vector<idx> subsys,
                                    idx d = 2)
{
    return MatrixView<Derived>(A, subsys, d);
}


} /* namespace experimental */
} /* namespace qpp */

#endif /* EXPERIMENTAL_EXPERIMENTAL_H_ */
