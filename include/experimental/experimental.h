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

// TODO: exception checking
// TODO: check for noexcept
// TODO: write the documentation
// TODO: test

namespace qpp
{
/**
* \namespace qpp::experimental
* \brief Experimental/test functions/classes, do not use or modify
*/
namespace experimental
{


/**
* \class qpp::experimental::MatrixViewBase
* \brief Matrix view base class for all other views
*/

template<typename Derived>
class MatrixViewBase
{
protected:
    const Eigen::MatrixBase<Derived>& _viewA;
private:
    idx _rows, _cols;
public:
    MatrixViewBase(const Eigen::MatrixBase<Derived>& A) :
            _viewA(A),
            _rows(static_cast<idx>(A.rows())),
            _cols(static_cast<idx>(A.cols()))
    { }

    // disable temporaries
    MatrixViewBase(const Eigen::MatrixBase<Derived>&& A) = delete;

    // getters and setters
    idx rows() const noexcept
    {
        return _rows;
    }

    idx cols() const noexcept
    {
        return _cols;
    }

    // raw reference
    const Eigen::MatrixBase<Derived>& get_ref() const noexcept
    {
        return _viewA;
    }

    // get a copy of the reference
    Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
    get_copy() const
    {
        return static_cast<Eigen::Matrix<typename Derived::Scalar,
                Eigen::Dynamic, Eigen::Dynamic>>(*this);
    }

    // conversions

    // allow only explicit conversions
    explicit operator
    Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>()
    const
    {
        Eigen::Matrix<
                typename Derived::Scalar,
                Eigen::Dynamic,
                Eigen::Dynamic
        > result(_rows, _cols);

#pragma omp parallel for collapse(2)
        for (idx j = 0; j < _cols; ++j)
            for (idx i = 0; i < _rows; ++i)
                result(i, j) = this->operator()(i, j);

        return result;
    }

    // check out of bounds access, slower than raw operator()(idx, idx)
    typename Derived::Scalar at(idx i, idx j = 0) const
    {
        // EXCEPTION CHECKS

        // END EXCEPTION CHECKS
        return this->operator()(i, j);
    }

    // abstract functions

    // doesn't throw in case of out of bound access
    virtual typename Derived::Scalar operator()(idx i, idx j = 0) const = 0;

    virtual ~MatrixViewBase() = default;
};

// default implementation
template<typename Derived>
typename Derived::Scalar
MatrixViewBase<Derived>::operator()(idx i, idx j) const
{
    assert(i < this->_rows && j < this->_cols);
    return this->_viewA(i, j);
}

/**
* \class qpp::experimental::MatrixView
* \brief Matrix view class, maps between a matrix and a multi-dimensional array
*/

// The qubit \a perm[\a i] is permuted to the location \a i.
template<typename Derived>
class MatrixView : public MatrixViewBase<Derived>
{
private:
    const std::vector<idx> _perm;
    const std::vector<idx> _dims;
public:
    MatrixView(const Eigen::MatrixBase<Derived>& A,
               const std::vector<idx> perm,
               const std::vector<idx> dims) :
            MatrixViewBase<Derived>(A), _perm(perm), _dims(dims)
    { }

    MatrixView(const Eigen::MatrixBase<Derived>& A,
               const std::vector<idx> perm,
               idx d = 2) :
            MatrixView(A, perm, std::vector<idx>(perm.size(), d))
    { }

    // additional getters
    std::vector<idx> perm() const noexcept
    {
        return _perm;
    }

    std::vector<idx> dims() const noexcept
    {
        return _dims;
    }

    // interface implementation
    typename Derived::Scalar operator()(idx i, idx j = 0) const override
    {
        idx Crowmidx[maxn], Ccolmidx[maxn];
        idx Crowmidx_shuffled[maxn], Ccolmidx_shuffled[maxn];
        idx Cdims_shuffled[maxn];

        internal::_n2multiidx(i, _dims.size(), _dims.data(), Crowmidx);
        internal::_n2multiidx(j, _dims.size(), _dims.data(), Ccolmidx);

        // TODO: check this loop
        for (idx k = 0; k < _dims.size(); ++k)
        {
            Cdims_shuffled[_perm[k]] = _dims[k];
            Crowmidx_shuffled[_perm[k]] = Crowmidx[k];
            Ccolmidx_shuffled[_perm[k]] = Ccolmidx[k];
        }

        i = internal::_multiidx2n(Crowmidx_shuffled, _dims.size(),
                                  Cdims_shuffled);
        j = internal::_multiidx2n(Ccolmidx_shuffled, _dims.size(),
                                  Cdims_shuffled);

        return this->_viewA(i, j);
        /*return MatrixViewBase<Derived>::operator()(i, j);*/
    }
};

template<typename Derived>
MatrixView<Derived> make_MatrixView(const Eigen::MatrixBase<Derived>& A,
                                    const std::vector<idx> perm,
                                    const std::vector<idx> dims)
{
    return MatrixView<Derived>(A, perm, dims);
}


template<typename Derived>
MatrixView<Derived> make_MatrixView(const Eigen::MatrixBase<Derived>& A,
                                    const std::vector<idx> perm,
                                    idx d = 2)
{
    return MatrixView<Derived>(A, perm, d);
}

// disable temporaries
template<typename Derived>
MatrixView<Derived> make_MatrixView(const Eigen::MatrixBase<Derived>&& A,
                                    const std::vector<idx> perm,
                                    const std::vector<idx> dims) = delete;

// disable temporaries
template<typename Derived>
MatrixView<Derived> make_MatrixView(const Eigen::MatrixBase<Derived>&& A,
                                    const std::vector<idx> perm,
                                    idx d = 2) = delete;

} /* namespace experimental */
} /* namespace qpp */

#endif /* EXPERIMENTAL_EXPERIMENTAL_H_ */
