/*
 * This file is part of Quantum++.
 *
 * Copyright (c) 2017 - 2025 softwareQ Inc. All rights reserved.
 *
 * MIT License
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

/**
 * \file qpp/internal/classes/labelled_vector_proxy.hpp
 * \brief qpp::internal::LabelledVectorProxy
 */

#ifndef QPP_INTERNAL_CLASSES_LABELLED_VECTOR_PROXY_HPP_
#define QPP_INTERNAL_CLASSES_LABELLED_VECTOR_PROXY_HPP_

#include <cstddef>
#include <numeric>
#include <type_traits>
#include <utility>
#include <vector>

#include <iostream>

#include "qpp/classes/exception.hpp"

namespace qpp {
namespace internal {

///< Mutable view into a vector via a labelling, i.e., v[label[i]]
template <class T, bool is_const>
class LabelledVectorProxy {
    using VecType =
        std::conditional_t<is_const, const std::vector<T>, std::vector<T>>;
    VecType& data_{};
    std::vector<std::size_t> label_{};

  public:
    LabelledVectorProxy(std::vector<T>& data, std::vector<std::size_t> label)
        : data_{data}, label_{std::move(label)} {}

    LabelledVectorProxy(std::vector<T>& data) : data_{data} {
        label_.resize(data_.size());
        std::iota(label_.begin(), label_.end(), 0);
    }

    const T& operator[](std::size_t i) const {
        // EXCEPTION CHECKS
        if (i + 1 > label_.size()) {
            throw exception::OutOfRange{
                "qpp::internal::LabelledVectorProxy::operator[]() const", "i"};
        }
        if (label_[i] + 1 > data_.size()) {
            throw exception::OutOfRange{
                "qpp::internal::LabelledVectorProxy::operator[]() const", "i"};
        }
        // END EXCEPTION CHECKS

        return data_[label_[i]];
    }

    // mutable operator[], enabled only if is_const == false
    template <bool B = is_const, typename = std::enable_if_t<!B>>
    T& operator[](std::size_t i) {
        // FIXME: operator[] exceptions in LabelledVectorProxy
        std::cout << "i: " << i << '\n';
        std::cout << "label[i]: " << label_[i] << '\n';
        std::cout << "label size: " << label_.size() << '\n';
        if (i + 1 > label_.size()) {
            throw exception::OutOfRange{
                "qpp::internal::LabelledVectorProxy::operator[]()", "i"};
        }
        if (label_[i] + 1 > data_.size()) {
            throw exception::OutOfRange{
                "qpp::internal::LabelledVectorProxy::operator[]()", "i"};
        }
        return data_[label_[i]];
    }

    // void set_val(std::size_t i, T val) { this->operator[](i) = val; }
};

} /* namespace internal */
} /* namespace qpp */

#endif /* QPP_INTERNAL_CLASSES_LABELLED_VECTOR_PROXY_HPP_ */
