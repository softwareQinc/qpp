/*
 * Covariance and smart pointers. Adapted from
 * https://github.com/CppCodeReviewers/Covariant-Return-Types-and-Smart-Pointers
 *
 * The MIT License (MIT)
 *
 * Copyright (c) 2014 C++ Code Revievers
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

#pragma once

#include <memory>

namespace qasmtools {
namespace ast {

template <typename T>
using ptr = std::unique_ptr<T>;

namespace object {
template <typename T>
inline ptr<T> clone(const T& object) {
    using base_type = typename T::base_type;
    static_assert(std::is_base_of<base_type, T>::value,
                  "T object has to derived from T::base_type");
    auto ptrr = static_cast<const base_type&>(object).clone();
    return ptr<T>(static_cast<T*>(ptrr));
}

template <typename T>
struct cloneable {
    using base_type = T;

    virtual ~cloneable() = default;

  protected:
    virtual T* clone() const = 0;

    template <typename X>
    friend ptr<X> object::clone(const X&);
};
} // namespace object

} // namespace ast
} // namespace qasmtools
