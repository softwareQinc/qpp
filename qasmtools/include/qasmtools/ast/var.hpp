/*
 * This file is part of qasmtools.
 *
 * Copyright (c) 2019 - 2021 softwareQ Inc. All rights reserved.
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
 * \file qasmtools/ast/var.hpp
 * \brief openQASM variable utilities
 */

#pragma once

#include "base.hpp"

#include <cstddef>
#include <optional>

namespace qasmtools {
namespace ast {

/**
 * \class qasmtools::ast::VarAccess
 * \brief Class for variable accesses
 *
 * Represents accesses into a register by the register name and an optional
 * offset or index into the register. If the offset is empty, the entire
 * register is the access -- e.g. in gates applied in parallel across registers.
 *
 * As leaf nodes that do not usually need to be used in polymorphic contexts,
 * variable accesses are the only nodes which by convention are **NOT**
 * allocated on the heap.
 */
class VarAccess final : public ASTNode {
    symbol var_;                ///< the identifier
    std::optional<int> offset_; ///< optional offset into a register variable

  public:
    friend std::hash<VarAccess>; ///< Hash function

    /**
     * \brief Construct a variable access
     *
     * \param pos The source position
     * \param var The register name
     * \param offset Optional integer offset into the register (default =
     * std::nullopt)
     */
    VarAccess(parser::Position pos, symbol var,
              std::optional<int> offset = std::nullopt)
        : ASTNode(pos), var_(var), offset_(offset) {}

    /**
     * \brief Copy constructor
     */
    VarAccess(const VarAccess& va)
        : ASTNode(va.pos_), var_(va.var_), offset_(va.offset_) {}

    /**
     * \brief Get the register name
     *
     * return Const reference to the register name
     */
    const symbol& var() const { return var_; }

    /**
     * \brief Get the offset
     *
     * return std::optional integer offset
     */
    std::optional<int> offset() const { return offset_; }

    /**
     * \brief Copy assignment overload
     */
    VarAccess& operator=(const VarAccess& v) {
        var_ = v.var_;
        offset_ = v.offset_;
        return *this;
    }

    /**
     * \brief Equal operator overload
     */
    bool operator==(const VarAccess& v) const {
        return var_ == v.var_ && offset_ == v.offset_;
    }

    /**
     * \brief Less operator overload
     *
     * Used to allow variable accesses as keys in ordered maps
     */
    bool operator<(const VarAccess& v) const {
        if (var_ == v.var_)
            return offset_ < v.offset_;
        else
            return var_ < v.var_;
    }

    /**
     * \brief Check whether the variable access contains another
     *
     * A variable access u contains v if u == v or if u is a register
     * and v is an offset into that register. Mainly useful for determining
     * dependencies between gates in the "sugared" source
     *
     * \param v Const reference to a variable access
     * \return true if the variable access contains v
     */
    bool contains(const VarAccess& v) const {
        if (offset_)
            return *this == v;
        else
            return v.var_ == var_;
    }

    /**
     * \brief Return the root of a variable access
     *
     * Strips any dereferences and returns a new variable access.
     * Satisfies root(v).contains(v) == true
     *
     * \param v Const reference to a variable access
     * \return var access for the root variable
     */
    VarAccess root() const { return VarAccess(pos_, var_); }

    friend std::size_t hash_value(const VarAccess& v) {
        std::size_t lhs = std::hash<symbol>{}(v.var_);
        lhs ^= std::hash<std::optional<int>>{}(v.offset_) + 0x9e3779b9 +
               (lhs << 6) + (lhs >> 2);
        return lhs;
    }

    void accept(Visitor& visitor) override { visitor.visit(*this); }
    std::ostream& pretty_print(std::ostream& os) const override {
        os << var_;
        if (offset_)
            os << "[" << *offset_ << "]";
        return os;
    }
  protected:
    VarAccess* clone() const override {
        return new VarAccess(pos_, var_, offset_);
    }
};

} // namespace ast
} // namespace qasmtools

namespace std {
/**
 * \brief Hash function for variable accesses
 *
 * Allows variable accesses to be used as keys in std::unordered_map.
 * Implementation and magic numbers taken from boost::hash_combine.
 */
template <>
struct hash<qasmtools::ast::VarAccess> {
    std::size_t operator()(const qasmtools::ast::VarAccess& v) const {
        std::size_t lhs = std::hash<std::string>{}(v.var_);
        lhs ^= std::hash<std::optional<int>>{}(v.offset_) + 0x9e3779b9 +
               (lhs << 6) + (lhs >> 2);
        return lhs;
    }
};
} // namespace std
