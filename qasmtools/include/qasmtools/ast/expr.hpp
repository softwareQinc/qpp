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
 * \file qasmtools/ast/expr.hpp
 * \brief openQASM expressions
 */

#pragma once

#include "../utils/angle.hpp"
#include "base.hpp"

#include <cmath>
#include <iomanip>

namespace qasmtools {
namespace ast {

/**
 * \brief Enum of binary operators
 */
enum class BinaryOp { Plus, Minus, Times, Divide, Pow };
inline std::ostream& operator<<(std::ostream& os, const BinaryOp& bop) {
    switch (bop) {
        case BinaryOp::Plus:
            os << "+";
            break;
        case BinaryOp::Minus:
            os << "-";
            break;
        case BinaryOp::Times:
            os << "*";
            break;
        case BinaryOp::Divide:
            os << "/";
            break;
        case BinaryOp::Pow:
            os << "^";
            break;
    }
    return os;
}

/**
 * \brief Enum of unary operators
 */
enum class UnaryOp { Neg, Sin, Cos, Tan, Ln, Sqrt, Exp };
inline std::ostream& operator<<(std::ostream& os, const UnaryOp& uop) {
    switch (uop) {
        case UnaryOp::Neg:
            os << "-";
            break;
        case UnaryOp::Sin:
            os << "sin";
            break;
        case UnaryOp::Cos:
            os << "cos";
            break;
        case UnaryOp::Tan:
            os << "tan";
            break;
        case UnaryOp::Ln:
            os << "ln";
            break;
        case UnaryOp::Sqrt:
            os << "sqrt";
            break;
        case UnaryOp::Exp:
            os << "exp";
            break;
    }
    return os;
}

/**
 * \class qasmtools::ast::Expr
 * \brief Base class for openQASM expressions
 */
class Expr : public ASTNode {
  public:
    Expr(parser::Position pos) : ASTNode(pos) {}
    virtual ~Expr() = default;

    /**
     * \brief Evaluate constant expressions
     *
     * All sub-classes must override this
     *
     * \return Returns the value of the expression if it
     *         is constant, or nullopt otherwise
     */
    virtual std::optional<double> constant_eval() const = 0;

    /**
     * \brief Internal pretty-printer with associative context
     *
     * \param ctx Whether the current associative context is ambiguous
     */
    virtual std::ostream& pretty_print(std::ostream& os, bool ctx) const = 0;
    std::ostream& pretty_print(std::ostream& os) const override {
        return pretty_print(os, false);
    }

  protected:
    virtual Expr* clone() const override = 0;
};

/**
 * \class qasmtools::ast::BExpr
 * \brief Class for binary operator expressions
 * \see qasmtools::ast::Expr
 */
class BExpr final : public Expr {
    ptr<Expr> lexp_; ///< the left sub-expression
    BinaryOp op_;    ///< the binary operator
    ptr<Expr> rexp_; ///< the right sub-expression

  public:
    /**
     * \brief Constructs a Binary expression
     *
     * \param pos The source position
     * \param lexp The left sub-expression
     * \param op The binary operator
     * \param rexp The right sub-expression
     */
    BExpr(parser::Position pos, ptr<Expr> lexp, BinaryOp op, ptr<Expr> rexp)
        : Expr(pos), lexp_(std::move(lexp)), op_(op), rexp_(std::move(rexp)) {}

    /**
     * \brief Protected heap-allocated construction
     */
    static ptr<BExpr> create(parser::Position pos, ptr<Expr> lexp, BinaryOp op,
                             ptr<Expr> rexp) {
        return std::make_unique<BExpr>(pos, std::move(lexp), op,
                                       std::move(rexp));
    }

    /**
     * \brief Get the binary operator
     *
     * \return A binary operator enum
     */
    BinaryOp op() const { return op_; }

    /**
     * \brief Get the left sub-expression
     *
     * \return A reference to the left sub-expression
     */
    Expr& lexp() { return *lexp_; }

    /**
     * \brief Get the right sub-expression
     *
     * \return A reference to the right sub-expression
     */
    Expr& rexp() { return *rexp_; }

    /**
     * \brief Set the left sub-expression
     *
     * \param exp The new left sub-expression
     */
    void set_lexp(ptr<Expr> exp) { lexp_ = std::move(exp); }

    /**
     * \brief Set the right sub-expression
     *
     * \param exp The new right sub-expression
     */
    void set_rexp(ptr<Expr> exp) { rexp_ = std::move(exp); }

    std::optional<double> constant_eval() const override {
        auto lexp = lexp_->constant_eval();
        auto rexp = rexp_->constant_eval();

        if (!lexp || !rexp)
            return std::nullopt;

        switch (op_) {
            case BinaryOp::Plus:
                return *lexp + *rexp;
            case BinaryOp::Minus:
                return *lexp - *rexp;
            case BinaryOp::Times:
                return *lexp * *rexp;
            case BinaryOp::Divide:
                return *lexp / *rexp;
            case BinaryOp::Pow:
                return pow(*lexp, *rexp);
            default:
                return 0; // inaccessible
        }
    }
    void accept(Visitor& visitor) override { visitor.visit(*this); }
    std::ostream& pretty_print(std::ostream& os, bool ctx) const override {
        if (ctx) {
            os << "(";
            lexp_->pretty_print(os, true);
            os << op_;
            rexp_->pretty_print(os, true);
            os << ")";
        } else {
            lexp_->pretty_print(os, true);
            os << op_;
            rexp_->pretty_print(os, true);
        }

        return os;
    }

  protected:
    BExpr* clone() const override {
        return new BExpr(pos_, object::clone(*lexp_), op_,
                         object::clone(*rexp_));
    }
};

/**
 * \class qasmtools::ast::UExpr
 * \brief Class for unary operator expressions
 * \see qasmtools::ast::Expr
 */
class UExpr final : public Expr {
    UnaryOp op_;    ///< the unary operator
    ptr<Expr> exp_; ///< the sub-expression

  public:
    /**
     * \brief Constructs a Unary expression
     *
     * \param pos The source position
     * \param op The unary operator
     * \param exp The sub-expression
     */
    UExpr(parser::Position pos, UnaryOp op, ptr<Expr> exp)
        : Expr(pos), op_(op), exp_(std::move(exp)) {}

    /**
     * \brief Protected heap-allocated construction
     */
    static ptr<UExpr> create(parser::Position pos, UnaryOp op, ptr<Expr> exp) {
        return std::make_unique<UExpr>(pos, op, std::move(exp));
    }

    /**
     * \brief Get the operator
     *
     * \return A unary operator enum
     */
    UnaryOp op() const { return op_; }

    /**
     * \brief Get the sub-expression
     *
     * \return A reference to the sub-expression
     */
    Expr& subexp() { return *exp_; }

    /**
     * \brief Set the sub-expression
     *
     * \param exp The new sub-expression
     */
    void set_subexp(ptr<Expr> exp) { exp_ = std::move(exp); }

    std::optional<double> constant_eval() const override {
        auto expr = exp_->constant_eval();

        if (!expr)
            return std::nullopt;

        switch (op_) {
            case UnaryOp::Neg:
                return -(*expr);
            case UnaryOp::Sin:
                return sin(*expr);
            case UnaryOp::Cos:
                return cos(*expr);
            case UnaryOp::Tan:
                return tan(*expr);
            case UnaryOp::Sqrt:
                return sqrt(*expr);
            case UnaryOp::Ln:
                return log(*expr);
            case UnaryOp::Exp:
                return exp(*expr);
            default:
                return 0; // inaccessible
        }
    }
    void accept(Visitor& visitor) override { visitor.visit(*this); }
    std::ostream& pretty_print(std::ostream& os, bool ctx) const override {
        (void) ctx;

        os << op_;
        if (op_ == UnaryOp::Neg)
            exp_->pretty_print(os, true);
        else {
            os << "(";
            exp_->pretty_print(os, false);
            os << ")";
        }

        return os;
    }

  protected:
    UExpr* clone() const override {
        return new UExpr(pos_, op_, object::clone(*exp_));
    }
};

/**
 * \class qasmtools::ast::PiExpr
 * \brief Class for pi constants
 * \see qasmtools::ast::Expr
 */
class PiExpr final : public Expr {

  public:
    /**
     * \brief Construct a Pi expression
     *
     * \param pos The source position
     */
    PiExpr(parser::Position pos) : Expr(pos) {}

    /**
     * \brief Protected heap-allocated construction
     */
    static ptr<PiExpr> create(parser::Position pos) {
        return std::make_unique<PiExpr>(pos);
    }

    std::optional<double> constant_eval() const override { return utils::pi; }
    void accept(Visitor& visitor) override { visitor.visit(*this); }
    std::ostream& pretty_print(std::ostream& os, bool ctx) const override {
        (void) ctx;

        os << "pi";
        return os;
    }

  protected:
    PiExpr* clone() const override { return new PiExpr(pos_); }
};

/**
 * \class qasmtools::ast::IntExpr
 * \brief Class for integer literal expressions
 * \see qasmtools::ast::Expr
 */
class IntExpr final : public Expr {
    int value_; ///< the integer value

  public:
    /**
     * \brief Construct an integer expression
     *
     * \param pos The source position
     * \param val The integer value
     */
    IntExpr(parser::Position pos, int value) : Expr(pos), value_(value) {}

    /**
     * \brief Protected heap-allocated construction
     */
    static ptr<IntExpr> create(parser::Position pos, int value) {
        return std::make_unique<IntExpr>(pos, value);
    }

    /**
     * \brief Get the integer value
     *
     * \return The integer value
     */
    int value() const { return value_; }

    std::optional<double> constant_eval() const override {
        return (double) value_;
    }
    void accept(Visitor& visitor) override { visitor.visit(*this); }
    std::ostream& pretty_print(std::ostream& os, bool ctx) const override {
        (void) ctx;

        os << value_;
        return os;
    }

  protected:
    IntExpr* clone() const override { return new IntExpr(pos_, value_); }
};

/**
 * \class qasmtools::ast::RealExpr
 * \brief Class for floating point literal expressions
 * \see qasmtools::ast::Expr
 */
class RealExpr final : public Expr {
    double value_; ///< the floating point value

  public:
    /**
     * \brief Construct a real-value expression
     *
     * \param pos The source position
     * \param val The floating point value
     */
    RealExpr(parser::Position pos, double value) : Expr(pos), value_(value) {}

    /**
     * \brief Protected heap-allocated construction
     */
    static ptr<RealExpr> create(parser::Position pos, double value) {
        return std::make_unique<RealExpr>(pos, value);
    }

    /**
     * \brief Get the real value
     *
     * \return The floating point value
     */
    double value() const { return value_; }

    std::optional<double> constant_eval() const override { return value_; }
    void accept(Visitor& visitor) override { visitor.visit(*this); }
    std::ostream& pretty_print(std::ostream& os, bool ctx) const override {
        (void) ctx;

        std::streamsize ss = os.precision();
        os << std::setprecision(15) << value_ << std::setprecision(ss);
        return os;
    }

  protected:
    RealExpr* clone() const override { return new RealExpr(pos_, value_); }
};

/**
 * \class qasmtools::ast::VarExpr
 * \brief Class for variable expressions
 * \see qasmtools::ast::Expr
 */
class VarExpr final : public Expr {
    symbol var_; ///< the identifier

  public:
    /**
     * \brief Construct a variable expression
     *
     * \param pos The source position
     * \param var The variable name
     */
    VarExpr(parser::Position pos, symbol var) : Expr(pos), var_(var) {}

    /**
     * \brief Protected heap-allocated construction
     */
    static ptr<VarExpr> create(parser::Position pos, symbol var) {
        return std::make_unique<VarExpr>(pos, var);
    }

    /**
     * \brief Get the variable name
     *
     * \return Constant reference to the name
     */
    const symbol& var() const { return var_; }

    std::optional<double> constant_eval() const override {
        return std::nullopt;
    }
    void accept(Visitor& visitor) override { visitor.visit(*this); }
    std::ostream& pretty_print(std::ostream& os, bool ctx) const override {
        (void) ctx;
        os << var_;
        return os;
    }

  protected:
    VarExpr* clone() const override { return new VarExpr(pos_, var_); }
};

/**
 * \brief Returns an Expr representing the given angle
 *
 * \param theta The angle
 * \return The equivalent QASM expression
 */
inline ptr<Expr> angle_to_expr(const utils::Angle& theta) {
    parser::Position pos;

    if (theta.is_symbolic()) {
        // Angle is of the form pi*(a/b) for a & b integers
        auto [a, b] = *(theta.symbolic_value());

        if (a == 0) {
            return std::make_unique<IntExpr>(IntExpr(pos, 0));
        } else if (a == 1) {
            return std::make_unique<BExpr>(
                pos, std::make_unique<PiExpr>(PiExpr(pos)), BinaryOp::Divide,
                std::make_unique<IntExpr>(IntExpr(pos, b)));
        } else {
            auto subexpr = std::make_unique<BExpr>(
                pos, std::make_unique<PiExpr>(PiExpr(pos)), BinaryOp::Times,
                std::make_unique<IntExpr>(IntExpr(pos, a)));

            return std::make_unique<BExpr>(
                pos, std::move(subexpr), BinaryOp::Divide,
                std::make_unique<IntExpr>(IntExpr(pos, b)));
        }
    } else {
        // Angle is real-valued
        return std::make_unique<RealExpr>(RealExpr(pos, theta.numeric_value()));
    }
}

} // namespace ast
} // namespace qasmtools
