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
 * \file qasmtools/ast/stmt.hpp
 * \brief openQASM statements
 */

#pragma once

#include "base.hpp"
#include "expr.hpp"
#include "var.hpp"

#include <functional>
#include <vector>

namespace qasmtools {
namespace ast {

/**
 * \class qasmtools::ast::Stmt
 * \brief Base class for openQASM statements
 */
class Stmt : public ASTNode {
  public:
    Stmt(parser::Position pos) : ASTNode(pos) {}
    virtual ~Stmt() = default;

    /**
     * \brief Internal pretty-printer which can suppress the output of the
     * stdlib
     *
     * \param suppress_std Whether to suppress output of the standard library
     */
    virtual std::ostream& pretty_print(std::ostream& os,
                                       bool suppress_std) const = 0;

    std::ostream& pretty_print(std::ostream& os) const override {
        return pretty_print(os, false);
    }

  protected:
    virtual Stmt* clone() const override = 0;
};

/**
 * \class qasmtools::ast::MeasureStmt
 * \brief Class for measurement statements
 * \see qasmtools::ast::Stmt
 */
class MeasureStmt final : public Stmt {
    VarAccess q_arg_; ///< the quantum bit|register
    VarAccess c_arg_; ///< the classical bit|register

  public:
    /**
     * \brief Constructs a measurement statement
     *
     * \param pos The source position
     * \param q_arg Rvalue reference to the quantum argument
     * \param c_arg Rvalue reference to the classical argument
     */
    MeasureStmt(parser::Position pos, VarAccess&& q_arg, VarAccess&& c_arg)
        : Stmt(pos), q_arg_(std::move(q_arg)), c_arg_(std::move(c_arg)) {}

    /**
     * \brief Protected heap-allocated construction
     */
    static ptr<MeasureStmt> create(parser::Position pos, VarAccess&& q_arg,
                                   VarAccess&& c_arg) {
        return std::make_unique<MeasureStmt>(pos, std::move(q_arg),
                                             std::move(c_arg));
    }

    /**
     * \brief Get the quantum argument
     *
     * \return Reference to the quantum argument
     */
    VarAccess& q_arg() { return q_arg_; }

    /**
     * \brief Get the classical argument
     *
     * \return Reference to the classical argument
     */
    VarAccess& c_arg() { return c_arg_; }

    /**
     * \brief Set the quantum argument
     *
     * \param arg Const reference to a new argument
     */
    void set_qarg(const VarAccess& arg) { q_arg_ = arg; }

    /**
     * \brief Set the classical argument
     *
     * \param arg Const reference to a new argument
     */
    void set_carg(const VarAccess& arg) { c_arg_ = arg; }

    void accept(Visitor& visitor) override { visitor.visit(*this); }
    std::ostream& pretty_print(std::ostream& os, bool) const override {
        os << "measure " << q_arg_ << " -> " << c_arg_ << ";\n";
        return os;
    }

  protected:
    MeasureStmt* clone() const override {
        return new MeasureStmt(pos_, VarAccess(q_arg_), VarAccess(c_arg_));
    }
};

/**
 * \class qasmtools::ast::ResetStmt
 * \brief Class for reset statements
 * \see qasmtools::ast::Stmt
 */
class ResetStmt final : public Stmt {
    VarAccess arg_; ///< the qbit|qreg

  public:
    /**
     * \brief Constructs a reset statement
     *
     * \param pos The source position
     * \param arg Rvalue reference to the argument
     */
    ResetStmt(parser::Position pos, VarAccess&& arg)
        : Stmt(pos), arg_(std::move(arg)) {}

    /**
     * \brief Protected heap-allocated construction
     */
    static ptr<ResetStmt> create(parser::Position pos, VarAccess&& arg) {
        return std::make_unique<ResetStmt>(pos, std::move(arg));
    }

    /**
     * \brief Get the argument
     *
     * \return Reference to the argument
     */
    VarAccess& arg() { return arg_; }

    /**
     * \brief Set the argument
     *
     * \param arg Const reference to a new argument
     */
    void set_arg(const VarAccess& arg) { arg_ = arg; }

    void accept(Visitor& visitor) override { visitor.visit(*this); }
    std::ostream& pretty_print(std::ostream& os, bool) const override {
        os << "reset " << arg_ << ";\n";
        return os;
    }

  protected:
    ResetStmt* clone() const override {
        return new ResetStmt(pos_, VarAccess(arg_));
    }
};

/**
 * \class qasmtools::ast::IfStmt
 * \brief Class for if statements
 * \see qasmtools::ast::Stmt
 */
class IfStmt final : public Stmt {
    symbol var_;     ///< classical register name
    int cond_;       ///< value to check against
    ptr<Stmt> then_; ///< statement to be executed if true

  public:
    /**
     * \brief Constructs an if statement
     *
     * \param pos The source position
     * \param var The variable (classical register) being tested
     * \param cond The integer value to test against
     * \param then The statement to execute in the then branch
     */
    IfStmt(parser::Position pos, symbol var, int cond, ptr<Stmt> then)
        : Stmt(pos), var_(var), cond_(cond), then_(std::move(then)) {}

    /**
     * \brief Protected heap-allocated construction
     */
    static ptr<IfStmt> create(parser::Position pos, symbol var, int cond,
                              ptr<Stmt> then) {
        return std::make_unique<IfStmt>(pos, var, cond, std::move(then));
    }

    /**
     * \brief Get the tested variable
     *
     * \return Const reference to the variable name
     */
    const symbol& var() const { return var_; }

    /**
     * \brief Get the integer condition
     *
     * \return The integer tested against
     */
    int cond() const { return cond_; }

    /**
     * \brief Get the then branch
     *
     * \return Reference to the "then" statement
     */
    Stmt& then() { return *then_; }

    /**
     * \brief Set the then branch
     *
     * \param then The new statement
     */
    void set_then(ptr<Stmt> then) { then_ = std::move(then); }

    void accept(Visitor& visitor) override { visitor.visit(*this); }
    std::ostream& pretty_print(std::ostream& os, bool) const override {
        os << "if (" << var_ << "==" << cond_ << ") " << *then_;
        return os;
    }

  protected:
    IfStmt* clone() const override {
        return new IfStmt(pos_, var_, cond_, object::clone(*then_));
    }
};

/**
 * \class qasmtools::ast::Gate
 * \brief Statement sub-class for gate
 */
class Gate : public Stmt {
  public:
    Gate(parser::Position pos) : Stmt(pos) {}
    virtual ~Gate() = default;

  protected:
    virtual Gate* clone() const = 0;
};

/**
 * \class qasmtools::ast::UGate
 * \brief Class for U gates
 * \see qasmtools::ast::Gate
 */
class UGate final : public Gate {
    ptr<Expr> theta_;  ///< theta angle
    ptr<Expr> phi_;    ///< phi angle
    ptr<Expr> lambda_; ///< lambda angle

    VarAccess arg_; ///< quantum bit|register

  public:
    /**
     * \brief Constructs a U gate
     *
     * \param pos The source position
     * \param theta The theta angle
     * \param phi The phi angle
     * \param lambda The lambda angle
     * \param arg Rvalue reference to the quantum argument
     */
    UGate(parser::Position pos, ptr<Expr> theta, ptr<Expr> phi,
          ptr<Expr> lambda, VarAccess&& arg)
        : Gate(pos), theta_(std::move(theta)), phi_(std::move(phi)),
          lambda_(std::move(lambda)), arg_(std::move(arg)) {}

    /**
     * \brief Protected heap-allocated construction
     */
    static ptr<UGate> create(parser::Position pos, ptr<Expr> theta,
                             ptr<Expr> phi, ptr<Expr> lambda, VarAccess&& arg) {
        return std::make_unique<UGate>(pos, std::move(theta), std::move(phi),
                                       std::move(lambda), std::move(arg));
    }

    /**
     * \brief Get the theta angle
     *
     * \return Reference to the angle expression
     */
    Expr& theta() { return *theta_; }

    /**
     * \brief Get the phi angle
     *
     * \return Reference to the angle expression
     */
    Expr& phi() { return *phi_; }

    /**
     * \brief Get the lambda angle
     *
     * \return Reference to the angle expression
     */
    Expr& lambda() { return *lambda_; }

    /**
     * \brief Get the argument
     *
     * \return Reference to the quantum argument
     */
    VarAccess& arg() { return arg_; }

    /**
     * \brief Set the theta angle
     *
     * \param theta The new angle expression
     */
    void set_theta(ptr<Expr> theta) { theta_ = std::move(theta); }

    /**
     * \brief Set the phi angle
     *
     * \param theta The new angle expression
     */
    void set_phi(ptr<Expr> phi) { phi_ = std::move(phi); }

    /**
     * \brief Set the lambda angle
     *
     * \param theta The new angle expression
     */
    void set_lambda(ptr<Expr> lambda) { lambda_ = std::move(lambda); }

    /**
     * \brief Set the argument
     *
     * \param arg The new argument
     */
    void set_arg(const VarAccess& arg) { arg_ = arg; }

    void accept(Visitor& visitor) override { visitor.visit(*this); }
    std::ostream& pretty_print(std::ostream& os, bool) const override {
        os << "U(" << *theta_ << "," << *phi_ << "," << *lambda_ << ") " << arg_
           << ";\n";
        return os;
    }

  protected:
    UGate* clone() const override {
        return new UGate(pos_, object::clone(*theta_), object::clone(*phi_),
                         object::clone(*lambda_), VarAccess(arg_));
    }
};

/**
 * \class qasmtools::ast::CNOTGate
 * \brief Class for CX gates
 * \see qasmtools::ast::Gate
 */
class CNOTGate final : public Gate {
    VarAccess ctrl_; ///< control qubit|qreg
    VarAccess tgt_;  ///< target qubit|qreg

  public:
    /**
     * \brief Constructs a CNOT gate
     *
     * \param pos The source position
     * \param ctrl Rvalue reference to the control argument
     * \param tgt Rvalue reference to the target argument
     */
    CNOTGate(parser::Position pos, VarAccess&& ctrl, VarAccess&& tgt)
        : Gate(pos), ctrl_(std::move(ctrl)), tgt_(std::move(tgt)) {}

    /**
     * \brief Protected heap-allocated construction
     */
    static ptr<CNOTGate> create(parser::Position pos, VarAccess&& ctrl,
                                VarAccess&& tgt) {
        return std::make_unique<CNOTGate>(pos, std::move(ctrl), std::move(tgt));
    }

    /**
     * \brief Get the control argument
     *
     * \return Reference to the quantum argument
     */
    VarAccess& ctrl() { return ctrl_; }

    /**
     * \brief Get the target argument
     *
     * \return Reference to the quantum argument
     */
    VarAccess& tgt() { return tgt_; }

    /**
     * \brief Set the control argument
     *
     * \param ctrl The new argument
     */
    void set_ctrl(const VarAccess& ctrl) { ctrl_ = ctrl; }

    /**
     * \brief Set the target argument
     *
     * \param tgt The new argument
     */
    void set_tgt(const VarAccess& tgt) { tgt_ = tgt; }

    void accept(Visitor& visitor) override { visitor.visit(*this); }
    std::ostream& pretty_print(std::ostream& os, bool) const override {
        os << "CX " << ctrl_ << "," << tgt_ << ";\n";
        return os;
    }

  protected:
    CNOTGate* clone() const override {
        return new CNOTGate(pos_, VarAccess(ctrl_), VarAccess(tgt_));
    }
};

/**
 * \class qasmtools::ast::BarrierGate
 * \brief Class for barrier gates
 * \see qasmtools::ast::Gate
 */
class BarrierGate final : public Gate {
    std::vector<VarAccess> args_; ///< list of quantum bits|registers

  public:
    /**
     * \brief Constructs a barrier gate
     *
     * \param pos The source position
     * \param args Rvalue reference to a list of arguments
     */
    BarrierGate(parser::Position pos, std::vector<VarAccess>&& args)
        : Gate(pos), args_(std::move(args)) {}

    /**
     * \brief Protected heap-allocated construction
     */
    static ptr<BarrierGate> create(parser::Position pos,
                                   std::vector<VarAccess>&& args) {
        return std::make_unique<BarrierGate>(pos, std::move(args));
    }

    /**
     * \brief Get the number of arguments
     *
     * \return The number of arguments
     */
    int num_args() const { return static_cast<int>(args_.size()); }

    /**
     * \brief Get the list of arguments
     *
     * \return Reference to the list of arguments
     */
    std::vector<VarAccess>& args() { return args_; }

    /**
     * \brief Get the ith argument
     *
     * \param i The number of the argument (0 indexed)
     * \return Reference to the ith argument
     */
    VarAccess& arg(int i) { return args_[i]; }

    /**
     * \brief Apply a function to each argument
     *
     * \param f Void function accepting a reference to the argument
     */
    void foreach_arg(std::function<void(VarAccess&)> f) {
        for (auto it = args_.begin(); it != args_.end(); it++)
            f(*it);
    }

    /**
     * \brief Set the ith argument
     *
     * \param i The number of the argument (0 indexed)
     * \param arg The new argument
     */
    void set_arg(int i, const VarAccess& arg) { args_[i] = arg; }

    void accept(Visitor& visitor) override { visitor.visit(*this); }
    std::ostream& pretty_print(std::ostream& os, bool) const override {
        os << "barrier ";
        for (auto it = args_.begin(); it != args_.end(); it++) {
            os << (it == args_.begin() ? "" : ",") << *it;
        }
        os << ";\n";
        return os;
    }

  protected:
    BarrierGate* clone() const override {
        return new BarrierGate(pos_, std::vector<VarAccess>(args_));
    }
};

/**
 * \class qasmtools::ast::DeclaredGate
 * \brief Class for declared gate applications
 * \see qasmtools::ast::Gate
 */
class DeclaredGate final : public Gate {
    symbol name_;                   ///< gate identifier
    std::vector<ptr<Expr>> c_args_; ///< list of classical arguments
    std::vector<VarAccess> q_args_; ///< list of quantum arguments

  public:
    /**
     * \brief Constructs an application of a declared gate
     *
     * \param pos The source position
     * \param name The gate name
     * \param c_args Rvalue reference to a list of classical arguments
     * \param q_args Rvalue reference to a list of quantum arguments
     */
    DeclaredGate(parser::Position pos, symbol name,
                 std::vector<ptr<Expr>>&& c_args,
                 std::vector<VarAccess>&& q_args)
        : Gate(pos), name_(name), c_args_(std::move(c_args)),
          q_args_(std::move(q_args)) {}

    /**
     * \brief Protected heap-allocated construction
     */
    static ptr<DeclaredGate> create(parser::Position pos, symbol name,
                                    std::vector<ptr<Expr>>&& c_args,
                                    std::vector<VarAccess>&& q_args) {
        return std::make_unique<DeclaredGate>(pos, name, std::move(c_args),
                                              std::move(q_args));
    }

    /**
     * \brief Get the gate name
     *
     * \return Const reference to the gate name
     */
    const symbol& name() const { return name_; }

    /**
     * \brief Get the number of classical arguments
     *
     * \return The number of arguments
     */
    int num_cargs() const { return static_cast<int>(c_args_.size()); }

    /**
     * \brief Get the number of quantum arguments
     *
     * \return The number of arguments
     */
    int num_qargs() const { return static_cast<int>(q_args_.size()); }

    /**
     * \brief Get the ith classical argument
     *
     * \param i The number of the argument, 0-indexed
     * \return Reference to an expression
     */
    Expr& carg(int i) { return *(c_args_[i]); }

    /**
     * \brief Get the ith quantum argument
     *
     * \param i The number of the argument, 0-indexed
     * \return Reference to the argument
     */
    VarAccess& qarg(int i) { return q_args_[i]; }

    /**
     * \brief Get the list of quantum arguments
     *
     * \return Reference to the list of arguments
     */
    std::vector<VarAccess>& qargs() { return q_args_; }

    /**
     * \brief Apply a function to each classical argument
     *
     * \param f Void function accepting an expression reference
     */
    void foreach_carg(std::function<void(Expr&)> f) {
        for (auto it = c_args_.begin(); it != c_args_.end(); it++)
            f(**it);
    }

    /**
     * \brief Apply a function to each quantum argument
     *
     * \param f Void function accepting a reference to an argument
     */
    void foreach_qarg(std::function<void(VarAccess&)> f) {
        for (auto it = q_args_.begin(); it != q_args_.end(); it++)
            f(*it);
    }

    /**
     * \brief Set the ith classical argument
     *
     * \param i The number of the argument, 0-indexed
     * \param expr An expression giving the new argument
     */
    void set_carg(int i, ptr<Expr> expr) { c_args_[i] = std::move(expr); }

    /**
     * \brief Set the ith quantum argument
     *
     * \param i The number of the argument, 0-indexed
     * \param arg The new argument
     */
    void set_qarg(int i, const VarAccess& arg) { q_args_[i] = arg; }

    void accept(Visitor& visitor) override { visitor.visit(*this); }
    std::ostream& pretty_print(std::ostream& os, bool) const override {
        os << name_;
        if (c_args_.size() > 0) {
            os << "(";
            for (auto it = c_args_.begin(); it != c_args_.end(); it++) {
                os << (it == c_args_.begin() ? "" : ",") << **it;
            }
            os << ")";
        }
        os << " ";
        for (auto it = q_args_.begin(); it != q_args_.end(); it++) {
            os << (it == q_args_.begin() ? "" : ",") << *it;
        }
        os << ";\n";
        return os;
    }

  protected:
    DeclaredGate* clone() const override {
        std::vector<ptr<Expr>> c_tmp;
        for (auto it = c_args_.begin(); it != c_args_.end(); it++) {
            c_tmp.emplace_back(object::clone(**it));
        }

        return new DeclaredGate(pos_, name_, std::move(c_tmp),
                                std::vector<VarAccess>(q_args_));
    }
};

} // namespace ast
} // namespace qasmtools
