/*
 * This file is part of qasmtools.
 *
 * Copyright (c) 2019 - 2024 softwareQ Inc. All rights reserved.
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
 * \file qasmtools/ast/replacer.hpp
 * \brief Node replacement for syntax trees
 */

#ifndef QASMTOOLS_AST_REPLACER_HPP_
#define QASMTOOLS_AST_REPLACER_HPP_

#include <unordered_map>

#include "program.hpp"
#include "visitor.hpp"

namespace qasmtools {
namespace ast {
/**
 * \class qasmtools::ast::Replacer
 * \brief Generic complete traversal with post-order node replacement
 * \see qasmtools::ast::Visitor
 *
 * The replacer provides a visitor-style interface where the visitor
 * of a node optionally returns a node of the same base type (in the
 * case of statements or gates, a possibly empty list of statements or
 * gates). The visitor logic then replaces the visited node with the
 * returned node, if non-empty, or otherwise leaves the node unchanged.
 *
 * Standard usage is to derive from Replacer and override the replace
 * methods only for the relevant nodes. The replace method is called
 * post-order -- after visiting all children -- and overriding a replace
 * method does not kill traversal to the node's children. To stop
 * descending into the children of a node, the node's visit overload
 * can be overridden.
 */
class Replacer : public Visitor {
    std::optional<VarAccess> replacement_var_;
    std::optional<ptr<Expr>> replacement_expr_;
    std::optional<std::list<ptr<Stmt>>> replacement_stmts_;
    std::optional<std::list<ptr<Gate>>> replacement_gates_;

  public:
    // Variables
    virtual std::optional<VarAccess> replace(VarAccess&) {
        return std::nullopt;
    }
    // Expressions
    virtual std::optional<ptr<Expr>> replace(BExpr&) { return std::nullopt; }
    virtual std::optional<ptr<Expr>> replace(UExpr&) { return std::nullopt; }
    virtual std::optional<ptr<Expr>> replace(PiExpr&) { return std::nullopt; }
    virtual std::optional<ptr<Expr>> replace(IntExpr&) { return std::nullopt; }
    virtual std::optional<ptr<Expr>> replace(RealExpr&) { return std::nullopt; }
    virtual std::optional<ptr<Expr>> replace(VarExpr&) { return std::nullopt; }
    // Statements
    virtual std::optional<std::list<ptr<Stmt>>> replace(MeasureStmt&) {
        return std::nullopt;
    }
    virtual std::optional<std::list<ptr<Stmt>>> replace(ResetStmt&) {
        return std::nullopt;
    }
    virtual std::optional<std::list<ptr<Stmt>>> replace(IfStmt&) {
        return std::nullopt;
    }
    // Gates
    virtual std::optional<std::list<ptr<Gate>>> replace(UGate&) {
        return std::nullopt;
    }
    virtual std::optional<std::list<ptr<Gate>>> replace(CNOTGate&) {
        return std::nullopt;
    }
    virtual std::optional<std::list<ptr<Gate>>> replace(BarrierGate&) {
        return std::nullopt;
    }
    virtual std::optional<std::list<ptr<Gate>>> replace(DeclaredGate&) {
        return std::nullopt;
    }
    // Declarations
    virtual std::optional<std::list<ptr<Stmt>>> replace(GateDecl&) {
        return std::nullopt;
    }
    virtual std::optional<std::list<ptr<Stmt>>> replace(OracleDecl&) {
        return std::nullopt;
    }
    virtual std::optional<std::list<ptr<Stmt>>> replace(RegisterDecl&) {
        return std::nullopt;
    }
    virtual std::optional<std::list<ptr<Gate>>> replace(AncillaDecl&) {
        return std::nullopt;
    }

    /* Visitor overrides */
    void visit(VarAccess& var) override { replacement_var_ = replace(var); }

    void visit(BExpr& expr) override {
        expr.lexp().accept(*this);
        if (replacement_expr_) {
            expr.set_lexp(std::move(*replacement_expr_));
            replacement_expr_ = std::nullopt;
        }

        expr.rexp().accept(*this);
        if (replacement_expr_) {
            expr.set_rexp(std::move(*replacement_expr_));
            replacement_expr_ = std::nullopt;
        }

        replacement_expr_ = replace(expr);
    }

    void visit(UExpr& expr) override {
        expr.subexp().accept(*this);
        if (replacement_expr_) {
            expr.set_subexp(std::move(*replacement_expr_));
            replacement_expr_ = std::nullopt;
        }

        replacement_expr_ = replace(expr);
    }

    void visit(PiExpr& expr) override { replacement_expr_ = replace(expr); }
    void visit(IntExpr& expr) override { replacement_expr_ = replace(expr); }
    void visit(RealExpr& expr) override { replacement_expr_ = replace(expr); }
    void visit(VarExpr& expr) override { replacement_expr_ = replace(expr); }

    void visit(MeasureStmt& stmt) override {
        stmt.q_arg().accept(*this);
        if (replacement_var_) {
            stmt.set_qarg(std::move(*replacement_var_));
            replacement_var_ = std::nullopt;
        }

        stmt.c_arg().accept(*this);
        if (replacement_var_) {
            stmt.set_carg(std::move(*replacement_var_));
            replacement_var_ = std::nullopt;
        }

        replacement_stmts_ = replace(stmt);
    }

    void visit(ResetStmt& stmt) override {
        stmt.arg().accept(*this);
        if (replacement_var_) {
            stmt.set_arg(std::move(*replacement_var_));
            replacement_var_ = std::nullopt;
        }

        replacement_stmts_ = replace(stmt);
    }

    // Vanilla QASM only allows a single statement in
    // the "then" branch, so we need to clone the statement
    // for each gate in the result
    void visit(IfStmt& stmt) override {
        stmt.then().accept(*this);
        if (replacement_stmts_) {
            std::optional<std::list<ptr<Stmt>>> ret = std::nullopt;
            for (auto& rep : *replacement_stmts_) {
                auto tmp = object::clone(stmt);
                tmp->set_then(std::move(rep));
                auto stmts = replace(*tmp);
                if (!ret) {
                    ret = std::move(stmts);
                } else if (stmts) {
                    ret->splice(ret->end(), *stmts);
                }
            }
            replacement_stmts_ = std::move(ret);
        } else if (replacement_gates_) {
            std::optional<std::list<ptr<Stmt>>> ret = std::nullopt;
            for (auto& rep : *replacement_gates_) {
                auto tmp = object::clone(stmt);
                tmp->set_then(std::move(rep));
                auto stmts = replace(*tmp);
                if (!ret) {
                    ret = std::move(stmts);
                } else if (stmts) {
                    ret->splice(ret->end(), *stmts);
                }
            }
            replacement_gates_ = std::nullopt;
            replacement_stmts_ = std::move(ret);
        } else {
            replacement_stmts_ = replace(stmt);
        }
    }

    void visit(UGate& gate) override {
        gate.theta().accept(*this);
        if (replacement_expr_) {
            gate.set_theta(std::move(*replacement_expr_));
            replacement_expr_ = std::nullopt;
        }

        gate.phi().accept(*this);
        if (replacement_expr_) {
            gate.set_phi(std::move(*replacement_expr_));
            replacement_expr_ = std::nullopt;
        }

        gate.lambda().accept(*this);
        if (replacement_expr_) {
            gate.set_lambda(std::move(*replacement_expr_));
            replacement_expr_ = std::nullopt;
        }

        gate.arg().accept(*this);
        if (replacement_var_) {
            gate.set_arg(std::move(*replacement_var_));
            replacement_var_ = std::nullopt;
        }

        replacement_gates_ = replace(gate);
    }

    void visit(CNOTGate& gate) override {
        gate.ctrl().accept(*this);
        if (replacement_var_) {
            gate.set_ctrl(std::move(*replacement_var_));
            replacement_var_ = std::nullopt;
        }

        gate.tgt().accept(*this);
        if (replacement_var_) {
            gate.set_tgt(std::move(*replacement_var_));
            replacement_var_ = std::nullopt;
        }

        replacement_gates_ = replace(gate);
    }

    void visit(BarrierGate& gate) override {
        for (int i = 0; i < gate.num_args(); i++) {
            gate.arg(i).accept(*this);
            if (replacement_var_) {
                gate.set_arg(i, std::move(*replacement_var_));
                replacement_var_ = std::nullopt;
            }
        }

        replacement_gates_ = replace(gate);
    }

    void visit(DeclaredGate& gate) override {
        for (int i = 0; i < gate.num_cargs(); i++) {
            gate.carg(i).accept(*this);
            if (replacement_expr_) {
                gate.set_carg(i, std::move(*replacement_expr_));
                replacement_expr_ = std::nullopt;
            }
        }

        for (int i = 0; i < gate.num_qargs(); i++) {
            gate.qarg(i).accept(*this);
            if (replacement_var_) {
                gate.set_qarg(i, std::move(*replacement_var_));
                replacement_var_ = std::nullopt;
            }
        }

        replacement_gates_ = replace(gate);
    }

    void visit(GateDecl& decl) override {
        auto it = decl.body().begin();

        while (it != decl.body().end()) {
            (**it).accept(*this);
            if (replacement_gates_) {
                it = decl.body().erase(it);
                decl.body().splice(it, std::move(*replacement_gates_));
                replacement_gates_ = std::nullopt;
            } else {
                ++it;
            }
        }

        replacement_stmts_ = replace(decl);
    }

    void visit(OracleDecl& decl) override {
        replacement_stmts_ = replace(decl);
    }
    void visit(RegisterDecl& decl) override {
        replacement_stmts_ = replace(decl);
    }

    void visit(AncillaDecl& decl) override {
        replacement_gates_ = replace(decl);
    }

    void visit(Program& prog) override {
        auto it = prog.body().begin();

        while (it != prog.end()) {
            (**it).accept(*this);
            if (replacement_stmts_) {
                it = prog.body().erase(it);
                prog.body().splice(it, std::move(*replacement_stmts_));
                replacement_stmts_ = std::nullopt;
            } else if (replacement_gates_) {
                it = prog.body().erase(it);
                for (auto ti = replacement_gates_->begin();
                     ti != replacement_gates_->end(); ti++) {
                    prog.body().emplace(it, std::move(*ti));
                }
                replacement_gates_ = std::nullopt;
            } else {
                ++it;
            }
        }
    }
};

/**
 * \class qasmtools::ast::GateReplacer
 * \brief Bulk gate replacement
 * \see qasmtools::ast::Replacer
 *
 * Implements bulk replacement of gates given by a hash map. Use the
 * functional interface qasmtools::ast::replace_gates rather than
 * the class.
 */
class GateReplacer final : public Replacer {
  public:
    GateReplacer(std::unordered_map<int, std::list<ptr<Gate>>>&& replacements)
        : replacements_(std::move(replacements)) {}

    std::optional<std::list<ptr<Gate>>> replace(UGate& g) {
        return replace_gate(g);
    }
    std::optional<std::list<ptr<Gate>>> replace(CNOTGate& g) {
        return replace_gate(g);
    }
    std::optional<std::list<ptr<Gate>>> replace(BarrierGate& g) {
        return replace_gate(g);
    }
    std::optional<std::list<ptr<Gate>>> replace(DeclaredGate& g) {
        return replace_gate(g);
    }

  private:
    std::unordered_map<int, std::list<ptr<Gate>>> replacements_;

    std::optional<std::list<ptr<Gate>>> replace_gate(Gate& gate) {
        auto it = replacements_.find(gate.uid());

        if (it != replacements_.end()) {
            return std::move(it->second);
        } else {
            return std::nullopt;
        }
    }
};

/**
 * \brief Replaces the specified gates within an AST
 *
 * Used to perform a list of gate replacements in one traversal. All keys in the
 * hash map should refer to gates (U, CNOT, barrier or a declared gate).
 * For replacement of other types of nodes, use the Replacer class.
 *
 * \param node Reference to the root of the AST in which replacement will take
 * place
 * \param replacements Hash map from gate UID's to a list of gates which
 * should replace it
 */
inline void
replace_gates(ASTNode& node,
              std::unordered_map<int, std::list<ptr<Gate>>>&& replacements) {
    GateReplacer replacer(std::move(replacements));
    node.accept(replacer);
}

} /* namespace ast */
} /* namespace qasmtools */

#endif /* QASMTOOLS_AST_REPLACER_HPP_ */
