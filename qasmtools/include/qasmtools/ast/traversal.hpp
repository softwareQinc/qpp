/*
 * This file is part of qasmtools.
 *
 * Copyright (c) 2019 - 2023 softwareQ Inc. All rights reserved.
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
 * \file qasmtools/ast/traversal.hpp
 * \brief Node traversal for syntax trees
 */

#pragma once

#include "program.hpp"
#include "visitor.hpp"

namespace qasmtools {
namespace ast {

/**
 * \class qasmtools::ast::Traverse
 * \brief Generic complete traversal of ASTs
 * \see qasmtools::ast::Visitor
 *
 * Implements a generic, pass-through traversal of the entire
 * AST. Standard usage is to derive from this class and override only
 * the nodes desired.
 *
 * Note that overriding a node kills traversal to
 * children of that node. This can be useful for cutting off traversal
 * to certain subtrees. Traversal logic can be accessed through the
 * parent class, e.g. by calling Traverse::visit. In this way, the
 * Traverse class can be used to implement post-order, pre-order, or
 * mixed pre/post algorithms by directly calling the traversal logic.
 */
class Traverse : public Visitor {
  public:
    void visit(VarAccess& var) override {}
    void visit(BExpr& expr) override {
        expr.lexp().accept(*this);
        expr.rexp().accept(*this);
    }
    void visit(UExpr& expr) override { expr.subexp().accept(*this); }
    void visit(PiExpr& expr) override {}
    void visit(IntExpr& expr) override {}
    void visit(RealExpr& expr) override {}
    void visit(VarExpr& expr) override {}
    void visit(MeasureStmt& stmt) override {
        stmt.q_arg().accept(*this);
        stmt.c_arg().accept(*this);
    }
    void visit(ResetStmt& stmt) override { stmt.arg().accept(*this); }
    void visit(IfStmt& stmt) override { stmt.then().accept(*this); }
    void visit(UGate& gate) override {
        gate.theta().accept(*this);
        gate.phi().accept(*this);
        gate.lambda().accept(*this);
        gate.arg().accept(*this);
    }
    void visit(CNOTGate& gate) override {
        gate.ctrl().accept(*this);
        gate.tgt().accept(*this);
    }
    void visit(BarrierGate& gate) override {
        for (int i = 0; i < gate.num_args(); i++)
            gate.arg(i).accept(*this);
    }
    void visit(DeclaredGate& gate) override {
        for (int i = 0; i < gate.num_cargs(); i++)
            gate.carg(i).accept(*this);
        for (int i = 0; i < gate.num_qargs(); i++)
            gate.qarg(i).accept(*this);
    }

    void visit(GateDecl& decl) override {
        for (auto it = decl.begin(); it != decl.end(); it++)
            (**it).accept(*this);
    }

    void visit(OracleDecl& decl) override {}
    void visit(RegisterDecl& decl) override {}
    void visit(AncillaDecl& decl) override {}
    void visit(Program& prog) override {
        for (auto it = prog.begin(); it != prog.end(); it++)
            (**it).accept(*this);
    }
};

} // namespace ast
} // namespace qasmtools
