/*
 * This file is part of qasmtools.
 *
 * Copyright (c) 2019 - 2022 softwareQ Inc. All rights reserved.
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
 * \file qasmtools/tools/ast_printer.hpp
 * \brief Direct AST printer for debugging
 */

#pragma once

#include "../ast/ast.hpp"

namespace qasmtools {
namespace tools {

class ASTPrinter final : public ast::Visitor {
    std::ostream& os_;
    std::string prefix_;

  public:
    ASTPrinter(std::ostream& os) : os_(os), prefix_("") {}

    // Variables
    void visit(ast::VarAccess& ap) {
        os_ << prefix_ << "|- Var(" << ap.var();
        if (ap.offset())
            os_ << "[" << *ap.offset() << "]";
        os_ << ")\n";
    }

    // Expressions
    void visit(ast::BExpr& expr) {
        os_ << prefix_ << "|- BExpr(" << expr.op() << ")\n";

        prefix_ += "  ";
        expr.lexp().accept(*this);
        expr.rexp().accept(*this);
        prefix_.resize(prefix_.length() - 2);
    }

    void visit(ast::UExpr& expr) {
        os_ << prefix_ << "|- UExpr(" << expr.op() << ")\n";

        prefix_ += "  ";
        expr.subexp().accept(*this);
        prefix_.resize(prefix_.length() - 2);
    }

    void visit(ast::PiExpr&) { os_ << prefix_ << "|- Pi\n"; }

    void visit(ast::IntExpr& expr) {
        os_ << prefix_ << "|- Int(" << expr.value() << ")\n";
    }

    void visit(ast::RealExpr& expr) {
        os_ << prefix_ << "|- Real(" << expr.value() << ")\n";
    }

    void visit(ast::VarExpr& expr) {
        os_ << prefix_ << "|- Var(" << expr.var();
        os_ << ")\n";
    }

    // Statements
    void visit(ast::MeasureStmt& stmt) {
        os_ << prefix_ << "|- Measure\n";

        prefix_ += "  ";
        stmt.q_arg().accept(*this);
        stmt.c_arg().accept(*this);
        prefix_.resize(prefix_.length() - 2);
    }

    void visit(ast::ResetStmt& stmt) {
        os_ << prefix_ << "|- Reset\n";

        prefix_ += "  ";
        stmt.arg().accept(*this);
        prefix_.resize(prefix_.length() - 2);
    }

    void visit(ast::IfStmt& stmt) {
        os_ << prefix_ << "|- If(" << stmt.var() << "==" << stmt.cond()
            << ")\n";

        prefix_ += "  ";
        stmt.then().accept(*this);
        prefix_.resize(prefix_.length() - 2);
    }

    // Gates
    void visit(ast::UGate& gate) {
        os_ << prefix_ << "|- UGate\n";

        prefix_ += "  ";
        gate.theta().accept(*this);
        gate.phi().accept(*this);
        gate.lambda().accept(*this);
        gate.arg().accept(*this);
        prefix_.resize(prefix_.length() - 2);
    }

    void visit(ast::CNOTGate& gate) {
        os_ << prefix_ << "|- CXGate\n";

        prefix_ += "  ";
        gate.ctrl().accept(*this);
        gate.tgt().accept(*this);
        prefix_.resize(prefix_.length() - 2);
    }

    void visit(ast::BarrierGate& gate) {
        os_ << prefix_ << "|- Barrier\n";

        prefix_ += "  ";
        gate.foreach_arg([this](auto& arg) { arg.accept(*this); });
        prefix_.resize(prefix_.length() - 2);
    }

    void visit(ast::DeclaredGate& gate) {
        os_ << prefix_ << "|- Declared(" << gate.name() << ")\n";

        prefix_ += "  ";
        gate.foreach_qarg([this](auto& arg) { arg.accept(*this); });
        gate.foreach_carg([this](auto& arg) { arg.accept(*this); });
        prefix_.resize(prefix_.length() - 2);
    }

    // Declarations
    void visit(ast::GateDecl& decl) {
        os_ << prefix_ << "|- Gate Decl(" << decl.id() << "(";
        for (auto& param : decl.c_params())
            os_ << param << ",";
        os_ << ")[";
        for (auto& param : decl.q_params())
            os_ << param << ",";
        os_ << "]";
        if (decl.is_opaque())
            os_ << ", opaque";
        os_ << ")\n";

        prefix_ += "  ";
        decl.foreach_stmt([this](auto& stmt) { stmt.accept(*this); });
        prefix_.resize(prefix_.length() - 2);
    }

    void visit(ast::OracleDecl& decl) {
        os_ << prefix_ << "|- Oracle Decl(" << decl.id() << "[";
        for (auto& param : decl.params())
            os_ << param << ",";
        os_ << "] = " << decl.fname() << ")\n";
    }

    void visit(ast::RegisterDecl& decl) {
        os_ << prefix_ << "|- Register Decl(" << decl.id() << "[" << decl.size()
            << "]";
        if (decl.is_quantum())
            os_ << ", quantum";
        os_ << ")\n";
    }

    void visit(ast::AncillaDecl& decl) {
        os_ << prefix_ << "|- Ancilla Decl(" << decl.id() << "[" << decl.size()
            << "]";
        if (decl.is_dirty())
            os_ << ", dirty";
        os_ << ")\n";
    }

    // Program
    void visit(ast::Program& prog) {
        os_ << prefix_ << "|- Program\n";

        prefix_ += "  ";
        prog.foreach_stmt([this](auto& stmt) { stmt.accept(*this); });
        prefix_.resize(prefix_.length() - 2);
    }
};

void print_tree(ast::ASTNode& node, std::ostream& os = std::cout) {
    ASTPrinter printer(os);
    node.accept(printer);
}

} // namespace tools
} // namespace qasmtools
