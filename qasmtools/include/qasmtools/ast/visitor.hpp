/*
 * This file is part of qasmtools.
 *
 * Copyright (c) 2019 - 2025 softwareQ Inc. All rights reserved.
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
 * \file qasmtools/ast/visitor.hpp
 * \brief Visitor interface for syntax trees
 */

#ifndef QASMTOOLS_AST_VISITOR_HPP_
#define QASMTOOLS_AST_VISITOR_HPP_

namespace qasmtools {
namespace ast {

/* Forward declarations */
class VarAccess;
class BExpr;
class UExpr;
class PiExpr;
class IntExpr;
class RealExpr;
class VarExpr;
class MeasureStmt;
class ResetStmt;
class IfStmt;
class UGate;
class CNOTGate;
class BarrierGate;
class DeclaredGate;
class GateDecl;
class OracleDecl;
class RegisterDecl;
class AncillaDecl;
class Program;

/**
 * \class qasmtools::ast::Visitor
 * \brief Base visitor interface
 *
 * Classic visitor via (virtual) double dispatch. Standard usage is to
 * derive from this class and provide implementations of visit for **every**
 * node type.
 *
 * Traversal to sub-nodes is handled by the particular visitor, not the
 * node class. For a visitor that automatically handles traversal and also
 * allows picking and choosing the particular visit overloads, see
 * qasmtools::ast::Traverse.
 */
class Visitor {
  public:
    // Variables
    virtual void visit(VarAccess&) = 0;
    // Expressions
    virtual void visit(BExpr&) = 0;
    virtual void visit(UExpr&) = 0;
    virtual void visit(PiExpr&) = 0;
    virtual void visit(IntExpr&) = 0;
    virtual void visit(RealExpr&) = 0;
    virtual void visit(VarExpr&) = 0;
    // Statements
    virtual void visit(MeasureStmt&) = 0;
    virtual void visit(ResetStmt&) = 0;
    virtual void visit(IfStmt&) = 0;
    // Gates
    virtual void visit(UGate&) = 0;
    virtual void visit(CNOTGate&) = 0;
    virtual void visit(BarrierGate&) = 0;
    virtual void visit(DeclaredGate&) = 0;
    // Declarations
    virtual void visit(GateDecl&) = 0;
    virtual void visit(OracleDecl&) = 0;
    virtual void visit(RegisterDecl&) = 0;
    virtual void visit(AncillaDecl&) = 0;
    // Program
    virtual void visit(Program&) = 0;
    // Destructor
    virtual ~Visitor() = default;
};

} /* namespace ast */
} /* namespace qasmtools */

#endif /* QASMTOOLS_AST_VISITOR_HPP_ */
