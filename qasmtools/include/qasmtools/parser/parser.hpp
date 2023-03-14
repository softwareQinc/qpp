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
 *
 * Adapted from Bruno Schmitt's Tweedledee library
 */

/**
 * \file qasmtools/parser/parser.hpp
 * \brief OpenQASM parsing
 */

#ifndef QASMTOOLS_PARSER_PARSER_HPP_
#define QASMTOOLS_PARSER_PARSER_HPP_

#include "../ast/ast.hpp"
#include "preprocessor.hpp"

#include <list>

namespace qasmtools {
namespace parser {

/**
 * \class qasmtools::parser::ParseError
 * \brief Exception class for parse errors
 */
class ParseError : public std::exception {
  public:
    ParseError() noexcept = default;
    ~ParseError() = default;
    const char* what() const noexcept { return "Parse error(s)"; }
};

/**
 * \class qasmtools::parser::Parser
 * \brief OpenQASM parser class
 * \see qasmtools::parser::Preprocessor
 */
class Parser {
    Preprocessor& pp_lexer_; ///< preprocessed, tokenized input stream

    bool error_ = false;          ///< whether a parse error has occured
    bool supress_errors_ = false; ///< whether to supress errors
    Token current_token_;         ///< current token
    int bits_ = 0;                ///< number of bits
    int qubits_ = 0;              ///< number of qubits

  public:
    /**
     * \brief Constructs a parser for a given preprocessor stream
     *
     * \param pp_lexer The preprocessed lexer stream to be parsed
     */
    Parser(Preprocessor& pp_lexer)
        : pp_lexer_(pp_lexer), current_token_(pp_lexer_.next_token()) {}

    /**
     * \brief Parses the tokenized stream as a QCircuit object
     *
     * \return A unique pointer to a QCircuit object
     */
    ast::ptr<ast::Program> parse(bool check = true) {
        // Parse the program
        auto result = parse_program();
        if (error_)
            throw ParseError();

        // Perform semantic analysis before returning
        if (check)
            ast::check_source(*result);

        return result;
    }

  private:
    /**
     * \brief Consume a token and retrieve the next one
     *
     * \param reset Whether to unsupress errors (optional, default is false)
     */
    void consume_token(bool reset = false) {
        current_token_ = pp_lexer_.next_token();
        if (reset)
            supress_errors_ = false;
    }

    /**
     * \brief Consume a particular type of token
     *
     * Checks that the current token is of type given by the argument
     * and consumes the token, setting the error flag if it is not
     * of the correct type
     *
     * \param expected The type of token to be consumed
     * \return The current token
     */
    Token expect_and_consume_token(Token::Kind expected) {
        if (current_token_.is_not(expected)) {
            error_ = true;
            if (!supress_errors_) {
                std::cerr << current_token_.position();
                std::cerr << ": expected " << expected;
                std::cerr << " but got " << current_token_.kind() << "\n";
                ;
                supress_errors_ = true;
            }
            return current_token_;
        }

        auto return_token = current_token_;
        consume_token();
        return return_token;
    }

    /**
     * \brief Consume tokens until a particular token is found
     *
     * Repeatedly skips tokens, setting the error flag if necessary,
     * until the given token or eof is found
     *
     * \param expected The type of token to be consumed
     * \return The next expected token, or eof
     */
    Token consume_until(Token::Kind expected) {
        while (current_token_.is_not(expected) &&
               current_token_.is_not(Token::Kind::eof)) {
            error_ = true;
            if (!supress_errors_) {
                std::cerr << current_token_.position();
                std::cerr << ": expected " << expected;
                std::cerr << " but got " << current_token_.kind() << "\n";
                ;
                supress_errors_ = true;
            }
            consume_token();
        }

        auto return_token = current_token_;
        consume_token(true);
        return return_token;
    }

    /**
     * \brief Try to consume a particular type of token
     *
     * Attempts to parse a particular type of token, only consuming
     * the token if it is of the correct type
     *
     * \param expected The type of token to be consumed
     * \return True if and only if a token of type expected was consumed
     */
    bool try_and_consume_token(Token::Kind expected) {
        if (current_token_.is_not(expected)) {
            return false;
        }
        consume_token();
        return true;
    }

    /**
     * \brief Parse an OpenQASM 2.0 program
     *
     * <mainprogram> = OPENQASM <real> ; <program>
     * <program>     = <statement> | <program> <statement>
     * <statement>   = <regdecl>
     *               | <gatedecl> <goplist> }
     *               | <gatedecl> }
     *               | <opaquedecl> ;
     *               | <qop>
     *               | if ( <id> == <nninteger> ) <qop>
     *               | barrier <anylist> ;
     *
     * \return A QASM AST object
     */
    ast::ptr<ast::Program> parse_program() {
        auto pos = current_token_.position();
        std::list<ast::ptr<ast::Stmt>> ret;

        // The first (non-comment) line of an Open QASM program must be
        // OPENQASM M.m; indicating a major version M and minor version m.
        parse_header();

        while (!current_token_.is(Token::Kind::eof)) {
            switch (current_token_.kind()) {
                    // Parse declarations (<decl>)
                case Token::Kind::kw_creg:
                    ret.emplace_back(parse_reg_decl(false));
                    break;
                case Token::Kind::kw_qreg:
                    ret.emplace_back(parse_reg_decl(true));
                    break;

                case Token::Kind::kw_gate:
                    ret.emplace_back(parse_gate_decl());
                    break;

                case Token::Kind::kw_opaque:
                    ret.emplace_back(parse_opaque_decl());
                    break;
                case Token::Kind::kw_oracle:
                    ret.emplace_back(parse_oracle_decl());
                    break;

                    // Parse quantum operations (<qop>)
                case Token::Kind::identifier:
                case Token::Kind::kw_cx:
                case Token::Kind::kw_measure:
                case Token::Kind::kw_reset:
                case Token::Kind::kw_u:
                    ret.emplace_back(parse_qop());
                    break;

                case Token::Kind::kw_barrier:
                    ret.emplace_back(parse_barrier());
                    break;

                case Token::Kind::kw_if:
                    ret.emplace_back(parse_if());
                    break;

                default:
                    error_ = true;
                    if (!supress_errors_) {
                        std::cerr << current_token_.position();
                        std::cerr
                            << ": expected a global declaration or statement";
                        std::cerr << " but got " << current_token_.kind()
                                  << "\n";
                        ;
                        supress_errors_ = true;
                    }

                    consume_until(Token::Kind::semicolon);
                    break;
            }
        }

        return ast::Program::create(pos, pp_lexer_.includes_stdlib(),
                                    std::move(ret), bits_, qubits_);
    }

    /**
     * \brief Parse an OpenQASM 2.0 header
     *
     * OPENQASM <nninteger> ;
     */
    void parse_header() {
        consume_token();
        expect_and_consume_token(Token::Kind::kw_openqasm);
        expect_and_consume_token(Token::Kind::real);
        consume_until(Token::Kind::semicolon);
    }

    /**
     * \brief Parse a register declaration
     * \note <regdecl>
     *
     * <regdecl> = qreg <id> [ <nninteger> ] ;
     *           | creg <id> [ <nninteger> ] ;
     *
     * \return Unique pointer to a statement object
     */
    ast::ptr<ast::Stmt> parse_reg_decl(bool quantum) {
        auto pos = current_token_.position();

        consume_token();
        auto id = expect_and_consume_token(Token::Kind::identifier);
        expect_and_consume_token(Token::Kind::l_square);
        auto size = expect_and_consume_token(Token::Kind::nninteger);
        expect_and_consume_token(Token::Kind::r_square);
        consume_until(Token::Kind::semicolon);

        quantum ? qubits_ += size.as_int() : bits_ += size.as_int();
        return ast::RegisterDecl::create(pos, id.as_string(), quantum,
                                         size.as_int());
    }

    /**
     * \brief Parse an ancilla declaration
     * \note <ancdecl>
     *
     * <ancdecl> = ancilla <id> [ <nninteger> ] ;
     *           | dirty ancilla <id> [ <nninteger> ] ;
     *
     * \return Unique pointer to a (gate) statement object
     */
    ast::ptr<ast::Gate> parse_ancilla_decl() {
        bool dirty = false;
        auto pos = current_token_.position();

        if (try_and_consume_token(Token::Kind::kw_dirty)) {
            dirty = true;
        }

        expect_and_consume_token(Token::Kind::kw_ancilla);
        auto id = expect_and_consume_token(Token::Kind::identifier);
        expect_and_consume_token(Token::Kind::l_square);
        auto size = expect_and_consume_token(Token::Kind::nninteger);
        expect_and_consume_token(Token::Kind::r_square);
        consume_until(Token::Kind::semicolon);

        return ast::AncillaDecl::create(pos, id.as_string(), dirty,
                                        size.as_int());
    }

    /**
     * \brief Parse a gate declaration
     * \note <gatedecl>
     *
     * <gatedecl> = gate <id> <idlist> {
     *            | gate <id> ( ) <idlist> {
     *            | gate <id> ( <idlist> ) <idlist> {
     *
     * \return Unique pointer to a statement object
     */
    ast::ptr<ast::Stmt> parse_gate_decl() {
        auto pos = current_token_.position();

        consume_token();
        auto id = expect_and_consume_token(Token::Kind::identifier);

        std::vector<ast::symbol> c_params;
        if (try_and_consume_token(Token::Kind::l_paren)) {
            c_params = parse_idlist();
            expect_and_consume_token(Token::Kind::r_paren);
        }

        auto q_params = parse_idlist();

        expect_and_consume_token(Token::Kind::l_brace);
        std::list<ast::ptr<ast::Gate>> body;
        if (!try_and_consume_token(Token::Kind::r_brace)) {
            body = parse_goplist();
            expect_and_consume_token(Token::Kind::r_brace);
        }

        return ast::GateDecl::create(pos, id.as_string(), false, c_params,
                                     q_params, std::move(body));
    }

    /**
     * \brief Parse an opaque gate declaration
     * \note <opaquedecl>
     *
     * <opaquedecl> = opaque <id> <idlist> ;
     *              | opaque <id> ( ) <idlist> ;
     *              | opaque <id> ( <idlist> ) <idlist> ;
     *
     * \return Unique pointer to a statement object
     */
    ast::ptr<ast::Stmt> parse_opaque_decl() {
        auto pos = current_token_.position();

        consume_token();
        auto id = expect_and_consume_token(Token::Kind::identifier);

        std::vector<ast::symbol> c_params;
        if (try_and_consume_token(Token::Kind::l_paren)) {
            c_params = parse_idlist();
            expect_and_consume_token(Token::Kind::r_paren);
        }

        auto q_params = parse_idlist();

        consume_until(Token::Kind::semicolon);
        return ast::GateDecl::create(pos, id.as_string(), true, c_params,
                                     q_params, {});
    }

    /**
     * \brief Parse an oracle declaration
     * \note <oracledecl>
     *
     * <oracledecl> = oracle <id> <idlist> { <string> }
     *
     * \return Unique pointer to a statement object
     */
    ast::ptr<ast::Stmt> parse_oracle_decl() {
        auto pos = current_token_.position();

        consume_token();
        auto id = expect_and_consume_token(Token::Kind::identifier);
        auto params = parse_idlist();
        expect_and_consume_token(Token::Kind::l_brace);
        auto file = expect_and_consume_token(Token::Kind::string);
        expect_and_consume_token(Token::Kind::r_brace);

        return ast::OracleDecl::create(pos, id.as_string(), params,
                                       file.as_string());
    }

    /**
     * \brief Parse a list of gate operations
     * \note <goplist>
     *
     * <goplist> = <gop>
     *           | <goplist> <gop>
     *
     * \return Vector of gate objects
     */
    std::list<ast::ptr<ast::Gate>> parse_goplist() {
        std::list<ast::ptr<ast::Gate>> ret;
        bool finished = false;

        while (!finished) {
            switch (current_token_.kind()) {
                case Token::Kind::kw_dirty:
                case Token::Kind::kw_ancilla:
                case Token::Kind::kw_cx:
                case Token::Kind::kw_u:
                case Token::Kind::identifier:
                case Token::Kind::kw_barrier:
                    ret.emplace_back(parse_gop());
                    break;

                default:
                    finished = true;
                    break;
            }
        }

        return ret;
    }

    /**
     * \brief Parse a quantum operation
     * \note <qop>
     *
     * <qop> = <gop>
     *       | measure <argument> -> <argument> ;
     *       | reset <argument> ;
     *
     * \return Unique pointer to a statement object
     */
    ast::ptr<ast::Stmt> parse_qop() {
        switch (current_token_.kind()) {
            case Token::Kind::kw_cx:
            case Token::Kind::kw_u:
            case Token::Kind::identifier:
                return parse_gop();

            case Token::Kind::kw_measure:
                return parse_measure();

            case Token::Kind::kw_reset:
                return parse_reset();

            default:
                error_ = true;
                if (!supress_errors_) {
                    std::cerr << current_token_.position();
                    std::cerr << ": expected a quantum operation, but got ";
                    std::cerr << current_token_.kind() << "\n";
                    ;
                    supress_errors_ = true;
                }

                return nullptr;
        }
    }

    /**
     * \brief Parse a gate operation
     * \note <gop>
     *
     * <gop> = U ( <explist> ) <argument> ;
     *       | CX <argument> , <argument> ;
     *       | <id> <anylist> ;
     *       | <id> ( ) <anylist> ;
     *       | <id> ( <explist> ) <anylist> ;
     *       | barrier <idlist> ;
     *
     * \return Unique pointer to a gate object
     */
    ast::ptr<ast::Gate> parse_gop() {
        switch (current_token_.kind()) {
            case Token::Kind::kw_dirty:
            case Token::Kind::kw_ancilla:
                return parse_ancilla_decl();

            case Token::Kind::kw_u:
                return parse_unitary();

            case Token::Kind::kw_cx:
                return parse_cnot();

            case Token::Kind::identifier:
                return parse_gate_statement();

            case Token::Kind::kw_barrier:
                return parse_barrier();

            default:
                error_ = true;
                if (!supress_errors_) {
                    std::cerr << current_token_.position();
                    std::cerr << ": expected a gate operation but got ";
                    std::cerr << current_token_.kind() << "\n";
                    ;
                    supress_errors_ = true;
                }

                return nullptr;
        }
    }

    /**
     * \brief Parse a list of identifiers
     * \note <idlist>
     *
     * <idlist> = <id>
     *          | <idlist> , <id>
     *
     * \return Vector of identifiers
     */
    std::vector<ast::symbol> parse_idlist() {
        std::vector<ast::symbol> ret;

        // doesn't accept empty lists
        while (true) {
            auto id = expect_and_consume_token(Token::Kind::identifier);
            ret.push_back(id.as_string());
            if (!try_and_consume_token(Token::Kind::comma)) {
                break;
            }
        }

        return ret;
    }

    /**
     * \brief Parse a variable expression
     * \note <var>
     *
     * <var> = <id>
     *       | <id> [ <nninteger> ]
     *
     * \return Variable object
     */
    ast::VarAccess parse_argument() {
        auto pos = current_token_.position();

        auto id = expect_and_consume_token(Token::Kind::identifier);
        if (try_and_consume_token(Token::Kind::l_square)) {
            auto offset = expect_and_consume_token(Token::Kind::nninteger);
            expect_and_consume_token(Token::Kind::r_square);

            return ast::VarAccess(pos, id.as_string(), offset.as_int());
        }

        return ast::VarAccess(pos, id.as_string());
    }

    /**
     * \brief Parse a list of variable expressions
     * \note <varlist>
     *
     * <varlist> = <var>
     *           | <varlist> , <var>
     *
     * \return Vector of variable objects
     */
    std::vector<ast::VarAccess> parse_anylist() {
        std::vector<ast::VarAccess> ret;

        // doesn't accept empty lists
        while (true) {
            ret.emplace_back(parse_argument());
            if (!try_and_consume_token(Token::Kind::comma)) {
                break;
            }
        }

        return ret;
    }

    /**
     * \brief Parse a list of expressions
     * \note <explist>
     *
     * <explist> = <exp>
     *           | <explist> , <exp>
     *
     * \return Vector of expressions
     */
    std::vector<ast::ptr<ast::Expr>> parse_explist() {
        // doesn't accept empty lists
        std::vector<ast::ptr<ast::Expr>> ret;

        while (true) {
            ret.emplace_back(parse_exp());
            if (!try_and_consume_token(Token::Kind::comma)) {
                break;
            }
        }

        return ret;
    }

    /**
     * \brief Parse an expression
     * \note <exp>
     *
     * <exp> = <real> | <nninteger> | pi | <id>
     *       | <exp> + <exp> | <exp> - <exp>
     *       | <exp> * <exp> | <exp> / <exp>
     *       | - <exp>
     *       | <exp> ^ <exp>
     *       | ( <exp> )
     *       | <unaryop> ( <exp> )
     *
     * \return Unique pointer to an expression object
     */
    ast::ptr<ast::Expr> parse_exp(int min_precedence = 1) {
        auto pos = current_token_.position();

        auto lexp = parse_atom();
        while (1) {
            auto next_min_precedence = min_precedence;

            switch (current_token_.kind()) {
                case Token::Kind::plus:
                    if (min_precedence > 1) {
                        return lexp;
                    }
                    next_min_precedence = 2;
                    break;

                case Token::Kind::minus:
                    if (min_precedence > 1) {
                        return lexp;
                    }
                    next_min_precedence = 2;
                    break;

                case Token::Kind::star:
                    if (min_precedence > 2) {
                        return lexp;
                    }
                    next_min_precedence = 3;
                    break;

                case Token::Kind::slash:
                    if (min_precedence > 2) {
                        return lexp;
                    }
                    next_min_precedence = 3;
                    break;

                case Token::Kind::caret:
                    if (min_precedence > 3) {
                        return lexp;
                    }
                    next_min_precedence = 4;
                    break;

                default:
                    return lexp;
            }

            auto bop = parse_binaryop();
            auto rexp = parse_exp(next_min_precedence);
            lexp =
                ast::BExpr::create(pos, std::move(lexp), bop, std::move(rexp));
        }

        return lexp;
    }

    /**
     * \brief Parse an atomic expression
     * \see qpp:qasmtools::Parser::parse_exp()
     *
     * \return Unique pointer to an expression object
     */
    ast::ptr<ast::Expr> parse_atom() {
        auto pos = current_token_.position();

        switch (current_token_.kind()) {
            case Token::Kind::l_paren: {
                consume_token();
                auto exp = parse_exp();
                expect_and_consume_token(Token::Kind::r_paren);
                return exp;
            }

            case Token::Kind::minus: {
                consume_token();
                auto exp = parse_atom();
                return ast::UExpr::create(pos, ast::UnaryOp::Neg,
                                          std::move(exp));
            }

            case Token::Kind::plus: {
                consume_token();
                return parse_atom();
            }

            case Token::Kind::identifier: {
                auto id = current_token_;
                consume_token();
                return ast::VarExpr::create(pos, id.as_string());
            }
            case Token::Kind::nninteger: {
                auto integer = current_token_;
                consume_token();
                return ast::IntExpr::create(pos, integer.as_int());
            }
            case Token::Kind::kw_pi:
                consume_token();
                return ast::PiExpr::create(pos);
            case Token::Kind::real: {
                auto real = current_token_;
                consume_token();
                return ast::RealExpr::create(pos, real.as_real());
            }

            case Token::Kind::kw_sin:
            case Token::Kind::kw_cos:
            case Token::Kind::kw_tan:
            case Token::Kind::kw_exp:
            case Token::Kind::kw_ln:
            case Token::Kind::kw_sqrt: {
                auto op = parse_unaryop();
                expect_and_consume_token(Token::Kind::l_paren);
                auto exp = parse_exp();
                expect_and_consume_token(Token::Kind::r_paren);
                return ast::UExpr::create(pos, op, std::move(exp));
            }

            default:
                error_ = true;
                if (!supress_errors_) {
                    std::cerr << current_token_.position();
                    std::cerr << ": expected an atomic expression but got ";
                    std::cerr << current_token_.kind() << "\n";
                    ;
                    supress_errors_ = true;
                }

                return nullptr;
        }
    }

    /**
     * \brief Parse a binary operator
     * \see qpp:qasmtools::Parser::parse_exp()
     *
     * \return Binary operator
     */
    ast::BinaryOp parse_binaryop() {
        switch (current_token_.kind()) {
            case Token::Kind::plus:
                consume_token();
                return ast::BinaryOp::Plus;
            case Token::Kind::minus:
                consume_token();
                return ast::BinaryOp::Minus;
            case Token::Kind::star:
                consume_token();
                return ast::BinaryOp::Times;
            case Token::Kind::slash:
                consume_token();
                return ast::BinaryOp::Divide;
            case Token::Kind::caret:
                consume_token();
                return ast::BinaryOp::Pow;
            default:
                error_ = true;
                if (!supress_errors_) {
                    std::cerr << current_token_.position();
                    std::cerr << ": expected a binary operator but got ";
                    std::cerr << current_token_.kind() << "\n";
                    ;
                    supress_errors_ = true;
                }

                return ast::BinaryOp::Plus;
        }
    }

    /**
     * \brief Parse a unary operator
     * \see qpp:qasmtools::Parser::parse_exp()
     *
     * \return Unary operator
     */
    ast::UnaryOp parse_unaryop() {
        switch (current_token_.kind()) {
            case Token::Kind::kw_sin:
                consume_token();
                return ast::UnaryOp::Sin;
            case Token::Kind::kw_cos:
                consume_token();
                return ast::UnaryOp::Cos;
            case Token::Kind::kw_tan:
                consume_token();
                return ast::UnaryOp::Tan;
            case Token::Kind::kw_exp:
                consume_token();
                return ast::UnaryOp::Exp;
            case Token::Kind::kw_ln:
                consume_token();
                return ast::UnaryOp::Ln;
            case Token::Kind::kw_sqrt:
                consume_token();
                return ast::UnaryOp::Sqrt;
            default:
                error_ = true;
                if (!supress_errors_) {
                    std::cerr << current_token_.position();
                    std::cerr << ": expected a unary operator but got ";
                    std::cerr << current_token_.kind() << "\n";
                    ;
                    supress_errors_ = true;
                }

                return ast::UnaryOp::Neg;
        }
    }

    /**
     * \brief Parse a CNOT gate
     *
     * CX <var> , <var> ;
     *
     * \return Unique pointer to a gate object
     */
    ast::ptr<ast::Gate> parse_cnot() {
        auto pos = current_token_.position();

        expect_and_consume_token(Token::Kind::kw_cx);
        auto ctrl = parse_argument();
        expect_and_consume_token(Token::Kind::comma);
        auto tgt = parse_argument();
        consume_until(Token::Kind::semicolon);

        return ast::CNOTGate::create(pos, std::move(ctrl), std::move(tgt));
    }

    /**
     * \brief Parse a single qubit U gate
     *
     * U ( <exp> , <exp> , <exp> ) <var> ;
     *
     * \return Unique pointer to a gate object
     */
    ast::ptr<ast::Gate> parse_unitary() {
        auto pos = current_token_.position();

        expect_and_consume_token(Token::Kind::kw_u);
        expect_and_consume_token(Token::Kind::l_paren);
        auto theta = parse_exp();
        expect_and_consume_token(Token::Kind::comma);
        auto phi = parse_exp();
        expect_and_consume_token(Token::Kind::comma);
        auto lambda = parse_exp();
        expect_and_consume_token(Token::Kind::r_paren);
        auto arg = parse_argument();
        consume_until(Token::Kind::semicolon);

        return ast::UGate::create(pos, std::move(theta), std::move(phi),
                                  std::move(lambda), std::move(arg));
    }

    /**
     * \brief Parse a declared gate application
     *
     * <id> <varlist> ;
     * <id> ( <explist> ) <varlist> ;
     *
     * \return Unique pointer to a gate object
     */
    ast::ptr<ast::Gate> parse_gate_statement() {
        auto pos = current_token_.position();

        auto id = expect_and_consume_token(Token::Kind::identifier);
        std::vector<ast::ptr<ast::Expr>> c_args;
        if (try_and_consume_token(Token::Kind::l_paren)) {
            c_args = parse_explist();
            expect_and_consume_token(Token::Kind::r_paren);
        }

        auto q_args = parse_anylist();
        consume_until(Token::Kind::semicolon);

        return ast::DeclaredGate::create(pos, id.as_string(), std::move(c_args),
                                         std::move(q_args));
    }

    /**
     * \brief Parse a measurement
     *
     * measure <var> -> <var> ;
     *
     * \return Unique pointer to a statement object
     */
    ast::ptr<ast::Stmt> parse_measure() {
        auto pos = current_token_.position();

        expect_and_consume_token(Token::Kind::kw_measure);
        auto q_arg = parse_argument();
        expect_and_consume_token(Token::Kind::arrow);
        auto c_arg = parse_argument();
        consume_until(Token::Kind::semicolon);

        return ast::MeasureStmt::create(pos, std::move(q_arg),
                                        std::move(c_arg));
    }

    /**
     * \brief Parse a reset statement
     *
     * reset <var> ;
     *
     * \return Unique pointer to a statement object
     */
    ast::ptr<ast::Stmt> parse_reset() {
        auto pos = current_token_.position();

        expect_and_consume_token(Token::Kind::kw_reset);
        auto arg = parse_argument();
        consume_until(Token::Kind::semicolon);

        return ast::ResetStmt::create(pos, std::move(arg));
    }

    /**
     * \brief Parse a barrier statement
     *
     * barrier <varlist> ;
     *
     * \return Unique pointer to a gate object
     */
    ast::ptr<ast::Gate> parse_barrier() {
        auto pos = current_token_.position();

        expect_and_consume_token(Token::Kind::kw_barrier);
        auto args = parse_anylist();
        consume_until(Token::Kind::semicolon);

        return ast::BarrierGate::create(pos, std::move(args));
    }

    /**
     * \brief Parse an if statement
     *
     * if ( <id> == <nninteger> ) <qop>
     *
     * \return Unique pointer to a statement object
     */
    ast::ptr<ast::Stmt> parse_if() {
        auto pos = current_token_.position();

        expect_and_consume_token(Token::Kind::kw_if);
        expect_and_consume_token(Token::Kind::l_paren);
        auto id = expect_and_consume_token(Token::Kind::identifier);
        expect_and_consume_token(Token::Kind::equalequal);
        auto integer = expect_and_consume_token(Token::Kind::nninteger);
        expect_and_consume_token(Token::Kind::r_paren);
        auto op = parse_qop();

        return ast::IfStmt::create(pos, id.as_string(), integer.as_int(),
                                   std::move(op));
    }
};

/**
 * \brief Parse a specified file
 */
inline ast::ptr<ast::Program> parse_file(std::string fname) {
    Preprocessor pp;
    Parser parser(pp);

    std::shared_ptr<std::ifstream> ifs(new std::ifstream);

    ifs->open(fname, std::ifstream::in);
    if (!ifs->good()) {
        ifs->close();
        std::cerr << "File \"" << fname << "\" not found!\n";
        throw ParseError();
    }

    pp.add_target_stream(ifs, fname);

    return parser.parse();
}

/**
 * \brief Parse input from stdin
 */
inline ast::ptr<ast::Program> parse_stdin(std::string name = "") {
    Preprocessor pp;
    Parser parser(pp);

    // This is a bad idea, but it's necessary for automatic bookkeeping
    // accross all different forms and sources of source streams
    pp.add_target_stream(
        std::shared_ptr<std::istream>(&std::cin, [](std::istream*) {}), name);

    return parser.parse();
}

/**
 * \brief Parse input stream
 */
inline ast::ptr<ast::Program> parse_stream(std::istream& stream) {
    Preprocessor pp;
    Parser parser(pp);

    // do not manage the stream, use [](std::istream*){} as shared_ptr deleter
    pp.add_target_stream(
        std::shared_ptr<std::istream>(&stream, [](std::istream*) {}));

    return parser.parse();
}

/**
 * \brief Parse a string
 * \note For small programs
 */
inline ast::ptr<ast::Program> parse_string(const std::string& str,
                                           std::string name = "") {
    Preprocessor pp;
    Parser parser(pp);
    std::shared_ptr<std::istream> is =
        std::make_shared<std::istringstream>(str);

    pp.add_target_stream(is, name);

    return parser.parse();
}

} /* namespace parser */
} /* namespace qasmtools */

#endif /* QASMTOOLS_PARSER_PARSER_HPP_ */
