/*
 * This file is part of Quantum++.
 *
 * MIT License
 *
 * Copyright (c) 2013 - 2019 Vlad Gheorghiu (vgheorgh@gmail.com)
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
 * \file qasm/parser.h
 * \brief openQASM 2.0 parsing
 */

#ifndef QASM_PARSER_H_
#define QASM_PARSER_H_

namespace qpp {
namespace qasm {

/**
 * \class qpp::qasm::Parser
 * \brief openQASM parser class
 * \see qpp::qasm::Lexer
 */
class Parser {
    Preprocessor& pp_lexer_; ///< preprocessed, tokenized input stream

    bool error_ = false;            ///< whether a parse error has occured
    bool supress_errors_ = false;   ///< whether to supress errors
    Token current_token_ = Token(); ///< current token
    int bits_ = 0;                  ///< number of bits
    int qubits_ = 0;                ///< number of qubits

  public:
    /**
     * \brief Constructs a parser for a given preprocessor stream
     *
     * \param pp_lexer The preprocessed lexer stream to be parsed
     */
    Parser(Preprocessor& pp_lexer) : pp_lexer_(pp_lexer) {}

    /**
     * \brief Parses the tokenized stream as a QCircuit object
     *
     * \return A unique pointer to a QCircuit object
     */
    std::unique_ptr<QCircuit> parse() {
        // Parse the program
        auto result = parse_program();
        if (error_)
            throw exception::ParseError("qpp::qasm::Parser::parse()");

        // Generate the QCircuit
        return result.to_QCircuit();
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
                std::cerr << current_token_.location();
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
                std::cerr << current_token_.location();
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
     * \brief Parse an openQASM 2.0 program
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
    QASM parse_program() {
        std::vector<StatementPtr> ret;

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
                        std::cerr << current_token_.location();
                        std::cerr << ": expected a declaration or statement";
                        std::cerr << " but got " << current_token_.kind() << "\n";
                        ;
                        supress_errors_ = true;
                    }

                    consume_until(Token::Kind::semicolon);
                    break;
            }
        }

        return QASM(bits_, qubits_, std::move(ret));
    }

    /**
     * \brief Parse an openQASM 2.0 header
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
    StatementPtr parse_reg_decl(bool quantum) {
        auto loc = current_token_.location();

        consume_token();
        auto identifier = expect_and_consume_token(Token::Kind::identifier);
        expect_and_consume_token(Token::Kind::l_square);
        auto size = expect_and_consume_token(Token::Kind::nninteger);
        expect_and_consume_token(Token::Kind::r_square);
        consume_until(Token::Kind::semicolon);

        quantum ? qubits_ += (int) size : bits_ += (int) size;
        return StatementPtr(new RegisterDecl(loc, identifier, quantum, size));
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
    StatementPtr parse_gate_decl() {
        auto loc = current_token_.location();

        consume_token();
        auto identifier = expect_and_consume_token(Token::Kind::identifier);

        std::vector<ident> c_params;
        if (try_and_consume_token(Token::Kind::l_paren)) {
            c_params = parse_idlist();
            expect_and_consume_token(Token::Kind::r_paren);
        }

        auto q_params = parse_idlist();

        expect_and_consume_token(Token::Kind::l_brace);
        std::vector<GatePtr> body;
        if (!try_and_consume_token(Token::Kind::r_brace)) {
            body = parse_goplist();
            expect_and_consume_token(Token::Kind::r_brace);
        }

        return StatementPtr(new GateDecl(loc, identifier, false, c_params,
                                         q_params, std::move(body)));
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
    StatementPtr parse_opaque_decl() {
        auto loc = current_token_.location();

        consume_token();
        auto identifier = expect_and_consume_token(Token::Kind::identifier);

        std::vector<ident> c_params;
        if (try_and_consume_token(Token::Kind::l_paren)) {
            c_params = parse_idlist();
            expect_and_consume_token(Token::Kind::r_paren);
        }

        auto q_params = parse_idlist();

        return StatementPtr(
            new GateDecl(loc, identifier, true, c_params, q_params, {}));
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
    std::vector<GatePtr> parse_goplist() {
        std::vector<GatePtr> ret;
        bool finished = false;

        while (!finished) {
            switch (current_token_.kind()) {
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
    StatementPtr parse_qop() {
        switch (current_token_.kind()) {
            case Token::Kind::kw_cx:
            case Token::Kind::kw_u:
            case Token::Kind::identifier:
                return StatementPtr(parse_gop());

            case Token::Kind::kw_measure:
                return parse_measure();

            case Token::Kind::kw_reset:
                return parse_reset();

            default:
                error_ = true;
                if (!supress_errors_) {
                    std::cerr << current_token_.location();
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
    GatePtr parse_gop() {
        switch (current_token_.kind()) {
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
                    std::cerr << current_token_.location();
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
    std::vector<ident> parse_idlist() {
        std::vector<ident> ret;

        // doesn't accept empty lists
        while (true) {
            auto identifier = expect_and_consume_token(Token::Kind::identifier);
            ret.push_back(identifier);
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
    Varinfo parse_argument() {
        auto loc = current_token_.location();

        auto id = expect_and_consume_token(Token::Kind::identifier);
        if (try_and_consume_token(Token::Kind::l_square)) {
            auto offset = expect_and_consume_token(Token::Kind::nninteger);
            expect_and_consume_token(Token::Kind::r_square);

            return Varinfo(loc, id, offset);
        }

        return Varinfo(loc, id);
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
    std::vector<Varinfo> parse_anylist() {
        std::vector<Varinfo> ret;

        // doesn't accept empty lists
        while (true) {
            auto arg = parse_argument();
            ret.push_back(arg);
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
    std::vector<ExprPtr> parse_explist() {
        // doesn't accept empty lists
        std::vector<ExprPtr> ret;

        while (true) {
            auto exp = parse_exp();
            ret.push_back(std::move(exp));
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
    ExprPtr parse_exp(int min_precedence = 1) {
        auto loc = current_token_.location();

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
                ExprPtr(new BExpr(loc, std::move(lexp), bop, std::move(rexp)));
        }

        return lexp;
    }

    /**
     * \brief Parse an atomic expression
     * \see qpp:qasm::Parser::parse_exp()
     *
     * \return Unique pointer to an expression object
     */
    ExprPtr parse_atom() {
        auto loc = current_token_.location();

        switch (current_token_.kind()) {
            case Token::Kind::l_paren: {
                consume_token();
                auto exp = parse_exp();
                expect_and_consume_token(Token::Kind::r_paren);
                return exp;
            }

            case Token::Kind::minus: {
                consume_token();
                auto exp = parse_exp();
                return ExprPtr(new UExpr(loc, UnaryOp::Neg, std::move(exp)));
            }

            case Token::Kind::identifier: {
                auto identifier = current_token_;
                consume_token();
                return ExprPtr(new VarExpr(loc, identifier));
            }
            case Token::Kind::nninteger: {
                auto integer = current_token_;
                consume_token();
                return ExprPtr(new IntExpr(loc, integer));
            }
            case Token::Kind::kw_pi:
                consume_token();
                return ExprPtr(new PiExpr(loc));
            case Token::Kind::real: {
                auto real = current_token_;
                consume_token();
                return ExprPtr(new RealExpr(loc, real));
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
                return ExprPtr(new UExpr(loc, op, std::move(exp)));
            }

            default:
                error_ = true;
                if (!supress_errors_) {
                    std::cerr << current_token_.location();
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
     * \see qpp:qasm::Parser::parse_exp()
     *
     * \return Binary operator
     */
    BinaryOp parse_binaryop() {
        switch (current_token_.kind()) {
            case Token::Kind::plus:
                consume_token();
                return BinaryOp::Plus;
            case Token::Kind::minus:
                consume_token();
                return BinaryOp::Minus;
            case Token::Kind::star:
                consume_token();
                return BinaryOp::Times;
            case Token::Kind::slash:
                consume_token();
                return BinaryOp::Divide;
            case Token::Kind::caret:
                consume_token();
                return BinaryOp::Pow;
            default:
                error_ = true;
                if (!supress_errors_) {
                    std::cerr << current_token_.location();
                    std::cerr << ": expected a binary operator but got ";
                    std::cerr << current_token_.kind() << "\n";
                    ;
                    supress_errors_ = true;
                }

                return BinaryOp::Plus;
        }
    }

    /**
     * \brief Parse a unary operator
     * \see qpp:qasm::Parser::parse_exp()
     *
     * \return Unary operator
     */
    UnaryOp parse_unaryop() {
        switch (current_token_.kind()) {
            case Token::Kind::kw_sin:
                consume_token();
                return UnaryOp::Sin;
            case Token::Kind::kw_cos:
                consume_token();
                return UnaryOp::Cos;
            case Token::Kind::kw_tan:
                consume_token();
                return UnaryOp::Tan;
            case Token::Kind::kw_exp:
                consume_token();
                return UnaryOp::Exp;
            case Token::Kind::kw_ln:
                consume_token();
                return UnaryOp::Ln;
            case Token::Kind::kw_sqrt:
                consume_token();
                return UnaryOp::Sqrt;
            default:
                error_ = true;
                if (!supress_errors_) {
                    std::cerr << current_token_.location();
                    std::cerr << ": expected a unary operator but got ";
                    std::cerr << current_token_.kind() << "\n";
                    ;
                    supress_errors_ = true;
                }

                return UnaryOp::Neg;
        }
    }

    /**
     * \brief Parse a CNOT gate
     *
     * CX <var> , <var> ;
     *
     * \return Unique pointer to a gate object
     */
    GatePtr parse_cnot() {
        auto loc = current_token_.location();

        expect_and_consume_token(Token::Kind::kw_cx);
        auto ctrl = parse_argument();
        expect_and_consume_token(Token::Kind::comma);
        auto tgt = parse_argument();
        consume_until(Token::Kind::semicolon);

        return GatePtr(new CNOTGate(loc, ctrl, tgt));
    }

    /**
     * \brief Parse a single qubit U gate
     *
     * U ( <exp> , <exp> , <exp> ) <var> ;
     *
     * \return Unique pointer to a gate object
     */
    GatePtr parse_unitary() {
        auto loc = current_token_.location();

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

        return GatePtr(new UGate(loc, std::move(theta), std::move(phi),
                                 std::move(lambda), arg));
    }

    /**
     * \brief Parse a declared gate application
     *
     * <id> <varlist> ;
     * <id> ( <explist> ) <varlist> ;
     *
     * \return Unique pointer to a gate object
     */
    GatePtr parse_gate_statement() {
        auto loc = current_token_.location();

        auto id = expect_and_consume_token(Token::Kind::identifier);
        std::vector<ExprPtr> c_args;
        if (try_and_consume_token(Token::Kind::l_paren)) {
            c_args = parse_explist();
            expect_and_consume_token(Token::Kind::r_paren);
        }

        auto q_args = parse_anylist();
        consume_until(Token::Kind::semicolon);

        return GatePtr(new DeclaredGate(loc, id, std::move(c_args), q_args));
    }

    /**
     * \brief Parse a measurement
     *
     * measure <var> -> <var> ;
     *
     * \return Unique pointer to a statement object
     */
    StatementPtr parse_measure() {
        auto loc = current_token_.location();

        expect_and_consume_token(Token::Kind::kw_measure);
        auto q_arg = parse_argument();
        expect_and_consume_token(Token::Kind::arrow);
        auto c_arg = parse_argument();
        consume_until(Token::Kind::semicolon);

        return StatementPtr(new MeasureStatement(loc, q_arg, c_arg));
    }

    /**
     * \brief Parse a reset statement
     *
     * reset <var> ;
     *
     * \return Unique pointer to a statement object
     */
    StatementPtr parse_reset() {
        auto loc = current_token_.location();

        expect_and_consume_token(Token::Kind::kw_reset);
        auto arg = parse_argument();
        consume_until(Token::Kind::semicolon);

        return StatementPtr(new ResetStatement(loc, arg));
    }

    /**
     * \brief Parse a barrier statement
     *
     * barrier <varlist> ;
     *
     * \return Unique pointer to a gate object
     */
    GatePtr parse_barrier() {
        auto loc = current_token_.location();

        expect_and_consume_token(Token::Kind::kw_barrier);
        auto args = parse_anylist();
        consume_until(Token::Kind::semicolon);

        return GatePtr(new BarrierGate(loc, args));
    }

    /**
     * \brief Parse an if statement
     *
     * if ( <id> == <nninteger> ) <qop>
     *
     * \return Unique pointer to a statement object
     */
    StatementPtr parse_if() {
        auto loc = current_token_.location();

        expect_and_consume_token(Token::Kind::kw_if);
        expect_and_consume_token(Token::Kind::l_paren);
        auto identifier = expect_and_consume_token(Token::Kind::identifier);
        expect_and_consume_token(Token::Kind::equalequal);
        auto integer = expect_and_consume_token(Token::Kind::nninteger);
        expect_and_consume_token(Token::Kind::r_paren);
        auto op = parse_qop();

        return StatementPtr(
            new IfStatement(loc, identifier, integer, std::move(op)));
    }
};

} /* namespace qasm */
} /* namespace qpp */

#endif
