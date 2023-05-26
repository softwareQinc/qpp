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
 * \file qasmtools/parser/token.hpp
 * \brief Tokens
 */

#ifndef QASMTOOLS_PARSER_TOKEN_HPP_
#define QASMTOOLS_PARSER_TOKEN_HPP_

#include <unordered_map>
#include <variant>

#include "position.hpp"

namespace qasmtools {
namespace parser {

/**
 * \class qasmtools::parser::Token
 * \brief OpenQASM token class
 * \see qasmtools::parser::Lexer
 */
class Token {
  public:
    /**
     * \brief Token types
     */
    enum class Kind {
        error,
        eof,
        comment,
        identifier,
        real,
        nninteger,
        string,
        l_square,
        r_square,
        l_paren,
        r_paren,
        l_brace,
        r_brace,
        period,
        star,
        plus,
        minus,
        arrow,
        slash,
        caret,
        semicolon,
        equalequal,
        comma,
        colon,
        kw_include,
        kw_barrier,
        kw_creg,
        kw_cx,
        kw_gate,
        kw_if,
        kw_measure,
        kw_pi,
        kw_opaque,
        kw_openqasm,
        kw_qreg,
        kw_reset,
        kw_u,
        kw_sin,
        kw_cos,
        kw_tan,
        kw_exp,
        kw_ln,
        kw_sqrt,
        kw_oracle,
        kw_dirty,
        kw_ancilla
    };

    /**
     * \brief Extraction operator overload
     *
     * \param os Output stream passed by reference
     * \param k qasmtools::parser::Kind enum class
     * \return Reference to the output stream
     */
    friend std::ostream& operator<<(std::ostream& os, const Kind& k) {
        switch (k) {
            case Kind::error:
                os << "ERROR";
                break;
            case Kind::eof:
                os << "EOF";
                break;
            case Kind::comment:
                os << "COMMENT";
                break;
            case Kind::identifier:
                os << "ID";
                break;
            case Kind::real:
                os << "REAL";
                break;
            case Kind::nninteger:
                os << "INT";
                break;
            case Kind::string:
                os << "STRING";
                break;
            case Kind::l_square:
                os << "LBRACKET";
                break;
            case Kind::r_square:
                os << "RBRACKET";
                break;
            case Kind::l_paren:
                os << "LPAREN";
                break;
            case Kind::r_paren:
                os << "RPAREN";
                break;
            case Kind::l_brace:
                os << "LBRACE";
                break;
            case Kind::r_brace:
                os << "RBRACE";
                break;
            case Kind::period:
                os << "PERIOD";
                break;
            case Kind::star:
                os << "STAR";
                break;
            case Kind::plus:
                os << "PLUS";
                break;
            case Kind::minus:
                os << "MINUS";
                break;
            case Kind::arrow:
                os << "RARROW";
                break;
            case Kind::slash:
                os << "SLASH";
                break;
            case Kind::caret:
                os << "CARAT";
                break;
            case Kind::semicolon:
                os << "SEMICOLON";
                break;
            case Kind::equalequal:
                os << "EQUALS";
                break;
            case Kind::comma:
                os << "COMMA";
                break;
            case Kind::colon:
                os << "COLON";
                break;
            case Kind::kw_include:
                os << "INCLUDE";
                break;
            case Kind::kw_barrier:
                os << "BARRIER";
                break;
            case Kind::kw_creg:
                os << "CREG";
                break;
            case Kind::kw_cx:
                os << "CX";
                break;
            case Kind::kw_gate:
                os << "GATE";
                break;
            case Kind::kw_if:
                os << "IF";
                break;
            case Kind::kw_measure:
                os << "MEASURE";
                break;
            case Kind::kw_pi:
                os << "PI";
                break;
            case Kind::kw_opaque:
                os << "OPAQUE";
                break;
            case Kind::kw_openqasm:
                os << "OPENQASM";
                break;
            case Kind::kw_qreg:
                os << "QREG";
                break;
            case Kind::kw_reset:
                os << "RESET";
                break;
            case Kind::kw_u:
                os << "UNITARY";
                break;
            case Kind::kw_sin:
                os << "SIN";
                break;
            case Kind::kw_cos:
                os << "COS";
                break;
            case Kind::kw_tan:
                os << "TAN";
                break;
            case Kind::kw_exp:
                os << "EXP";
                break;
            case Kind::kw_ln:
                os << "LN";
                break;
            case Kind::kw_sqrt:
                os << "SQRT";
                break;
            case Kind::kw_oracle:
                os << "ORACLE";
                break;
            case Kind::kw_ancilla:
                os << "ANCILLA";
                break;
            case Kind::kw_dirty:
                os << "DIRTY";
                break;
        }
        return os;
    }

    Token() = delete;
    Token(Position pos, Kind k, const std::string& str,
          const std::variant<int, double, std::string>& value = {})
        : pos_(pos), kind_(k), str_(str), value_(value) {}

    /**
     * \brief Contextually return the token kind
     *
     * Allows switching directly on the token
     */
    operator Kind() const { return kind_; }

    /**
     * \brief Get the type of token
     *
     * \return The token type
     */
    Kind kind() const { return kind_; }

    /**
     * \brief Check whether the token is a particular type
     *
     * \param k The token type to compare against
     * \return True if and only if the token is of type k
     */
    bool is(Kind k) const { return kind_ == k; }

    /**
     * \brief Check whether the token is not a particular type
     *
     * \param k The token type to compare against
     * \return True if and only if the token is not of type k
     */
    bool is_not(Kind k) const { return kind_ != k; }

    /**
     * \brief User-defined conversion to int
     *
     * \note Does not perform validity checks
     *
     * \return The value of the token as an integer
     */
    int as_int() const { return std::get<int>(value_); }

    /**
     * \brief Get the floating point value
     *
     * \note Does not perform validity checks
     *
     * \return The value of the token as a floating point number
     */
    double as_real() const { return std::get<double>(value_); }

    /**
     * \brief Get the string value
     *
     * \note Does not perform validity checks
     *
     * \return The value of the token as a string
     */
    std::string as_string() const { return std::get<std::string>(value_); }

    /**
     * \brief Return the position of the token
     *
     * \return Const reference to the token position
     */
    const Position& position() const { return pos_; }

    /**
     * \brief Get the raw string
     *
     * \return The raw source string
     */
    const std::string& raw() const { return str_; }

    /**
     * \brief Extraction operator override
     *
     * Writes to the output stream a textual representation of the token
     *
     * \param os Output stream passed by reference
     * \return Reference to the output stream
     */
    friend std::ostream& operator<<(std::ostream& os, const Token& t) {
        os << t.kind_;
        switch (t) {
            case Kind::identifier:
            case Kind::real:
            case Kind::nninteger:
            case Kind::string:
                std::visit([&os](const auto& v) { os << "(" << v << ")"; },
                           t.value_);
                break;
            default:
                break;
        }
        return os;
    }

  private:
    Position pos_;                                 ///< the token position
    Kind kind_;                                    ///< the token type
    std::string str_;                              ///< the raw string
    std::variant<int, double, std::string> value_; ///< the token value
};

/**
 * \brief Hash-map of OpenQASM keywords and their token type
 */
static const std::unordered_map<std::string, Token::Kind> keywords{
    {"include", Token::Kind::kw_include},
    {"barrier", Token::Kind::kw_barrier},
    {"creg", Token::Kind::kw_creg},
    {"CX", Token::Kind::kw_cx},
    {"gate", Token::Kind::kw_gate},
    {"if", Token::Kind::kw_if},
    {"measure", Token::Kind::kw_measure},
    {"pi", Token::Kind::kw_pi},
    {"opaque", Token::Kind::kw_opaque},
    {"OPENQASM", Token::Kind::kw_openqasm},
    {"qreg", Token::Kind::kw_qreg},
    {"reset", Token::Kind::kw_reset},
    {"U", Token::Kind::kw_u},
    {"sin", Token::Kind::kw_sin},
    {"cos", Token::Kind::kw_cos},
    {"tan", Token::Kind::kw_tan},
    {"exp", Token::Kind::kw_exp},
    {"ln", Token::Kind::kw_ln},
    {"sqrt", Token::Kind::kw_sqrt},
    {"oracle", Token::Kind::kw_oracle},
    {"dirty", Token::Kind::kw_dirty},
    {"ancilla", Token::Kind::kw_ancilla}};

} /* namespace parser */
} /* namespace qasmtools */

#endif /* QASMTOOLS_PARSER_TOKEN_HPP_ */
