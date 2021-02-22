/*
 * This file is part of Quantum++.
 *
 * Copyright (c) 2013 - 2021 softwareQ Inc. All rights reserved.
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
 * \file qasm/token.hpp
 * \brief openQASM tokens
 */

#ifndef QASM_TOKEN_HPP_
#define QASM_TOKEN_HPP_

namespace qpp {
namespace qasm {

/**
 * \class qpp::qasm::Location
 * \brief Source location class
 */
class Location : public IDisplay {
    std::string fname_{}; ///< name of the containing file
    idx line_ = 1;        ///< line number
    idx column_ = 1;      ///< column number

  public:
    /**
     * \brief Default constructor
     */
    Location() = default;

    /**
     * \brief Constructs a location within a file
     *
     * \param fname Filename
     * \param line Line number
     * \param column Column number
     */
    Location(std::string fname, idx line, idx column)
        : fname_(std::move(fname)), line_(line), column_(column) {}

    /**
     * \brief qpp::IDisplay::display() override
     *
     * Writes to the output stream a textual representation of the
     * source location
     *
     * \param os Output stream passed by reference
     * \return Reference to the output stream
     */
    std::ostream& display(std::ostream& os) const override {
        os << fname_ << ":" << line_ << ":" << column_;
        return os;
    }

    /**
     * \brief Default assignment operator
     *
     * \return Reference to the current instance
     */
    Location& operator=(const Location&) = default;

    /**
     * \brief The name of the containing file
     *
     * \return Const reference to the filename
     */
    const std::string& get_filename() const { return fname_; }

    /**
     * \brief The line of the location
     *
     * \return The line number
     */
    idx get_linenum() const { return line_; }

    /**
     * \brief The column of the location
     *
     * \return The column of the location
     */
    idx get_colnum() const { return column_; }

    /**
     * \brief Advance the line number by a specified amount
     *
     * \note Sets the column to 0
     *
     * \param num Number of lines to advance (optional, default is 1)
     */
    void advance_line(idx num = 1) {
        line_ += num;
        column_ = 0;
    }

    /**
     * \brief Advance the column number by a specified amount
     *
     * \param num Number of columns to advance (optional, default is 1)
     */
    void advance_column(idx num = 1) { column_ += num; }
};

/**
 * \class qpp::qasm::Token
 * \brief openQASM token class
 * \see qpp::qasm::Lexer
 */
class Token : public IDisplay {
  public:
    /**
     * \brief Token types
     */
    enum class Kind {
        unknown,
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
    };

    /**
     * \brief Default constructor
     */
    Token() = default;

    /**
     * \brief Constructs a token with a given type, source location, and string
     * value
     *
     * \param loc The location of the token within a source stream
     * \param k The type of token being created
     * \param value The string value of the token
     *
     */
    Token(const Location& loc, Kind k, std::string value)
        : loc_(loc), kind_(k), value_(std::move(value)) {}

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
     * \brief Check whether the token is one of a list of types
     *
     * \param k1 First token types to check
     * \param k2 Second token types to check
     * \param ks List of remaining token types to check
     * \return True if and only if the token is one of the types in ks
     */
    template <typename... Ts>
    bool is_one_of(Kind k1, Kind k2, Ts... ks) const {
        return is(k1) || is_one_of(k2, ks...);
    }

    /**
     * \brief User-defined conversion to double
     *
     * \note Does not perform validity checks
     *
     * \return The value of the token as a floating point number
     */
    operator double() { return std::stod(value_); }

    /**
     * \brief User-defined conversion to std::string
     *
     * \return The value of the token as a string
     */
    operator std::string() { return std::string(value_); }

    /**
     * \brief User-defined conversion to idx
     *
     * \note Does not perform validity checks
     *
     * \return The value of the token as an unsigned integer
     */
    operator idx() { return static_cast<idx>(std::stoi(value_)); }

    /**
     * \brief User-defined conversion to int
     *
     * \note Does not perform validity checks
     *
     * \return The value of the token as an integer
     */
    operator int() { return std::stoi(value_); }

    /**
     * \brief Return the location of the token
     *
     * \return Const reference to the token location
     */
    const Location& location() const { return loc_; }

    /**
     * \brief Extraction operator overload for qpp::qasm::Token::Kind enum class
     *
     * \param os Output stream passed by reference
     * \param k qpp::qasm::Token::Kind enum class
     * \return Reference to the output stream
     */
    friend std::ostream& operator<<(std::ostream& os, const Kind& k) {
        switch (k) {
            case Kind::unknown:
                os << "UNKNOWN";
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
        }
        return os;
    }

    /**
     * \brief qpp::IDisplay::display() override
     *
     * Writes to the output stream a textual representation of the token
     *
     * \param os Output stream passed by reference
     * \return Reference to the output stream
     */
    std::ostream& display(std::ostream& os) const override {
        os << kind_;
        switch (kind_) {
            case Kind::identifier:
                os << "(" << value_ << ")";
                break;
            case Kind::real:
                os << "(" << std::stof(value_) << ")";
                break;
            case Kind::nninteger:
                os << "(" << std::stoi(value_) << ")";
                break;
            case Kind::string:
                os << "(\"" << value_ << "\")";
                break;
            default:
                break;
        }
        return os;
    }

  private:
    Location loc_ = Location(); ///< the token location
    Kind kind_ = Kind::unknown; ///< the token type
    std::string value_{};       ///< the token value
};

/**
 * \brief Hash-map of openQASM keywords and their token type
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
    {"sqrt", Token::Kind::kw_sqrt}};

} /* namespace qasm */
} /* namespace qpp */

#endif /* QASM_TOKEN_HPP_ */
