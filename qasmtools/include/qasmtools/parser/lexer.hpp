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
 *
 * Adapted from Bruno Schmitt's Tweedledee library
 */

/**
 * \file qasmtools/parser/lexer.hpp
 * \brief Lexical analysis
 */

#pragma once

#include "token.hpp"

#include <cctype>
#include <memory>

namespace qasmtools {
namespace parser {

/**
 * \class qasmtools::parser::Lexer
 * \brief openPARSER lexer class
 *
 * The Lexer reads from (a shared_ptr to) an istream object given during
 * initialization. Rather than lex the entire buffer at once, tokens are
 * lexed and returned on-demand
 */
class Lexer {
  public:
    Lexer(const Lexer&) = delete;
    Lexer& operator=(const Lexer&) = delete;

    Lexer(std::shared_ptr<std::istream> buffer, const std::string& fname = "")
        : pos_(fname, 1, 1), buf_(buffer) {}

    /**
     * \brief Lex and return the next token
     *
     * \note Advances the buffer to the position after the consumed token
     *
     * \return The token that was lexed
     */
    Token next_token() { return lex(); }

  private:
    Position pos_; ///< current position in the source stream
    std::shared_ptr<std::istream> buf_; ///< stream buffer being lexed

    /**
     * \brief Skips the specified number of characters
     *
     * \param n The number of characters to skip ahead (optional, default is 1)
     */
    void skip_char(int n = 1) {
        buf_->ignore(n);
        pos_.advance_column(n);
    }

    /**
     * \brief Skips over whitespace
     *
     * \return True if and only if whitespace was actually consumed
     */
    bool skip_whitespace() {
        int consumed = 0;

        while (buf_->peek() == ' ' || buf_->peek() == '\t') {
            buf_->ignore();
            ++consumed;
        }

        pos_.advance_column(consumed);
        if (consumed != 0)
            return true;
        else
            return false;
    }

    /**
     * \brief Skips the rest of the line
     */
    void skip_line_comment() {
        int consumed = 0;

        int c = buf_->peek();
        while (c != 0 && c != '\n' && c != '\r' && c != EOF) {
            buf_->ignore();
            ++consumed;
            c = buf_->peek();
        }

        pos_.advance_column(consumed);
    }

    /**
     * \brief Lex a numeric constant
     *
     * \note [0-9]+(.[0-9]*)([eE][+-][0-9]+)?
     *
     * \param tok_start The position of the beginning of the token
     * \return An integer or real type token
     */
    Token lex_numeric_constant(Position tok_start) {
        std::string str;
        str.reserve(64); // Reserve space to avoid reallocation
        bool integral = true;

        while (std::isdigit(buf_->peek())) {
            str.push_back(buf_->peek());
            skip_char();
        }

        // lex decimal
        if (buf_->peek() == '.') {
            integral = false;
            str.push_back(buf_->peek());
            skip_char();

            while (std::isdigit(buf_->peek())) {
                str.push_back(buf_->peek());
                skip_char();
            }
        }

        // lex exponent
        if (buf_->peek() == 'e' || buf_->peek() == 'E') {
            integral = false;
            str.push_back(buf_->peek());
            skip_char();

            if (buf_->peek() == '-' || buf_->peek() == '+') {
                str.push_back(buf_->peek());
                skip_char();
            }

            while (std::isdigit(buf_->peek())) {
                str.push_back(buf_->peek());
                skip_char();
            }
        }

        if (integral) {
            return Token(tok_start, Token::Kind::nninteger, str,
                         std::stoi(str));
        } else {
            return Token(tok_start, Token::Kind::real, str, std::stof(str));
        }
    }

    /**
     * \brief Lex an identifier
     *
     * \note [A-Za-z][_A-Za-z0-9]*
     *
     * \param tok_start The position of the beginning of the token
     * \return An identifier type token
     */
    Token lex_identifier(Position tok_start) {
        std::string str;
        str.reserve(64); // Reserve space to avoid reallocation

        while (std::isalpha(buf_->peek()) || std::isdigit(buf_->peek()) ||
               buf_->peek() == '_') {
            str.push_back(buf_->peek());
            skip_char();
        }

        // Check if the identifier is a known keyword
        auto keyword = keywords.find(str);
        if (keyword != keywords.end()) {
            return Token(tok_start, keyword->second, str);
        }

        return Token(tok_start, Token::Kind::identifier, str, str);
    }

    /**
     * \brief Lex a string literal
     *
     * \note "[.]*"
     *
     * \param tok_start The position of the beginning of the token
     * \return A string type token
     */
    Token lex_string(Position tok_start) {
        std::string str;
        str.reserve(64); // Reserve space to avoid reallocation

        while (buf_->peek() != '"' && buf_->peek() != '\n' &&
               buf_->peek() != '\r') {
            str.push_back(buf_->peek());
            skip_char();
        }

        if (buf_->peek() != '"') {
            std::cerr << "Lexical error at " << tok_start << ": unmatched \"\n";
            return Token(tok_start, Token::Kind::error, str);
        }

        skip_char();
        return Token(tok_start, Token::Kind::string, str, str);
    }

    /**
     * \brief Lex a token
     *
     * \note See arXiv:1707.03429 for the full grammar
     *
     * \return The lexed token
     */
    Token lex() {
        Position tok_start = pos_;
        skip_whitespace();

        switch (buf_->peek()) {
            case EOF:
                skip_char();
                return Token(tok_start, Token::Kind::eof, "");

            case '\r':
                skip_char();
                if (buf_->peek() != '\n') {
                    buf_->ignore();
                    pos_.advance_line();
                    return lex();
                }
                // FALLTHROUGH
            case '\n':
                buf_->ignore();
                pos_.advance_line();
                return lex();

            case '/':
                skip_char();
                if (buf_->peek() == '/') {
                    skip_line_comment();
                    return lex();
                }
                return Token(tok_start, Token::Kind::slash, "/");

                // clang-format off
            case '0': case '1': case '2': case '3': case '4':
            case '5': case '6': case '7': case '8': case '9':
            case '.':
                return lex_numeric_constant(tok_start);
                // clang-format on

            case 'C':
                skip_char();
                if (buf_->peek() == 'X') {
                    skip_char();
                    return Token(tok_start, Token::Kind::kw_cx, "CX");
                }

                skip_char();
                std::cerr
                    << "Lexical error at " << tok_start
                    << ": identifiers must start with lowercase letters\n";
                return Token(tok_start, Token::Kind::error,
                             std::string({'C', (char) buf_->get()}));

            case 'U':
                skip_char();
                return Token(tok_start, Token::Kind::kw_u, "U");

                // clang-format off
            case 'O':
            case 'a': case 'b': case 'c': case 'd': case 'e': case 'f': case 'g':
            case 'h': case 'i': case 'j': case 'k': case 'l': case 'm': case 'n':
            case 'o': case 'p': case 'q': case 'r': case 's': case 't': case 'u':
            case 'v': case 'w': case 'x': case 'y': case 'z':
                return lex_identifier(tok_start);
                // clang-format on

            case '[':
                skip_char();
                return Token(tok_start, Token::Kind::l_square, "[");

            case ']':
                skip_char();
                return Token(tok_start, Token::Kind::r_square, "]");

            case '(':
                skip_char();
                return Token(tok_start, Token::Kind::l_paren, "(");

            case ')':
                skip_char();
                return Token(tok_start, Token::Kind::r_paren, ")");

            case '{':
                skip_char();
                return Token(tok_start, Token::Kind::l_brace, "{");

            case '}':
                skip_char();
                return Token(tok_start, Token::Kind::r_brace, "}");

            case '*':
                skip_char();
                return Token(tok_start, Token::Kind::star, "*");

            case '+':
                skip_char();
                return Token(tok_start, Token::Kind::plus, "+");

            case '-':
                skip_char();

                if (buf_->peek() == '>') {
                    skip_char();
                    return Token(tok_start, Token::Kind::arrow, "->");
                }

                return Token(tok_start, Token::Kind::minus, "-");

            case '^':
                skip_char();
                return Token(tok_start, Token::Kind::caret, "^");

            case ';':
                skip_char();
                return Token(tok_start, Token::Kind::semicolon, ";");

            case '=':
                skip_char();
                if (buf_->peek() == '=') {
                    skip_char();
                    return Token(tok_start, Token::Kind::equalequal, "==");
                }

                skip_char();
                std::cerr << "Lexical error at " << tok_start
                          << ": expected \"=\" after \"=\"\n";
                return Token(tok_start, Token::Kind::error,
                             std::string({'=', (char) buf_->get()}));

            case ',':
                skip_char();
                return Token(tok_start, Token::Kind::comma, ";");

            case '"':
                skip_char();
                return lex_string(tok_start);

            default:
                skip_char();
                return Token(tok_start, Token::Kind::error,
                             std::string({(char) buf_->get()}));
        }
    }
};

} // namespace parser
} // namespace qasmtools
