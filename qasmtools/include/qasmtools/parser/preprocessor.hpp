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
 *
 * Adapted from Bruno Schmitt's Tweedledee library
 */

/**
 * \file qasmtools/parser/preprocessor.hpp
 * \brief Manages includes for OpenQASM parsing
 */

#ifndef QASMTOOLS_PARSER_PREPROCESSOR_HPP_
#define QASMTOOLS_PARSER_PREPROCESSOR_HPP_

#include <fstream>
#include <sstream>
#include <vector>

#include "lexer.hpp"

namespace qasmtools {
namespace parser {

#if QASMTOOLS_QASM2_SPECS
/**
 * \brief OpenQASM 2.0 standard library (qelib1.inc) as a string constant
 */
static const std::string std_include =
    "gate u3(theta,phi,lambda) q { U(theta,phi,lambda) q; }\n"
    "gate u2(phi,lambda) q { U(pi/2,phi,lambda) q; }\n"
    "gate u1(lambda) q { U(0,0,lambda) q; }\n"
    "gate cx c,t { CX c,t; }\n"
    "gate id a { U(0,0,0) a; }\n"
    "gate u0(gamma) q { U(0,0,0) q; }\n"
    "gate x a { u3(pi,0,pi) a; }\n"
    "gate y a { u3(pi,pi/2,pi/2) a; }\n"
    "gate z a { u1(pi) a; }\n"
    "gate h a { u2(0,pi) a; }\n"
    "gate s a { u1(pi/2) a; }\n"
    "gate sdg a { u1(-pi/2) a; }\n"
    "gate t a { u1(pi/4) a; }\n"
    "gate tdg a { u1(-pi/4) a; }\n"
    "gate rx(theta) a { u3(theta, -pi/2,pi/2) a; }\n"
    "gate ry(theta) a { u3(theta,0,0) a; }\n"
    "gate rz(phi) a { u1(phi) a; }\n"
    "gate cz a,b { h b; cx a,b; h b; }\n"
    "gate cy a,b { sdg b; cx a,b; s b; }\n"
    "gate swap a,b { cx a,b; cx b,a; cx a,b; }\n"
    "gate ch a,b { h b; sdg b;cx a,b;h b; t b;cx a,b;t b; h b; s b; x b; s "
    "a;}\n"
    "gate ccx a,b,c{ h c; cx b,c; tdg c; cx a,c; t c; cx b,c; tdg c; cx a,c; t "
    "b; t c; h c; cx a,b; t a; tdg b; cx a,b;}\n"
    "gate crz(lambda) a,b{ u1(lambda/2) b; cx a,b; u1(-lambda/2) b; cx a,b;}\n"
    "gate cu1(lambda) a,b{ u1(lambda/2) a; cx a,b; u1(-lambda/2) b; cx a,b; "
    "u1(lambda/2) b;}\n"
    "gate cu3(theta,phi,lambda) c,t { u1((lambda-phi)/2) t; cx c,t; "
    "u3(-theta/2,0,-(phi+lambda)/2) t;  cx c,t;  u3(theta/2,phi,0) t;}\n";
#else
/**
 * \brief OpenQASM 2.0 standard library + r and cswap gates, as a string
 * constant, as defined by Qiskit in
 * https://github.com/Qiskit/qiskit-terra/tree/master/qiskit/circuit/library/standard_gates
 */
static const std::string std_include =
    "gate u3(theta,phi,lambda) q { U(theta,phi,lambda) q; }\n"
    "gate u2(phi,lambda) q { U(pi/2,phi,lambda) q; }\n"
    "gate u1(lambda) q { U(0,0,lambda) q; }\n"
    "gate cx c,t { CX c,t; }\n"
    "gate id a { U(0,0,0) a; }\n"
    "gate u0(gamma) q { U(0,0,0) q; }\n"
    "gate x a { u3(pi,0,pi) a; }\n"
    "gate y a { u3(pi,pi/2,pi/2) a; }\n"
    "gate z a { u1(pi) a; }\n"
    "gate h a { u2(0,pi) a; }\n"
    "gate s a { u1(pi/2) a; }\n"
    "gate sdg a { u1(-pi/2) a; }\n"
    "gate t a { u1(pi/4) a; }\n"
    "gate tdg a { u1(-pi/4) a; }\n"
    "gate r(theta, phi) a { u3(theta,phi-pi/2,-phi+pi/2) a; }\n"
    "gate rx(theta) a { r(theta,0) a; } \n"
    "gate ry(theta) a { r(theta,pi/2) a; }\n"
    /**
     * Rz = e^{-i*phi/2} U1
     * This differs from the docstring here:
     * https://github.com/Qiskit/qiskit-terra/blob/main/qiskit/circuit/library/standard_gates/rz.py#L65
     */
    "gate rz(phi) a { x a; u1(-phi/2) a; x a; u1(phi/2) a; } \n"
    "gate cz a,b { h b; cx a,b; h b; }\n"
    "gate cy a,b { sdg b; cx a,b; s b; }\n"
    "gate swap a,b { cx a,b; cx b,a; cx a,b; }\n"
    "gate ch a,b { s b; h b; t b; cx a,b; tdg b; h b; sdg b; }\n"
    "gate ccx a,b,c { h c; cx b,c; tdg c; cx a,c; t c; cx b,c; tdg c; cx a,c; "
    "t b; t c; h c; cx a,b; t a; tdg b; cx a,b; }\n"
    "gate cswap a,b,c { cx c,b; ccx a,b,c; cx c,b; }\n"
    "gate crz(lambda) a,b { u1(lambda/2) b; cx a,b; u1(-lambda/2) b; cx a,b; "
    "}\n"
    "gate cu1(lambda) a,b { u1(lambda/2) a; cx a,b; u1(-lambda/2) b; cx a,b; "
    "u1(lambda/2) b; }\n"
    "gate cu3(theta,phi,lambda) c,t { u1((lambda+phi)/2) c; "
    "u1((lambda-phi)/2) t; cx c,t; u3(-theta/2,0,-(phi+lambda)/2) t; cx c,t; "
    "u3(theta/2,phi,0) t; }\n";
#endif

/**
 * \class qasmtools::parser::Preprocessor
 * \brief OpenQASM preprocessor class
 * \see qasmtools::parser::Lexer
 *
 * The preprocessor acts as a wrapper around the lexer, providing a token stream
 * that matches the stream produced by explicitly inserting included code.
 * Effectively, the preprocessor acts as the lexer on preprocessed code.
 */
class Preprocessor {
    using LexerPtr = std::unique_ptr<Lexer>;

    std::vector<LexerPtr> lexer_stack_{}; ///< owning stack of lexers
    LexerPtr current_lexer_ = nullptr;    ///< current lexer

    bool std_include_ = false; ///< whether qelib1 has been included

  public:
    /**
     * \brief Default constructor
     */
    Preprocessor() = default;

    /**
     * \brief Inserts a file into the current lexing context
     *
     * Buffers the given file, then pushes the current buffer onto
     * the stack and sets the new buffer as the current buffer
     *
     * \param file_path The file to insert
     * \return True on success
     */
    bool add_target_file(const std::string& file_path) {
        std::shared_ptr<std::ifstream> ifs(new std::ifstream);

        ifs->open(file_path, std::ifstream::in);

        if (!ifs->good()) {
            return false;
        }

        if (current_lexer_ != nullptr) {
            lexer_stack_.push_back(std::move(current_lexer_));
        }
        current_lexer_ = std::unique_ptr<Lexer>(new Lexer(ifs, file_path));
        return true;
    }

    /**
     * \brief Inserts a buffer into the current lexing context
     *
     * Pushes the current buffer onto the stack and sets the new buffer as the
     * current buffer
     *
     * \param buffer Shared pointer to an input buffer
     * \param fname Filename associated with the buffer (optional)
     */
    void add_target_stream(std::shared_ptr<std::istream> buffer,
                           const std::string& fname = "") {
        if (current_lexer_ != nullptr) {
            lexer_stack_.push_back(std::move(current_lexer_));
        }
        current_lexer_ = std::unique_ptr<Lexer>(new Lexer(buffer, fname));
    }

    /**
     * \brief Gets the next token in the current buffer
     * \note Returns unknown token if there is no buffer to lex
     *
     * Lexes and returns the next token from the current buffer,
     * popping the next buffer off the stack if necessary
     *
     * \return The next lexed token
     */
    Token next_token() {
        if (current_lexer_ == nullptr) {
            return Token(Position(), Token::Kind::eof, "");
        }
        auto token = current_lexer_->next_token();
        if (token.is(Token::Kind::kw_include)) {
            handle_include();
            token = current_lexer_->next_token();
        } else if (token.is(Token::Kind::eof)) {
            if (!lexer_stack_.empty()) {
                current_lexer_ = std::move(lexer_stack_.back());
                lexer_stack_.pop_back();
                token = current_lexer_->next_token();
            } else {
                current_lexer_ = nullptr;
            }
        }
        return token;
    }

    /**
     * \brief Prints and consumes all tokens buffered
     */
    void print_all_tokens() {
        Token current = next_token();

        while (!current.is(Token::Kind::eof)) {
            std::cout << current.position() << ": " << current << " "
                      << current.raw() << "\n";
            current = next_token();
        }
    }

    bool includes_stdlib() { return std_include_; }

  private:
    /**
     * \brief Handles include statements
     * \note Expects STRING SEMICOLON
     *
     * Reads a filename from a string token and attempts to open and switch
     * lexing context to the file
     */
    void handle_include() {
        auto token = current_lexer_->next_token();
        if (token.is_not(Token::Kind::string)) {
            std::cerr << "Error: Include must be followed by a file name\n";
            return;
        }

        auto target = token.as_string();
        if (target == "qelib1.inc") {
            std_include_ = true;
        }

        token = current_lexer_->next_token();
        if (token.is_not(Token::Kind::semicolon)) {
            std::cerr << "Warning: Missing a ';'\n";
        }
        if (add_target_file(target)) {
            return;
        } else if (target == "qelib1.inc") {
            add_target_stream(std::shared_ptr<std::istream>(
                                  new std::stringstream(std_include)),
                              "qelib1.inc");
            return;
        } else {
            std::cerr << "Error: Couldn't open file " << target << "\n";
        }
    }
};

} /* namespace parser */
} /* namespace qasmtools */

#endif /* QASMTOOLS_PARSER_PREPROCESSOR_HPP_ */
