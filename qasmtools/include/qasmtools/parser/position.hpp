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
 * \file qasmtools/parser/position.hpp
 */

#ifndef QASMTOOLS_PARSER_POSITION_HPP_
#define QASMTOOLS_PARSER_POSITION_HPP_

#include <iostream>
#include <string>

namespace qasmtools {
namespace parser {

/**
 * \class qasmtools::parser::Position
 * \brief Positions in source code
 */
class Position {
    std::string fname_ = ""; ///< name of the containing file
    int line_ = 1;           ///< line number
    int column_ = 1;         ///< column number

  public:
    /**
     * \brief Default constructor
     */
    Position() = default;

    /**
     * \brief Constructs a position within a file
     *
     * \param fname Filename
     * \param line Line number
     * \param column Column number
     */
    Position(const std::string& fname, int line, int column)
        : fname_(fname), line_(line), column_(column) {}

    /**
     * \brief Extraction operator overload
     *
     * \param os Output stream
     * \param pos qasmtools::parser::Position
     * \return Reference to the output stream
     */
    friend std::ostream& operator<<(std::ostream& os, const Position& pos) {
        os << pos.fname_ << ":" << pos.line_ << ":" << pos.column_;
        return os;
    }

    /**
     * \brief The name of the containing file
     *
     * \return Const reference to the filename
     */
    const std::string& get_filename() const { return fname_; }

    /**
     * \brief The line of the position
     *
     * \return The line number
     */
    int get_linenum() const { return line_; }

    /**
     * \brief The column of the position
     *
     * \return The column number
     */
    int get_column() const { return column_; }

    /**
     * \brief Advance the line number by a specified amount
     *
     * \note Sets the column to 0
     *
     * \param num Number of lines to advance (optional, default is 1)
     */
    void advance_line(int num = 1) {
        line_ += num;
        column_ = 1;
    }

    /**
     * \brief Advance the column number by a specified amount
     *
     * \param num Number of columns to advance (optional, default is 1)
     */
    void advance_column(int num = 1) { column_ += num; }
};

} /* namespace parser */
} /* namespace qasmtools */

#endif /* QASMTOOLS_PARSER_POSITION_HPP_ */
