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
 * \file qasm/qasm.hpp
 * \brief openQASM to QCircuit interface
 */

#ifndef QASM_QASM_HPP_
#define QASM_QASM_HPP_

namespace qpp {
namespace qasm {

/**
 * \brief Reads a openQASM circuit from stdin and returns its qpp::QCircuit
 * representation
 *
 * \return qpp::QCircuit
 */
inline QCircuit read(std::istream& stream) {
    Preprocessor pp;
    Parser parser(pp);

    // do not manage the stream, use [](std::istream*){} as shared_ptr deleter
    pp.add_target_stream(
        std::shared_ptr<std::istream>(&stream, [](std::istream*) {}));

    return *parser.parse();
}

/**
 * \brief Reads a openQASM circuit from a file and returns its qpp::QCircuit
 * representation
 *
 * \return qpp::QCircuit
 */
inline QCircuit read_from_file(const std::string& fname) {
    Preprocessor pp;
    Parser parser(pp);

    std::shared_ptr<std::ifstream> ifs(new std::ifstream);

    // EXCEPTION CHECKS

    ifs->open(fname, std::ifstream::in);
    if (!ifs->good()) {
        ifs->close();
        throw exception::FilesystemError("qpp::qasm::read_from_file()");
    }
    // END EXCEPTION CHECKS

    pp.add_target_stream(ifs, fname);

    return *parser.parse();
}

} /* namespace qasm */
} /* namespace qpp */

#endif /* QASM_QASM_HPP_ */
