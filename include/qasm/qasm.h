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
 */

/**
 * \file qasm/qasm.h
 * \brief openQASM to QCircuit interface
 */

#ifndef QASM_QASM_H_
#define QASM_QASM_H_

namespace qpp {
namespace qasm {

/**
 * \brief Reads a openQASM circuit from stdin and returns its QCircuit
 * representation
 *
 * \return Unique pointer to a QCircuit object
 */
inline std::unique_ptr<QCircuit> read_from_stdin() {
    Preprocessor pp;
    Parser parser(pp);

    pp.add_target_stream(std::shared_ptr<std::istream>(&std::cin));

    return parser.parse();
}

/**
 * \brief Reads a openQASM circuit from a file and returns its QCircuit
 * representation
 *
 * \return Unique pointer to a QCircuit object
 */
inline std::unique_ptr<QCircuit> read_from_file(const std::string& fname_) {
    Preprocessor pp;
    Parser parser(pp);

    std::shared_ptr<std::ifstream> ifs(new std::ifstream);

    ifs->open(fname_, std::ifstream::in);
    if (!ifs->good()) {
        ifs->close();
        throw exception::FileNotFound("qpp::qasm::read_file");
    }

    pp.add_target_stream(ifs, fname_);

    return parser.parse();
}

} /* namespace qasm */
} /* namespace qpp */

#endif /* QASM_QASM_H_ */
