/*
 * This file is part of Quantum++.
 *
 * MIT License
 *
 * Copyright (c) 2013 - 2020 Vlad Gheorghiu.
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
 * \file experimental/experimental.h
 * \brief Experimental/test functions/classes
 */

#ifndef EXPERIMENTAL_EXPERIMENTAL_H_
#define EXPERIMENTAL_EXPERIMENTAL_H_

namespace qpp {
/**
 * \namespace qpp::experimental
 * \brief Experimental/test functions/classes, do not use or modify
 */
namespace experimental {
/**
 * \brief Quantumly-accessible Random Access Memory (qRAM) over classical data,
 * implements \f$\sum_j\alpha_j|j\rangle\stackrel{qRAM}{\longrightarrow}
   \sum_j\alpha_j|j\rangle|m_j\rangle\f$
 *
 * \see
 * <a href="https://iopscience.iop.org/article/10.1088/1367-2630/17/12/123010">
 * New Journal of Physics, Vol. 17, No. 12, Pp. 123010 (2015)</a>
 *
 * \param in Input amplitudes
 * \param data Vector storing the classical data
 * \param DqRAM qRAM subsystem dimension (optional, by default is set to 1 +
 * maximum value stored in the qRAM)
 * \return Superposition over the qRAM values
 */
ket qRAM(const ket& in, const qram& data, idx DqRAM = -1) {
    // EXCEPTION CHECKS

    // check zero-size
    if (!internal::check_nonzero_size(in))
        throw exception::ZeroSize("qpp::experimental::qRAM()");

    // check equal dimensions
    idx Din = static_cast<idx>(in.rows());
    if (data.size() < Din)
        throw exception::DimsInvalid("qpp::experimental::qRAM()");

    idx max_val_qRAM = *std::max_element(std::begin(data), std::end(data));
    if (DqRAM != static_cast<idx>(-1) && DqRAM <= max_val_qRAM)
        throw exception::DimsInvalid("qpp::experimental::qRAM()");
    // END EXCEPTION CHECKS

    if (DqRAM == static_cast<idx>(-1))
        DqRAM = max_val_qRAM + 1;
    idx Dout = Din * DqRAM;
    ket result(Dout);

#ifdef WITH_OPENMP_
#pragma omp parallel for
#endif // WITH_OPENMP_
    for (idx i = 0; i < Din; ++i)
        result.block(i * DqRAM, 0, DqRAM, 1) = in[i] * st.j(data[i], DqRAM);

    return result;
}
} /* namespace experimental */
} /* namespace qpp */

#endif /* EXPERIMENTAL_EXPERIMENTAL_H_ */
