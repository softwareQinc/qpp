/*
 * This file is part of pyqpp.
 *
 * Copyright (c) 2017 - 2026 softwareQ Inc. All rights reserved.
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
 * \file <pyqpp/number_theory_bind.hpp>
 * \brief Bindings for <qpp/number_theory.hpp>
 */

#ifndef PYQPP_NUMBER_THEORY_BIND_HPP_
#define PYQPP_NUMBER_THEORY_BIND_HPP_

#include "pyqpp/pyqpp_common.hpp"

inline void init_number_theory(py::module_& m) {
    using namespace qpp;

    // --- 1. Continued Fractions ---
    m.def("x2contfrac", &qpp::x2contfrac, py::arg("x"), py::arg("N"),
          py::arg("cut") = 10000, "Simple continued fraction expansion of x");

    m.def("contfrac2x", &qpp::contfrac2x, py::arg("cf"),
          py::arg("N") = std::numeric_limits<idx>::max(),
          "Real representation of a simple continued fraction");

    m.def("convergents",
          py::overload_cast<const std::vector<bigint>&>(&qpp::convergents),
          "Vector of convergent pairs (a_k, b_k) from a continued fraction");

    m.def("convergents", py::overload_cast<realT, idx>(&qpp::convergents),
          "Vector of convergent pairs (a_k, b_k) approximating real x");

    // --- 2. GCD and LCM ---
    // Note: Python's math.gcd only takes 2 args (pre-3.9), these handle lists
    // too.
    m.def("gcd", py::overload_cast<bigint, bigint>(&qpp::gcd),
          "Greatest common divisor of two integers");
    m.def("gcd", py::overload_cast<const std::vector<bigint>&>(&qpp::gcd),
          "GCD of a list of integers");

    m.def("lcm", py::overload_cast<bigint, bigint>(&qpp::lcm),
          "Least common multiple of two integers");
    m.def("lcm", py::overload_cast<const std::vector<bigint>&>(&qpp::lcm),
          "LCM of a list of integers");

    m.def(
        "egcd",
        [](bigint a, bigint b) {
            auto result = qpp::egcd(a, b);
            return std::make_tuple(std::get<0>(result), std::get<1>(result),
                                   std::get<2>(result));
        },
        "Extended GCD. Returns (m, n, gcd) such that ma + nb = gcd(a, b).");

    // --- 3. Modular Arithmetic ---
    m.def("modmul", &qpp::modmul,
          "Modular multiplication (a*b % p) without overflow");
    m.def("modpow", &qpp::modpow, "Fast modular exponentiation (a^n % p)");
    m.def("modinv", &qpp::modinv,
          "Modular inverse of a mod p (a and p must be co-prime)");

    // --- 4. Primes and Factors ---
    m.def("factors", &qpp::factors, "Prime factor decomposition of an integer");
    m.def("isprime", &qpp::isprime, py::arg("p"), py::arg("k") = 80,
          "Miller-Rabin primality test");
    m.def("randprime", &qpp::randprime, py::arg("a"), py::arg("b"),
          py::arg("N") = 1000,
          "Generates a random prime in the interval [a, b]");

    // --- 5. Permutations ---
    m.def("invperm", &qpp::invperm, "Inverse of a permutation");
    m.def("compperm", &qpp::compperm, "Composition of two permutations");
}

#endif /* PYQPP_NUMBER_THEORY_BIND_HPP_ */
