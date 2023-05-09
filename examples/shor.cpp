// Shor's algorithm
// Source: ./examples/shor.cpp
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <tuple>
#include <vector>

#include "qpp.h"

int main() {
    using namespace qpp;

    bigint N = 15;                   // number to factor
    auto a = rand<bigint>(3, N - 1); // random co-prime with N
    while (gcd(a, N) != 1) {
        a = rand<bigint>(3, N - 1);
    }
    // qubits required for half of the circuit, in total we need 2n qubits
    // if you know the order 'r' of 'a', then you can take the smallest 'n' s.t.
    // 2^n >= 2 * r^2, i.e., n = ceil(log2(2 * r^2))
    auto n = static_cast<idx>(std::ceil(2 * std::log2(N)));
    auto D = static_cast<idx>(std::llround(std::pow(2, n))); // dimension 2^n

    std::cout << ">> Factoring N = " << N << " with coprime a = " << a << '\n';
    std::cout << ">> n = " << n << ", D = 2^n = " << D << " and 2^(2n) = ";
    std::cout << D * D << '\n';

    // vector with labels of the first half of the qubits
    std::vector<idx> first_subsys(n);
    std::iota(std::begin(first_subsys), std::end(first_subsys), 0);

    // vector with labels of the second half of the qubits
    std::vector<idx> second_subsys(n);
    std::iota(std::begin(second_subsys), std::end(second_subsys), n);

    // QUANTUM STAGE
    // prepare the initial state |0>^\otimes n \otimes |0...01>
    ket psi = kron(st.zero(2 * n - 1), 1_ket);

    // apply Hadamards H^\otimes n on first half of the qubits
    for (idx i = 0; i < n; ++i) {
        psi = apply(psi, gt.H, {i});
    }

    // perform the modular exponentiation as a sequence of
    // modular multiplications
    for (idx i = 0; i < n; ++i) {
        // compute 2^(n-i-1) mod N
        bigint j = std::llround(std::pow(2, n - i - 1));
        // compute the a^(2^(n-i-1)) mod N
        bigint aj = modpow(a, static_cast<bigint>(j), N);
        // apply the controlled modular multiplication
        psi = applyCTRL(psi, gt.MODMUL(aj, N, n), {i}, second_subsys);
    }

    // apply inverse QFT on first half of the qubits
    psi = applyTFQ(psi, first_subsys);
    // END QUANTUM STAGE

    // FIRST MEASUREMENT STAGE
    auto measured1 = measure_seq(psi, first_subsys); // measure first n qubits
    std::vector<idx> vect_results1 = std::get<RES>(measured1); // results
    realT prob1 = prod(std::get<PROB>(measured1)); // probability of the result
    idx n1 = multiidx2n(vect_results1, std::vector<idx>(n, 2)); // binary to int
    auto x1 = static_cast<realT>(n1) / static_cast<realT>(D); // multiple of 1/r

    std::cout << ">> First measurement:  " << disp(vect_results1, " ") << " ";
    std::cout << "i.e., j = " << n1 << " with probability " << prob1;
    std::cout << '\n';

    bool failed = true;
    bigint r1 = 0, c1 = 0;
    for (auto&& elem : convergents(x1, 10)) {
        std::tie(c1, r1) = elem;
        auto c1r1 = static_cast<realT>(c1) / static_cast<realT>(r1);
        if (abs(x1 - c1r1) <
            1. / std::pow(2, (static_cast<realT>(n) - 1.) / 2.)) {
            failed = false;
            break;
        }
    }
    if (failed) {
        std::cout << ">> Factoring failed at stage 1, please try again!\n";
        std::exit(EXIT_FAILURE);
    }
    // END FIRST MEASUREMENT STAGE

    // SECOND MEASUREMENT STAGE
    auto measured2 = measure_seq(psi, first_subsys); // measure first n qubits
    std::vector<idx> vect_results2 = std::get<RES>(measured2); // results
    realT prob2 = prod(std::get<PROB>(measured2)); // probability of the result
    idx n2 = multiidx2n(vect_results2, std::vector<idx>(n, 2)); // binary to int
    auto x2 = static_cast<realT>(n2) / static_cast<realT>(D); // multiple of 1/r

    std::cout << ">> Second measurement: " << disp(vect_results2, " ") << " ";
    std::cout << "i.e., j = " << n2 << " with probability " << prob2;
    std::cout << '\n';

    failed = true;
    idx r2 = 0, c2 = 0;
    for (auto&& elem : convergents(x2, 10)) {
        std::tie(c2, r2) = elem;
        auto c2r2 = static_cast<realT>(c2) / static_cast<realT>(r2);
        if (abs(x2 - c2r2) <
            1. / std::pow(2, (static_cast<realT>(n) - 1.) / 2.)) {
            failed = false;
            break;
        }
    }
    if (failed) {
        std::cout << ">> Factoring failed at stage 2, please try again!\n";
        std::exit(EXIT_FAILURE);
    }
    // END SECOND MEASUREMENT STAGE

    // THIRD POST-PROCESSING STAGE
    idx r = lcm(static_cast<bigint>(r1),
                static_cast<bigint>(r2)); // candidate order of a mod N
    std::cout << ">> r = " << r << ", a^r mod N = "
              << modpow(a, static_cast<bigint>(r), static_cast<bigint>(N))
              << '\n';
    if (r % 2 == 0 && modpow(a, static_cast<bigint>(r / 2), N) !=
                          static_cast<bigint>(N - 1)) {
        // at least one of those is a non-trivial factor
        bigint p = gcd(modpow(a, static_cast<bigint>(r / 2), N) - 1, N);
        bigint q = gcd(modpow(a, static_cast<bigint>(r / 2), N) + 1, N);
        if (p == 1)
            p = N / q;
        if (q == 1)
            q = N / p;
        std::cout << ">> Factors: " << p << " " << q << '\n';
    } else {
        std::cout << ">> Factoring failed at stage 3, please try again!\n";
        std::exit(EXIT_FAILURE);
    }
    // END THIRD POST-PROCESSING STAGE
}
