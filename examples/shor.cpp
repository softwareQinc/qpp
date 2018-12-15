// Shor's algorithm
// Source: ./examples/shor.cpp
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>

#include "qpp.h"

int main() {
    using namespace qpp;
    idx N = 35; // number to factor
    idx a = 4;  // co-prime with N, we need an even order!
    // number of bits required to represent N
    // idx n = static_cast<idx>(std::ceil(2 * std::log2(N)));
    idx n = 7; // suffices for this case, as the order 'r' of 'a' 'mod 35' is 6
               // and 2^n >= 2 * r^2
    idx D = static_cast<idx>(std::llround(std::pow(2, n))); // dimension 2^n

    std::cout << ">> Factoring N = " << N << " with coprime a = " << a << '\n';
    std::cout << ">> n = " << n << " and 2^n = " << D << '\n';

    // vector with labels of the first half of the qubits
    std::vector<idx> first_subsys(n);
    std::iota(std::begin(first_subsys), std::end(first_subsys), 0);

    // vector with labels of the second half of the qubits
    std::vector<idx> second_subsys(n);
    std::iota(std::begin(second_subsys), std::end(second_subsys), n);

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
        idx j = static_cast<idx>(std::llround(std::pow(2, n - i - 1)));
        // compute the a^(2^(n-i-1)) mod N
        idx aj = modpow(a, j, N);
        // apply the controlled modular multiplication
        psi = applyCTRL(psi, gt.MODMUL(aj, N, n), {i}, second_subsys);
    }

    // apply inverse QFT on first half of the qubits
    psi = applyINVQFT(psi, first_subsys);

    // FIRST STAGE
    auto measured1 = measure_seq(psi, first_subsys); // measure first n qubits
    auto list_results1 = std::get<0>(measured1);     // measurement results
    auto prob1 = std::get<1>(measured1); // probability of the result
    idx n1 = multiidx2n(list_results1, std::vector<idx>(n, 2)); // binary to int
    double x1 = static_cast<double>(n1) / D;

    std::cout << ">> First measurement:  " << disp(list_results1, " ") << " ";
    std::cout << "i.e. j = " << n1 << " with probability " << prob1;
    std::cout << '\n';

    bool failed = true;
    idx r1, c1;
    for (auto&& elem : convergents(x1, 10)) {
        std::tie(c1, r1) = elem;
        double c1r1 = static_cast<double>(c1) / r1;
        if (abs(x1 - c1r1) < 1. / std::pow(2, (n - 1) / 2.)) {
            failed = false;
            break;
        }
    }
    if (failed) {
        std::cout << ">> Factoring failed at stage 1, please try again!\n";
        std::exit(EXIT_SUCCESS);
    }
    // END FIRST STAGE

    // SECOND STAGE
    auto measured2 = measure_seq(psi, first_subsys); // measure first n qubits
    auto list_results2 = std::get<0>(measured2);     // measurement results
    auto prob2 = std::get<1>(measured2); // probability of the result
    idx n2 = multiidx2n(list_results2, std::vector<idx>(n, 2)); // binary to int
    double x2 = static_cast<double>(n2) / D;

    std::cout << ">> Second measurement: " << disp(list_results2, " ") << " ";
    std::cout << "i.e. j = " << n2 << " with probability " << prob2;
    std::cout << '\n';

    failed = true;
    idx r2, c2;
    for (auto&& elem : convergents(x2, 10)) {
        std::tie(c2, r2) = elem;
        double c2r2 = static_cast<double>(c2) / r2;
        if (abs(x2 - c2r2) < 1. / std::pow(2, (n - 1) / 2.)) {
            failed = false;
            break;
        }
    }
    if (failed) {
        std::cout << ">> Factoring failed at stage 2, please try again!\n";
        std::exit(EXIT_SUCCESS);
    }
    // END SECOND STAGE

    // THIRD STAGE, POST-PROCESSING
    idx r = lcm(r1, r2);
    if (r % 2 == 0 && modpow(a, r / 2, N) != static_cast<bigint>(N - 1)) {
        std::cout << ">> We (possibly) found the order of a (mod N): ";
        std::cout << r << '\n';
        std::cout << ">> Possible factors: ";
        std::cout << gcd(modpow(a, r / 2, N) - 1, N) << " ";
        std::cout << gcd(modpow(a, r / 2, N) + 1, N) << '\n';
    } else {
        std::cout << ">> Factoring failed at stage 3, please try again!\n";
        std::exit(EXIT_SUCCESS);
    }
    // END THIRD STAGE
}
