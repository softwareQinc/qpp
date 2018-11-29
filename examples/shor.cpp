// Shor's algorithm
// Source: ./examples/shor.cpp
#include <cmath>
#include <iostream>
#include <vector>

// in progress...

#include "qpp.h"

int main() {
    using namespace qpp;
    idx N = 21;       // number to factor
    idx a = 5;        // co-prime factor
    idx p = 1, q = 1; // factors
    std::cout << ">> Shor's algorithm, factoring N = " << N << '\n';

    // number of bits required to represent N
    idx n = static_cast<idx>(std::ceil(std::log2(N)));

    // first half of the qubits
    std::vector<idx> subsys0(n);
    std::iota(std::begin(subsys0), std::end(subsys0), 0);

    // second half of the qubits
    std::vector<idx> subsys1(n);
    std::iota(std::begin(subsys1), std::end(subsys1), n);

    // prepare the initial state |0>^\otimes n \otimes |0...01>
    ket psi = kron(st.zero(2 * n - 1), 1_ket);

    // apply QFT on first half of the qubits
    psi = applyQFT(psi, subsys0);

    // perform the modular exponentiation
    for (idx i = 0; i < n; ++i) {
        idx j = static_cast<idx>(std::llround(std::pow(2, n - i - 1)));
        idx aj = modpow(a, j, N);
        psi = applyCTRL(psi, gt.MODMUL(aj, N), {i}, subsys1);
    }

    // apply inverse QFT on first half of the qubits
    psi = applyINVQFT(psi, subsys0);

    // measure the first half of the qubits
    auto measured = measure_seq(psi, subsys0);
    auto m = std::get<0>(measured);
    idx x = multiidx2n(m, std::vector<idx>(n, 2));
    auto prob = std::get<1>(measured);
    std::cout << disp(m, " ") << '\n';
    std::cout << x << '\n';
    std::cout << prob << '\n';

    auto cfrac = x2contfrac(static_cast<double>(x) / std::pow(2, n), 10);
    std::cout << disp(cfrac, ", ") << "\n";

    std::cout << ">> The factors of " << N << " likely are " << p << " and "
              << q << '\n';
}
