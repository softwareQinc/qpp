// Source: ./examples/shor.cpp
//
// Implementation of Shor's quantum algorithm for integer factorization.
// Shor's algorithm finds the prime factors of an integer N by reducing the
// factorization problem to the problem of order-finding (finding the period
// 'r' of a function f(x) = a^x mod N).

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <optional>
#include <vector>

#include <qpp/qpp.hpp>

// performs continued fraction expansion to recover a candidate period 'r' from
// the measured phase value x = j/D
std::optional<qpp::bigint> get_candidate_period(qpp::realT x,
                                                qpp::realT threshold) {
    using namespace qpp;
    for (auto [numerator, denominator] : convergents(x, 10)) {
        // skip trivial results where the denominator is 1 or less
        if (denominator <= 1) {
            continue;
        }
        realT approximation =
            static_cast<realT>(numerator) / static_cast<realT>(denominator);

        // check if the convergent is within the allowed heuristic threshold
        if (std::abs(static_cast<long double>(x) - approximation) < threshold) {
            return denominator;
        }
    }

    return std::nullopt;
}

// structure to hold data returned from a quantum register measurement
struct MeasurementResult {
    std::vector<qpp::idx> bit_string; // raw measurement outcome (qubit states)
    qpp::idx integer_value; // base-10 integer represented by the bit string
    qpp::realT probability; // probability of observing this specific outcome
    qpp::bigint period_candidate; // period 'r' derived from the measurement
};

// performs a quantum measurement on the first register and attempts to
// extract a candidate period using continued fractions
std::optional<MeasurementResult>
perform_measurement(const qpp::ket& psi, const std::vector<qpp::idx>& subsys,
                    qpp::idx n, qpp::idx D, qpp::realT threshold) {
    using namespace qpp;
    // perform measurement on the specified qubits
    auto measurement_data = measure_seq(psi, subsys);
    std::vector<idx> bit_string = std::get<measure_idx::res>(measurement_data);
    idx integer_value = multiidx2n(bit_string, std::vector<idx>(n, 2));
    realT probability = prod(std::get<measure_idx::prob>(measurement_data));

    // calculate the phase x = j / 2^n
    realT x = static_cast<realT>(integer_value) / static_cast<realT>(D);

    // attempt to recover the period 'r'
    if (auto r = get_candidate_period(x, threshold)) {
        return MeasurementResult{bit_string, integer_value, probability,
                                 r.value()};
    }

    return std::nullopt;
}

// classical post-processing: attempts to extract non-trivial factors of N
// using the candidate period 'r'
std::optional<std::pair<qpp::bigint, qpp::bigint>>
find_factors(qpp::bigint a, qpp::bigint N, qpp::idx r) {
    using namespace qpp;
    // if r is odd, the algorithm cannot proceed
    if (r % 2 != 0) {
        return std::nullopt;
    }

    // check if a^(r/2) + 1 is a multiple of N (which yields a trivial factor)
    bigint val = modpow(a, static_cast<bigint>(r / 2), N);
    if (val == static_cast<bigint>(N - 1)) {
        return std::nullopt;
    }

    // use GCD to find candidate factors
    bigint p = gcd(val - 1, N);
    bigint q = gcd(val + 1, N);

    // ensure p and q are properly resolved if one resulted in a trivial GCD
    if (p == 1) {
        p = N / q;
    }
    if (q == 1) {
        q = N / p;
    }

    // verify the factors are non-trivial (1 < factor < N)
    if (p > 1 && p < N && q > 1 && q < N) {
        return std::make_pair(p, q);
    }

    return std::nullopt;
}

int main() {
    using namespace qpp;
    bigint N = 21;                   // the number to factor
    auto a = rand<bigint>(3, N - 1); // select a random 'a' co-prime with N
    while (gcd(a, N) != 1) {
        a = rand<bigint>(3, N - 1);
    }

    // register size: n is number of qubits, total qubits used is 2n
    // we need 2^n >= 2 * r^2 to guarantee precision for period finding
    auto n = static_cast<idx>(std::ceil(2 * std::log2(N)));
    auto D = idx{1} << n;

    // heuristic threshold for continued fraction convergence
    // std::pow is used instead of bit shifting to preserve fractional powers
    // when the exponent calculation results in a non-integer value
    auto threshold = 1. / std::pow(2, (static_cast<realT>(n) - 1.) / 2.);

    std::cout << ">> Factoring N = " << N << " with coprime a = " << a << '\n';
    std::cout << ">> Using 2*n = " << 2 * n << " qubits, 2^n = " << D
              << " and 2^(2n) = " << D * D << '\n';

    // map qubits to the first and second registers
    std::vector<idx> first_subsys(n), second_subsys(n);
    std::iota(first_subsys.begin(), first_subsys.end(), 0);
    std::iota(second_subsys.begin(), second_subsys.end(), n);

    // QUANTUM STAGE
    // initialize: First register |0...0> (size n), Second register |0...01>
    // (size n)
    ket psi = kron(st.zero((2 * n) - 1), 1_ket);

    // apply Hadamards to the first register to create a uniform superposition
    for (idx i = 0; i < n; ++i) {
        psi = apply(psi, gt.H, {i});
    }

    // modular exponentiation: perform controlled-U^j operations
    for (idx i = 0; i < n; ++i) {
        bigint j = idx{1} << (n - i - 1);
        bigint aj = modpow(a, j, N);

        // apply controlled modular multiplication (U^j)
        // NOTE: in a production circuit, this would be decomposed into base
        // gates.
        psi = applyCTRL(psi, gt.MODMUL(aj, N, n), {i}, second_subsys);
    }

    // apply Inverse Quantum Fourier Transform to the first register
    psi = applyTFQ(psi, first_subsys);

    // MEASUREMENT STAGE 1
    // measure the first register and attempt to derive the period 'r'
    auto m1 = perform_measurement(psi, first_subsys, n, D, threshold);
    if (!m1) {
        std::cout << ">> Factoring failed at stage 1, please try again!\n";
        std::exit(EXIT_FAILURE);
    }
    std::cout << ">> First measurement:  "
              << disp(m1->bit_string, IOManipContainerOpts{}.set_sep(" "))
              << " (j = " << m1->integer_value << ") with probability "
              << m1->probability << '\n';

    // MEASUREMENT STAGE 2
    // repeat measurement to acquire a second candidate period
    auto m2 = perform_measurement(psi, first_subsys, n, D, threshold);
    if (!m2) {
        std::cout << ">> Factoring failed at stage 2, please try again!\n";
        std::exit(EXIT_FAILURE);
    }
    std::cout << ">> Second measurement: "
              << disp(m2->bit_string, IOManipContainerOpts{}.set_sep(" "))
              << " (j = " << m2->integer_value << ") with probability "
              << m2->probability << '\n';

    // POST-PROCESSING STAGE 3
    // combine periods using Least Common Multiple (LCM) to improve success
    // probability
    bigint r1 = m1->period_candidate;
    bigint r2 = m2->period_candidate;
    idx r = lcm(r1, r2);

    std::cout << ">> r = " << r
              << ", a^r mod N = " << modpow(a, static_cast<bigint>(r), N)
              << '\n';

    // use the period to identify the factors
    if (auto factors = find_factors(a, N, r)) {
        auto [p, q] = factors.value();
        std::cout << ">> Factors: " << p << " " << q << '\n';
    } else {
        std::cout << ">> Factoring failed at stage 3, please try again!\n";
        std::exit(EXIT_FAILURE);
    }
}
