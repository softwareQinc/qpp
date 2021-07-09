// BB84
// Source: ./examples/BB84.cpp
#include <iostream>
#include <tuple>
#include <utility>
#include <vector>

#include "qpp.h"

// (basis, state) pair collection type ; 0 -> Z basis, 1 -> X basis
using bases_states_T = std::vector<std::pair<short, short>>;
using key_T = std::vector<short>; // key type (bit string)

// display Alice's and Bob's bases choices and states
void display(const bases_states_T& Alice_bases_states,
             const bases_states_T& Bob_bases_states) {
    auto n = Alice_bases_states.size();
    std::cout << "Alice's states: ";
    for (std::size_t i = 0; i < n; ++i)
        std::cout << Alice_bases_states[i].second << ' ';
    std::cout << '\n';
    std::cout << "Alice's bases:  ";
    for (std::size_t i = 0; i < n; ++i)
        std::cout << (Alice_bases_states[i].first == 0 ? 'Z' : 'X') << ' ';
    std::cout << '\n';
    std::cout << "Bob's bases:    ";
    for (std::size_t i = 0; i < n; ++i)
        std::cout << (Bob_bases_states[i].first == 0 ? 'Z' : 'X') << ' ';
    std::cout << '\n';
    std::cout << "Bob's states:   ";
    for (std::size_t i = 0; i < n; ++i)
        std::cout << Bob_bases_states[i].second << ' ';
    std::cout << '\n';
}

// key sifting, removes the locations where the bases do not coincide
void key_sift(bases_states_T& Alice_bases_states,
              bases_states_T& Bob_bases_states) {
    auto n = Alice_bases_states.size();
    bases_states_T result_A, result_B;
    for (std::size_t i = 0; i < n; ++i) {
        if (Alice_bases_states[i].first != Bob_bases_states[i].first)
            continue;
        result_A.emplace_back(Alice_bases_states[i]);
        result_B.emplace_back(Bob_bases_states[i]);
    }
    Alice_bases_states = result_A;
    Bob_bases_states = result_B;
}

// retrieves the raw key
key_T raw_key(const bases_states_T& bases_states) {
    auto n = bases_states.size();
    key_T result(n);
    for (std::size_t i = 0; i < n; ++i)
        result[i] = bases_states[i].second;
    return result;
}

// computes the final key; technically, here is where we do error correction and
// privacy amplification; however, here we simply discard the bits that differ
key_T final_key(const key_T& Alice_raw_key, const key_T& Bob_raw_key) {
    auto n = Alice_raw_key.size();
    key_T result;
    for (std::size_t i = 0; i < n; ++i) {
        if (Alice_raw_key[i] != Bob_raw_key[i])
            continue;
        result.emplace_back(Alice_raw_key[i]);
    }
    return result;
}

int main() {
    using namespace qpp;

    std::size_t n = 100; // number of qubits Alice sends to Bob
    double p = 0.5; // probability of Eve intercepting (and altering) the qubits
    std::cout << ">> BB84, sending n = " << n << " qubits from Alice to Bob\n";
    std::cout << ">> With probability p = " << p
              << " Eve intercepts the qubits and randomly measures them in the "
                 "Z or X basis, then sends them to Bob\n";

    // Alice's basis and state pairs (basis, state); 0 -> Z basis, 1 -> X basis
    bases_states_T Alice_bases_states(n);
    for (auto& elem : Alice_bases_states) {
        // chose a random basis, 0 -> Z basis, 1 -> X basis
        short basis = bernoulli() ? 0 : 1;
        // chose a random state, |0> or |1> in the basis 'basis'
        short state = bernoulli() ? 0 : 1;
        elem = std::make_pair(basis, state);
    }

    // Bob's basis and state pairs (basis, state); 0 -> Z basis, 1 -> X basis
    bases_states_T Bob_bases_states(n);
    for (auto& elem : Bob_bases_states) {
        // chose a random basis, 0 -> Z basis, 1 -> X basis
        short basis = bernoulli() ? 0 : 1;
        elem.first = basis;
    }

    // Alice "prepares" the qubits and "sends" them to Bob one by one
    // Eve is in the middle and intercepts/resends the qubits with probability p
    // Bob measures the received qubits one by one
    for (std::size_t i = 0; i < n; ++i) {
        auto basis_A = Alice_bases_states[i].first;
        auto state_A = Alice_bases_states[i].second;
        ket psi = (state_A == 0) ? 0_ket : 1_ket;
        if (basis_A != 0) // if X basis
            psi = gt.H * psi;

        // Eve intercepts the qubit and randomly measures it in the Z or X
        // basis, then sends it to Bob
        if (bernoulli(p)) {
            // chose a random basis, 0 -> Z basis, 1 -> X basis
            short basis_E = bernoulli();
            cmat U_E = (basis_E == 0) ? gt.Z : gt.H;
            auto measure_E = measure(psi, U_E);
            auto m = std::get<RES>(measure_E); // measurement result
            psi = std::get<ST>(measure_E)[m];  // update the state accordingly
        }

        // Bob measures the qubit Eve re-sent
        auto basis_B = Bob_bases_states[i].first;
        // Bob's measurement eigenvectors
        cmat U_B = (basis_B == 0) ? gt.Z : gt.H;
        auto measure_B = measure(psi, U_B);
        Bob_bases_states[i].second = std::get<RES>(measure_B);
    }

    // display the results before bases sifting
    std::cout << ">> Before sifting\n";
    display(Alice_bases_states, Bob_bases_states);

    // filter on same bases
    key_sift(Alice_bases_states, Bob_bases_states);

    // display the results after basis sifting
    std::cout << ">> After sifting (raw key)\n";
    display(Alice_bases_states, Bob_bases_states);

    std::cout << ">> Established keys\n";
    // display the raw key on Alice's side
    auto raw_key_A = raw_key(Alice_bases_states);
    std::cout << "Alice's raw key: " << disp(raw_key_A, " ", "", "") << '\n';

    // display the raw key on Bob's side
    auto raw_key_B = raw_key(Bob_bases_states);
    std::cout << "Bob's raw key:   " << disp(raw_key_B, " ", "", "") << '\n';

    // display the final key and the corresponding rate
    auto key = final_key(raw_key_A, raw_key_B);
    std::cout << "Final key:       " << disp(key, " ", "", "") << '\n';

    double key_rate = static_cast<double>(key.size()) / static_cast<double>(n);

    std::cout << ">> Bits/keys sizes: " << n << '/' << raw_key_A.size() << '/'
              << key.size() << '\n';
    std::cout << ">> Final key rate: " << key_rate << '\n';
}
