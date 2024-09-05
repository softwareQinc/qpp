// Classical reversible circuits
// Source: ./examples/reversible2.cpp

#include <algorithm>
#include <iostream>
#include <random>
#include <vector>

#include "qpp/qpp.hpp"

int main() {
    using namespace qpp;

    std::cout << ">> Classical reversible circuits (2)\n";

    idx n = 32;          // number of bits
    idx num_trials = 10; // number of trials
    Bit_circuit bc{n};
    bc.rand(); // randomize the vector
    Bit_circuit bc_initial = bc;

    std::cout << ">> Initial randomized bit circuit state\n"
              << bc_initial.to_string() << '\n';

    // randomize the indices where Toffoli will be applied
    std::random_device rd;
    std::mt19937 gen{rd()};
    std::vector<std::vector<idx>> indices(num_trials);

    // generate the indices
    for (idx i = 0; i < num_trials; ++i) {
        std::vector<idx> v(n);
        std::iota(v.begin(), v.end(), 0);
        std::shuffle(v.begin(), v.end(), gen);
        std::vector<idx> tof(v.data(), v.data() + 3);
        indices[i] = tof;
    }

    // apply the Toffoli gate to random places
    std::cout << ">> Applying Toffoli gates to bits\n";
    for (idx i = 0; i < num_trials; ++i) {
        for (auto&& elem : indices[i]) {
            std::cout << '\t' << elem << " ";
        }
        std::cout << '\n';
        bc.TOF(indices[i][0], indices[i][1], indices[i][2]);
    }
    std::cout << ">> Intermediate bit circuit state\n";
    std::cout << bc.to_string() << '\n';
    std::cout << ">> Hamming distance (from the initial circuit): "
              << bc_initial - bc << '\n';
    std::cout << ">> NOT count (intermediate circuit, should be zero): "
              << bc.get_gate_count("NOT") << '\n';
    std::cout << ">> X count (should be same as NOT count, i.e., zero): "
              << bc.get_gate_count("X") << '\n';
    std::cout << ">> Toffoli count (intermediate circuit): "
              << bc.get_gate_count("TOF") << '\n';

    // apply again the same Toffoli gates in reverse order
    std::cout << ">> Applying again Toffoli gates to bits\n";
    for (idx i = num_trials; i-- > 0;) {
        for (auto&& elem : indices[i]) {
            std::cout << '\t' << elem << " ";
        }
        std::cout << '\n';
        bc.TOF(indices[i][0], indices[i][1], indices[i][2]);
    }

    std::cout << ">> Final bit circuit state\n" << bc.to_string() << '\n';
    std::cout << ">> Are the initial and final circuits equal? "
              << std::boolalpha << (bc_initial == bc) << std::noboolalpha
              << '\n';

    std::cout << ">> Changing the string representation of the bits\n";
    std::cout << bc.to_string('o', 'i') << '\n';
}
