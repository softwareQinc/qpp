// Reversible classical circuits
// Source: ./examples/reversible2.cpp
#include <algorithm>
#include <iostream>
#include <random>
#include <vector>

#include "qpp.h"

int main() {
    using namespace qpp;
    std::cout << ">> Classical reversible circuits\n";

    const idx N = 32;          // number of N
    const idx num_trials = 10; // number of num_trials
    Bit_circuit bit_circuit{N};
    bit_circuit.rand(); // randomize the vector
    Bit_circuit initial_bit_circuit = bit_circuit;

    std::cout << ">> Initial randomized bit circuit\n\t" << initial_bit_circuit
              << '\n';

    // randomize the indices where Toffoli will be applied
    std::random_device rd;
    std::mt19937 gen{rd()};
    std::vector<std::vector<idx>> indices(num_trials);

    // generate the indices
    for (idx i = 0; i < num_trials; ++i) {
        std::vector<idx> v(N);
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
        bit_circuit.TOF(indices[i]);
    }
    std::cout << ">> Initial/Intermediate bit circuit\n\t";
    std::cout << initial_bit_circuit << "\n\t" << bit_circuit << '\n';
    std::cout << ">> Hamming weight (intermediate circuit): "
              << bit_circuit.count() << '\n';
    std::cout << ">> Hamming distance (from the initial circuit): "
              << initial_bit_circuit - bit_circuit << '\n';
    std::cout << ">> NOT count (intermediate circuit, should be zero): "
              << bit_circuit.gate_count.NOT << '\n';
    std::cout << ">> X count (should be same as NOT count, i.e. zero): "
              << bit_circuit.gate_count.X << '\n';
    std::cout << ">> Toffoli count (intermediate circuit): "
              << bit_circuit.gate_count.TOF << '\n';

    // apply again the same Toffoli gates in reverse order
    std::cout << ">> Applying again Toffoli gates to bits\n";
    for (idx i = num_trials; i-- > 0;) {
        for (auto&& elem : indices[i])
            std::cout << '\t' << elem << " ";
        std::cout << '\n';
        bit_circuit.TOF(indices[i]);
    }

    std::cout << ">> Final bit circuit:\n\t" << bit_circuit << '\n';
    std::cout << ">> Hamming weight: " << bit_circuit.count() << '\n';

    std::cout << std::boolalpha;
    std::cout << ">> Are the initial and final circuits equal? "
              << (initial_bit_circuit == bit_circuit) << '\n';

    std::cout << ">> Changing the string representation of the bits:\n\t";
    std::cout << bit_circuit.to_string('o', 'i') << '\n';
}
