// Reversible classical circuits
// Source: ./examples/reversible.cpp
#include <iostream>

#include "qpp.h"

int main() {
    using namespace qpp;
    std::cout << ">> Classical reversible circuits. ";
    std::cout << "Bits are labeled from right to left,\n   ";
    std::cout << "i.e. bit zero is the least significant bit (rightmost).\n";

    Dynamic_bitset bits{4};                                // 4 classical bits
    std::cout << ">> Initial bitset:\n\t" << bits << '\n'; // display them

    bits.rand(); // randomize the bits
    std::cout << ">> After randomization:\n\t" << bits << '\n'; // display them

    Bit_circuit bit_circuit{bits}; // bit circuit

    std::cout << ">> Apply X_0, followed by CNOT_02, CNOT_13 and TOF_013\n";
    bit_circuit.X(0); // apply a NOT gate on first bit
    bit_circuit.CNOT({0, 2}).CNOT({1, 3}).TOF({0, 1, 3}); // sequence operations

    std::cout << ">> Final bit circuit:\n\t" << bit_circuit << '\n';
    std::cout << ">> 3rd bit: " << bit_circuit.get(2) << '\n';
    std::cout << ">> CNOT count: " << bit_circuit.gate_count.CNOT << '\n';

    bit_circuit.reset(); // resets the circuit
    std::cout << ">> Reseted circuit:\n\t" << bit_circuit << '\n';
    std::cout << ">> CNOT count: " << bit_circuit.gate_count.CNOT << '\n';
}
