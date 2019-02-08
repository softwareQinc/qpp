// Used for testing, do not use it as an example
#include <iostream>
#include <unordered_set>

#include "qpp.h"

void add_gate(std::unordered_set<qpp::cmat>& gate_set, const qpp::cmat& U) {
    auto it = gate_set.find(U);

    // gate does not exist
    if (it == gate_set.end()) {
        std::cout << "new gate\n";
        gate_set.insert(U);
    }
    // gate exists already
    else {
        //std::cout << "here\n";
        // hash collision
        if ((*it) != U) {
            throw 42;
        }
        // no hash collision
        else {
            std::cout << "gate exists\n";
            return;
        }
    }
}

int main() {
    using namespace qpp;
    /////////// testing ///////////

    std::unordered_set<cmat> gate_set;

    cmat x(1, 1);
    x << 1.;

    std::cout << hash_eigen_expression(gt.CNOT) << '\n';
    std::cout << hash_eigen_expression(gt.TOF) << '\n';
    std::cout << hash_eigen_expression(x) << '\n';

    add_gate(gate_set, x);
    std::cout << "\t END1\n";
    add_gate(gate_set, x);
    std::cout << "\t END2\n";
    add_gate(gate_set, gt.CNOT);
    std::cout << "\t END3\n";
    add_gate(gate_set, gt.TOF);
    std::cout << "\t END4\n";
    add_gate(gate_set, st.b00);
    std::cout << "\t END5\n";
    add_gate(gate_set, st.b00);
    std::cout << "\t END6\n";

    std::cout << *(gate_set.find(x)) << '\n';
    std::cout << "\t END7\n";
}
