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
        // hash collision
        if ((*it).cols() == U.cols() && (*it).rows() == U.rows()) {
            if ((*it) != U) {
                throw 42;
            }
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

    cmat x(1, 2);
    x << 1.,1;

    cmat y(1, 1);
    y << 2.000000000000001,1;

    std::cout << hash_eigen_expression(gt.CNOT) << '\n';
    std::cout << hash_eigen_expression(gt.TOF) << '\n';

    // std::cout << (x != y) << '\n';
    add_gate(gate_set, gt.Id2);
    add_gate(gate_set, gt.Z);
    add_gate(gate_set, gt.X);
    add_gate(gate_set, gt.X);
    add_gate(gate_set, gt.CNOT);
    add_gate(gate_set, gt.TOF);
    add_gate(gate_set, gt.Id2);

    add_gate(gate_set, x);
//    add_gate(gate_set, y);

    std::cout << *(gate_set.find(gt.x)) << '\n';



}
