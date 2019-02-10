// Used for testing, do not use it as an example
#include <iostream>

#include "qpp.h"
using gate_hash_t = std::unordered_set<qpp::cmat, qpp::internal::HashEigen,
                                       qpp::internal::KeyEqualEigen>;

void add_gate(gate_hash_t& gate_set, const qpp::cmat& U) {
    auto ret = gate_set.insert(U);
    if (ret.second == false)
        std::cout << "gate exists\n";
    else
        std::cout << "new gate\n";
}

int main() {
    using namespace qpp;
    /////////// testing ///////////

    gate_hash_t gate_set;

    cmat x(1, 1);
    x << 1.;

    cmat y(1, 1);
    y << 1.000000000000001;

    std::cout << hash_eigen_expression(gt.CNOT) << '\n';
    std::cout << hash_eigen_expression(gt.TOF) << '\n';
    std::cout << hash_eigen_expression(x) << '\n';

    internal::HashEigen h;
    std::cout << h(gt.X + gt.X) << '\n';
    std::cout << h(gt.X) << '\n';
    std::cout << h(gt.Z) << '\n';

    add_gate(gate_set, x);
    add_gate(gate_set, x);
    add_gate(gate_set, gt.CNOT);
    add_gate(gate_set, gt.TOF);
    add_gate(gate_set, st.b00);
    add_gate(gate_set, st.b00);
    add_gate(gate_set, x);
    add_gate(gate_set, gt.X);
    add_gate(gate_set, gt.Z);
    add_gate(gate_set, gt.Z);

    std::cout << *(gate_set.find(x)) << '\n';
    std::cout << *(gate_set.find(st.b00)) << '\n';
    std::cout << *(gate_set.find(gt.X)) << '\n';
    std::cout << *(gate_set.find(gt.Z)) << '\n';
    std::cout << (x == y) << '\n';

    cmat c(2, 2);
    dmat d(2, 2);
    c << 2, 1, 3, 4;
    d << 1, 2, 3, 4;
    std::cout << h(c) << '\n';
    std::cout << h(d) << '\n';
}
