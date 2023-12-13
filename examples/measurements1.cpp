// Measurements
// Source: ./examples/measurements1.cpp

#include <iostream>
#include <tuple>

#include "qpp/qpp.h"

int main() {
    using namespace qpp;

    ket psi = 00_ket;
    cmat U = gt.CNOT * kron(gt.H, gt.Id2);
    ket result = U * psi; // we have the Bell state (|00> + |11>) / sqrt(2)

    std::cout << ">> We just produced the Bell state:\n";
    std::cout << disp(result) << '\n';

    // apply a bit flip on the second qubit
    result = apply(result, gt.X, {1}); // we produced (|01> + |10>) / sqrt(2)
    std::cout << ">> We produced the Bell state:\n";
    std::cout << disp(result) << '\n';

    // measure the first qubit in the X basis
    auto measured = measure(result, gt.H, {0});
    std::cout << ">> Measurement result: " << std::get<RES>(measured) << '\n';
    std::cout << ">> Probabilities: ";
    std::cout << disp(std::get<PROB>(measured), ", ") << '\n';
    std::cout << ">> Resulting states:\n";
    for (auto&& it : std::get<ST>(measured))
        std::cout << disp(it) << "\n\n";
}
