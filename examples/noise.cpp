// Quantum noise
// Source: ./examples/noise.cpp
#include <iostream>

#include "qpp.h"

int main() {
    using namespace qpp;

    std::cout << ">> Depolarizing qubit noise acting on a random state\n";
    realT p = 0.75;      // depolarizing probability (fully depolarizing)
    idx N = 10000;       // number of trials
    ket psi = randket(); // initial state
    cmat rho = cmat::Zero(2, 2);        // final density matrix
    QubitDepolarizingNoise noise{0.75}; // constructs a noise instance
    // compute the noise action; see qpp::NoiseBase::operator() for other
    // overloads for multi-partite states
    for (idx i = 0; i < N; ++i) {
        // apply the noise
        rho += prj(noise(psi));
    }
    rho /= static_cast<realT>(N); // normalize the resulting density matrix
    std::cout << ">> Resulting state after " << N
              << " actions of depolarizing noise with p = " << p << ":\n";
    std::cout << disp(rho) << '\n';
}
