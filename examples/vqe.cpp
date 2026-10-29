// Source: ./examples/vqe.cpp
//
// Variational Quantum Eigensolver (VQE)

#include <cmath>
#include <iomanip>
#include <iostream>

#include <qpp/qpp.hpp>

int main() {
    using namespace qpp;

    std::cout << ">> VQE example\n";

    // Hamiltonian is (Z x I) + (I x Z) + (X x X)
    const cmat H = kron(gt.Z, gt.Id2) + kron(gt.Id2, gt.Z) + kron(gt.X, gt.X);

    // quantum part: prepare the trial state and evaluate its energy.
    const auto energy = [&H](realT theta) -> realT {
        QCircuit qc{2};
        qc.gate(gt.RY(theta), 0);
        qc.CTRL(gt.X, 0, 1);

        QEngine engine{qc}; // starts in |00>
        engine.execute();

        const ket psi = engine.get_state();
        return expect(H, psi);
    };

    // classical part: gradient descent using the parameter-shift rule.
    const realT pi = std::acos(-1.0);
    const realT learning_rate = 0.2;
    realT theta = 0.3;

    std::cout << std::fixed << std::setprecision(10);

    for (idx iteration = 0; iteration < 200; ++iteration) {
        const realT E = energy(theta);
        const realT gradient =
            0.5 * (energy(theta + (pi / 2)) - energy(theta - (pi / 2)));

        if (iteration % 5 == 0) {
            std::cout << "    iteration = " << std::setw(3) << iteration
                      << "  energy = " << E << "  theta = " << theta
                      << std::endl;
        }

        if (std::abs(gradient) < 1e-10) {
            break;
        }

        theta -= learning_rate * gradient;
    }

    std::cout << ">> VQE energy:   " << energy(theta)
              << "\n>> Exact energy: " << -std::sqrt(5.0) << "\n";
}
