// Measurements
// Source: ./examples/measurements2.cpp

#include <iostream>
#include <tuple>
#include <vector>

#include "qpp/qpp.hpp"

int main() {
    using namespace qpp;

    ket psi = st.b00;
    std::vector<idx> subsys = {0};
    idx result;
    std::vector<realT> probs;
    std::vector<ket> states;

    // measures the first subsystem of the Bell state (|00> + |11>) / sqrt(2)
    // in the X basis
    std::tie(result, probs, states) = measure(psi, gt.H, {0});
    std::cout << ">> Measuring part "
              << disp(subsys, IOManipContainerOpts{}.set_sep(" "))
              << " of the state:\n";
    std::cout << disp(psi) << '\n';
    std::cout << ">> Measurement result: " << result << '\n';
    std::cout << ">> Probabilities: "
              << disp(probs, IOManipContainerOpts{}.set_sep(", ")) << '\n';
    std::cout << ">> Resulting normalized post-measurement states:\n";

    for (auto&& it : states) {
        std::cout << disp(it) << "\n\n";
    }

    // measure 2 subsystems out of a 4-qubit random density matrix
    cmat rho = randrho(16);
    subsys = {1, 2};
    cmat U = randU(4); // random basis on 2 qubits

    std::cout << ">> Measuring qubits "
              << disp(subsys, IOManipContainerOpts{}.set_sep(" "))
              << " of a 4-qubit random state in the random basis:\n";
    std::cout << disp(U) << '\n';

    std::vector<cmat> density_operators;

    std::tie(result, probs, density_operators) = measure(rho, U, {1, 2});
    std::cout << ">> Measurement result: " << result << '\n';
    std::cout << ">> Probabilities: "
              << disp(probs, IOManipContainerOpts{}.set_sep(", ")) << '\n';
    std::cout << ">> Sum of the probabilities: "
              << sum(probs.begin(), probs.end()) << '\n';
    std::cout << ">> Resulting normalized post-measurement states:\n";

    for (auto&& it : density_operators) {
        std::cout << disp(it) << "\n\n";
    }

    // check now how the state after the measurement "looks"
    // on the leftover subsystems {0, 3}

    // it should be the same as the partial trace over {1, 2}
    // of the initial state (before the measurement), as local CPTP maps
    // do not influence the complementary subsystems

    cmat rho_bar = ptrace(rho, subsys);
    cmat rho_out_bar = cmat::Zero(4, 4);

    // compute the resulting mixed state after the measurement
    for (idx i = 0; i < static_cast<idx>(probs.size()); ++i) {
        rho_out_bar += probs[i] * density_operators[i];
    }

    // verification
    std::cout << ">> Norm difference: " << norm(rho_bar - rho_out_bar) << '\n';

    // random Kraus
    std::cout << ">> Random channel on part of the state\n";
    rho_bar = ptrace(rho, {1, 3});
    rho_out_bar = ptrace(apply(rho, randkraus(3, 4), {1, 3}), {1, 3});

    // verification
    std::cout << ">> Norm difference: " << norm(rho_bar - rho_out_bar) << '\n';

    std::cout << ">> Sequential measurements on the state/density matrix:\n";
    psi = 0.8 * 01_ket + 0.6 * 10_ket;
    rho = prj(psi); // psi * adjoint(psi)
    std::cout << disp(psi) << '\n';

    std::vector<idx> subsys_ket{0};
    std::vector<idx> subsys_rho{1};

    auto [ms_ket, probs_ket, state_ket] = measure_seq(psi, subsys_ket);
    auto [ms_rho, probs_rho, state_rho] = measure_seq(rho, subsys_rho);

    // ket
    std::cout << ">> Ket, measuring subsystem(s) ";
    std::cout << disp(subsys_ket, IOManipContainerOpts{}.set_sep(" ")) << '\n';
    std::cout << ">> Outcome(s): "
              << disp(ms_ket, IOManipContainerOpts{}.set_sep(" ")) << '\n';
    std::cout << ">> Probability:  " << prod(probs_ket) << '\n';
    std::cout << ">> Resulting state:\n";
    std::cout << disp(state_ket) << '\n';

    // density matrix
    std::cout << ">> Density matrix, measuring subsystem(s) ";
    std::cout << disp(subsys_rho, IOManipContainerOpts{}.set_sep(" ")) << '\n';
    std::cout << ">> Outcome(s): "
              << disp(ms_rho, IOManipContainerOpts{}.set_sep(" ")) << '\n';
    std::cout << ">> Probability:  " << prod(probs_rho) << '\n';
    std::cout << ">> Resulting state:\n";
    std::cout << disp(state_rho) << '\n';
}
