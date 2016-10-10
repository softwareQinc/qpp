// Measurements, second take
// Source: ./examples/measurements2.cpp
#include <iostream>
#include <vector>
#include <tuple>
#include "qpp.h"

using namespace qpp;

int main()
{
    ket psi = st.b00;
    std::vector<idx> subsys = {0};
    idx result;
    std::vector<double> probs;
    std::vector<cmat> states;

    // measures the first subsystem of the Bell state (|00> + |11>) / sqrt(2)
    // in the X basis
    std::tie(result, probs, states) = measure(psi, gt.H, {0});
    std::cout << ">> Measuring part " << disp(subsys, " ")
              << " of the state: " << std::endl;
    std::cout << disp(psi) << std::endl;
    std::cout << ">> Measurement result: " << result << std::endl;
    std::cout << ">> Probabilities: " << disp(probs, ", ") << std::endl;
    std::cout << ">> Resulting normalized post-measurement states: "
              << std::endl;

    for (auto&& it: states)
        std::cout << disp(it) << std::endl << std::endl;

    // measure 2 subsystems out of a 4-qubit random density matrix
    cmat rho = randrho(16);
    subsys = {1, 2};
    cmat U = randU(4); // random basis on 2 qubits

    std::cout << ">> Measuring qubits " << disp(subsys, " ")
              << " of a 4-qubit random state in the random basis:" << std::endl;
    std::cout << disp(U) << std::endl;

    std::tie(result, probs, states) = measure(rho, U, {1, 2});
    std::cout << ">> Measurement result: " << result << std::endl;
    std::cout << ">> Probabilities: " << disp(probs, ", ") << std::endl;
    std::cout << ">> Sum of the probabilities: "
              << sum(probs.begin(), probs.end()) << std::endl;
    std::cout << ">> Resulting normalized post-measurement states: "
              << std::endl;

    for (auto&& it: states)
        std::cout << disp(it) << std::endl << std::endl;

    // check now how the state after the measurement "looks"
    // on the left over subsystems {0, 3}

    // it should be the same as the partial trace over {1, 2}
    // of the initial state (before the measurement), as local CPTP maps
    // do not influence the complementary subsystems

    cmat rho_bar = ptrace(rho, subsys);
    cmat rho_out_bar = cmat::Zero(4, 4);

    // compute the resulting mixed state after the measurement
    for (idx i = 0; i < probs.size(); ++i)
        rho_out_bar += probs[i] * states[i];

    // verification
    std::cout << ">> Norm difference: " << norm(rho_bar - rho_out_bar)
              << std::endl;

    // random Kraus
    std::cout << ">> Random channel on part of the state " << std::endl;
    rho_bar = ptrace(rho, {1, 3});
    rho_out_bar = ptrace(apply(rho, randkraus(3, 4), {1, 3}), {1, 3});

    // verification
    std::cout << ">> Norm difference: " << norm(rho_bar - rho_out_bar)
              << std::endl;

    std::cout << ">> Sequential measurements on the state/density matrix:"
              << std::endl;
    psi = 0.8 * mket({0, 1}) + 0.6 * mket({1, 0});
    rho = psi * adjoint(psi);
    std::cout << disp(psi) << std::endl;

    std::vector<idx> subsys_ket{0};
    std::vector<idx> subsys_rho{1};

    auto measured_ket = measure_seq(psi, subsys_ket);
    auto measured_rho = measure_seq(rho, subsys_rho);

    // ket
    std::cout << ">> Ket, measuring subsystem(s) ";
    std::cout << disp(subsys_ket, " ") << std::endl;
    std::cout << ">> Outcome(s): " << disp(std::get<0>(measured_ket), " ")
              << std::endl;
    std::cout << ">> Probability:  " << std::get<1>(measured_ket) << std::endl;
    std::cout << ">> Resulting state:  " << std::endl;
    std::cout << disp(std::get<2>(measured_ket)) << std::endl;

    // density matrix
    std::cout << ">> Density matrix, measuring subsystem(s) ";
    std::cout << disp(subsys_rho, " ") << std::endl;
    std::cout << ">> Outcome(s): " << disp(std::get<0>(measured_rho), " ")
              << std::endl;
    std::cout << ">> Probability:  " << std::get<1>(measured_rho) << std::endl;
    std::cout << ">> Resulting state:  " << std::endl;
    std::cout << disp(std::get<2>(measured_rho)) << std::endl;
}
