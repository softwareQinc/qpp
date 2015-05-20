// Measurements, second take
// Source: ./examples/measurements2.cpp
#include <qpp.h>
using namespace qpp;
using std::cout;
using std::endl;

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
    cout << ">> Measuring part " << disp(subsys, " ")
            << " of the state: " << endl;
    cout << disp(psi) << endl;
    cout << ">> Measurement result: " << result << endl;
    cout << ">> Probabilities: " << disp(probs, ", ") << endl;
    cout << ">> Resulting normalized post-measurement states: " << endl;

    for (auto&& it: states)
        cout << disp(it) << endl << endl;

    // measure 2 subsystems out of a 4-qubit random density matrix
    cmat rho = randrho(16);
    subsys = {1, 2};
    cmat U = randU(4); // random basis on 2 qubits

    cout << ">> Measuring qubits " << disp(subsys, " ")
            << " of a 4-qubit random state in the random basis:" << endl;
    cout << disp(U) << endl;

    std::tie(result, probs, states) = measure(rho, U, {1, 2});
    cout << ">> Measurement result: " << result << endl;
    cout << ">> Probabilities: " << disp(probs, ", ") << endl;
    cout << ">> Sum of the probabilities: "
            << sum(probs.begin(), probs.end()) << endl;
    cout << ">> Resulting normalized post-measurement states: " << endl;

    for (auto&& it: states)
        cout << disp(it) << endl << endl;

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
    cout << ">> Norm difference: " << norm(rho_bar - rho_out_bar) << endl;

    // random Kraus
    cout << ">> Random channel on part of the state " << endl;
    rho_bar = ptrace(rho, {1, 3});
    rho_out_bar = ptrace(apply(rho, randkraus(3, 4), {1, 3}), {1, 3});

    // verification
    cout << ">> Norm difference: " << norm(rho_bar - rho_out_bar) << endl;

    cout << ">> Sequential measurements on the state/density matrix:" << endl;
    psi = 0.8 * mket({0, 1}) + 0.6 * mket({1, 0});
    rho = psi * adjoint(psi);
    cout << disp(psi) << endl;

    std::vector<idx> subsys_ket{0};
    std::vector<idx> subsys_rho{1};

    auto meas_ket = measure_seq(psi, subsys_ket);
    auto meas_rho = measure_seq(rho, subsys_rho);

    // ket
    cout << ">> Ket, measuring subsystem(s) ";
    cout << disp(subsys_ket, " ") << endl;
    cout << ">> Outcome(s): " << disp(std::get<0>(meas_ket), " ") <<
            endl;
    cout << ">> Probability:  " << std::get<1>(meas_ket) << endl;
    cout << ">> Resulting state:  " << endl;
    cout << disp(std::get<2>(meas_ket)) << endl;

    // density matrix
    cout << ">> Density matrix, measuring subsystem(s) ";
    cout << disp(subsys_rho, " ") << endl;
    cout << ">> Outcome(s): " << disp(std::get<0>(meas_rho), " ") <<
            endl;
    cout << ">> Probability:  " << std::get<1>(meas_rho) << endl;
    cout << ">> Resulting state:  " << endl;
    cout << disp(std::get<2>(meas_rho)) << endl;
}
