#include <chrono>
#include <iostream>
#include <random>
#include <vector>

#include "qpp/qpp.h"

int main(int arg, char** argv) {
    using namespace qpp;
    using std::chrono::duration;
    using std::chrono::duration_cast;
    using std::chrono::high_resolution_clock;
    using std::chrono::milliseconds;

    if (arg != 3) {
        std::cout
            << "pass an arguments for number of trials and number of qubits\n";
        exit(EXIT_FAILURE);
    }

    idx number_of_trials = std::stoi(argv[1]);
    idx number_of_qubits = std::stoi(argv[2]);
    std::vector<cmat> gate_set = {gt.H, gt.X, gt.Y, gt.Z, gt.S, gt.T};
    ket state = randket();
    for (idx j = 1; j < number_of_qubits; j++) {
        state = kron(state, randket());
    }

    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_int_distribution<idx> gate_dist{0, gate_set.size() - 1};
    std::uniform_int_distribution<idx> qubit_dist{0, number_of_qubits - 1};

    std::vector<std::pair<cmat, idx>> gates_to_apply;
    for (idx j = 0; j < number_of_trials; j++) {
        gates_to_apply.push_back(
            std::make_pair(gate_set[gate_dist(mt)], qubit_dist(mt)));
    }

    duration<double, std::milli> total_time;
    for (auto gate : gates_to_apply) {
        auto t1 = high_resolution_clock::now();
        qpp::apply(state, gate.first, {gate.second});
        auto t2 = high_resolution_clock::now();
        total_time += (t2 - t1);
    }

    std::cout << total_time.count() << std::endl;
}

