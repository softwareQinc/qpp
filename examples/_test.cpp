// Used for testing, do not use it as an example
#include <iostream>

#include "qpp.h"
#include "experimental/experimental.h"

int main() {
    using namespace qpp;
    using namespace qpp::experimental;

    QCirc<int> qc(10, 10);
    qc.apply_all(gt.H);
    std::cout << qc.get_num_active_qubits() << '\n';
    qc.measure({3, 1, 7});
    std::cout << qc.get_num_active_qubits() << '\n';
    std::cout << qc.get_num_measured_qubits() << '\n';
    std::cout << qc.get_num_active_qubits() << '\n';
    qc.measure({2, 4, 5, 6, 0});
    std::cout << disp(qc.get_results(), " ") << '\n';
    qc.apply_all(gt.H);
    qc.apply_all(gt.X);
    qc.measure_all();
    // qc.measure_all();
    std::cout << disp(qc.get_results(), "") << '\n';
    std::cout << qc.get_results_as_N() << '\n';
    std::cout << qc.get_nq() << '\n';
    std::cout << qc.get_num_measured_qubits() << '\n';
    std::cout << qc.get_num_active_qubits() << "\n\n";
    std::cout << "\n\n";

    QCirc<int> qc1(2, 0);
    qc1.apply(gt.H, {0});
    qc1.applyCTRL(gt.X, {0}, {1});
    qc1.apply(gt.CNOT, {0, 1});
    qc1.measure({0});

    std::cout << qc1.get_num_measured_qubits() << '\n';
    std::cout << qc1.get_num_active_qubits() << '\n';
    std::cout << disp(qc1.get_results(), " ") << '\n';
    std::cout << disp(qc1.get_psi()) << "\n\n";

    qc1.reset();
    std::cout << qc1.get_num_measured_qubits() << '\n';
    std::cout << qc1.get_num_active_qubits() << '\n';
    std::cout << disp(qc1.get_results(), " ") << '\n';
    std::cout << disp(qc1.get_psi()) << "\n\n";

    QCirc<int> qc2(2, 10, 3);
    qc2.apply(gt.Xd(3), {0});
    qc2.measure({1});
    std::cout << disp(qc2.get_psi()) << "\n\n";
    std::cout << qc2.get_size() << " " << qc2.d() << "\n\n";

    using mytype = std::tuple<cmat, std::vector<idx>, std::vector<idx>>;
    std::vector<mytype> gates;

    cmat X(2, 2);
    X << 0, 0, 0, 1;
    auto elem = std::make_tuple(std::cref(X), std::vector<idx>({0}),
                                std::vector<idx>({1, 2, 3}));
    std::cout << disp(std::get<0>(elem)) << '\n';
    X(1, 1) = 2;
    std::cout << disp(std::get<0>(elem)) << '\n';
    std::cout << disp(X) << '\n';

    LogicalCircuit lc;
    lc.add(gt.X, "X", {0}, {1, 2});
    lc.add(gt.Z, "Z", {}, {1, 2});
    lc.add(gt.CNOT, "CNOT", {0}, {2});
    lc.add(gt.TOF, "TOF", {0, 1}, {2});
    std::cout << lc << '\n';
    std::cout << lc.get_gate_count() << "\n\n";
}
