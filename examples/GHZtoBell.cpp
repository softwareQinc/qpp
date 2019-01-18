//
// Created by Vlad on 2019-01-17.
//

#include <iostream>

#include "qpp.h"
#include "experimental/experimental.h"

int main() {
    using namespace qpp;
    using namespace qpp::experimental;

    idx n = 5;
    QPP_UNUSED_ idx D = std::pow(2, n);
    std::vector<idx> target(n, 0);
    std::iota(target.begin(), target.end(), 0);

    QCircuitDescription qcd(n, n);
    qcd.gate(gt.H, 0);
    for (idx i = 1; i < n; ++i) {
        qcd.gate(gt.CNOT, 0, i);
    }

    qcd.gate(gt.H, 0);
    qcd.gate(gt.H, 1);

    for (idx i = 0; i < n - 2; ++i) {
        qcd.measureV(gt.H, i, i);
    }

    QCircuit qc{qcd};
    qc.run(idx_infty, true);
    idx correction = sum(qc.get_dits()) % 2;
    std::cout << qc;

    ket psi = qc.get_psi();
    ket corrected_psi = psi;
    if (correction == 1)
        corrected_psi = apply(psi, gt.Z, {0});
    std::cout << "psi:\n" << disp(psi) << '\n';
    std::cout << "corrected psi:\n" << disp(corrected_psi) << '\n';
}
