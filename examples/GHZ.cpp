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
    idx D = std::pow(2, n);
    std::vector<idx> target(n, 0);
    std::iota(target.begin(), target.end(), 0);

    QCircuitDescription qcd(n, 1);
    qcd.gate(gt.H, 0);
    for (idx i = 1; i < n; ++i) {
        qcd.gate(gt.CNOT, 0, i);
    }

    qcd.gate(gt.H, 1);
    qcd.gate(gt.H, 2);
    //qcd.gate(gt.H, 4);

    cmat M = cmat::Zero(D, 1);
    M(0, 0) = 1. / std::sqrt(2);
    M(D - 1, 0) = 1. / std::sqrt(2);

    qcd.measureV(M, target, 0);

    QCircuit qc{qcd};
    qc.run(idx_infty, true);
    std::cout << qc;
    std::cout << "psi:\n" << disp(qc.get_psi()) << '\n';
}
