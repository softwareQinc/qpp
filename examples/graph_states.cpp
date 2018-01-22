// Graph states
// Source: ./examples/graph_states.cpp
#include <iostream>

#include "qpp.h"

int main() {
    using namespace qpp;
    // adjacency matrix, triangle graph (LU equivalent to a GHZ state)
    idx Gamma[3][3] = {{0, 1, 1}, {1, 0, 1}, {1, 1, 0}};

    // start with 2 states in |000>
    ket G0 = 000_ket;
    ket G1 = 000_ket;

    // and their density matrices
    cmat rhoG0 = prj(G0);
    cmat rhoG1 = prj(G1);

    // then construct the graph state via 2 methods:
    // qpp::apply() and qpp::applyCTRL()
    // result should be the same, we check later
    cmat H3 = kronpow(gt.H, 3); // all |+>
    G0 = (H3 * G0).eval();
    G1 = G0;
    rhoG0 = (H3 * rhoG0 * adjoint(H3)).eval();
    rhoG1 = rhoG0;
    // apply pairwise Control-Phases
    for (idx i = 0; i < 3; ++i)
        for (idx j = i + 1; j < 3; ++j) {
            if (Gamma[i][j]) {
                G0 = apply(G0, gt.CZ, {i, j});
                G1 = applyCTRL(G1, gt.Z, {i}, {j});
                rhoG0 = apply(rhoG0, gt.CZ, {i, j});
                rhoG1 = applyCTRL(rhoG1, gt.Z, {i}, {j});
            }
        }
    // end construction

    std::cout << ">> Resulting graph states:\n";
    std::cout << disp(G0) << "\n\n";
    std::cout << disp(G1) << '\n';
    // verification
    std::cout << ">> Norm difference: " << norm(G0 - G1) << '\n';

    // check the corresponding density matrices
    std::cout << ">> Resulting density matrices:\n";
    std::cout << disp(rhoG0) << "\n\n";
    std::cout << disp(rhoG1) << '\n';
    std::cout << ">> Norm difference: " << norm(rhoG0 - rhoG1) << '\n';

    // check the X-Z rule
    // applying X to a vertex is equivalent to applying Z to its neighbors
    ket G0X0 = apply(G0, gt.X, {0});
    cmat rhoG0X0 = apply(rhoG0, gt.X, {0});
    ket G0Z1Z2 = apply(G0, kron(gt.Z, gt.Z), {1, 2});
    cmat rhoG0Z1Z2 = apply(rhoG0, kron(gt.Z, gt.Z), {1, 2});

    // verification
    std::cout << ">> Checking the X-Z rule\n";
    std::cout << ">> X-Z rule. Norm difference for the kets: ";
    std::cout << norm(G0X0 - G0Z1Z2) << '\n';
    std::cout << ">> X-Z rule. Norm difference for the corresponding "
              << "density matrices: ";
    std::cout << norm(rhoG0X0 - rhoG0Z1Z2) << '\n';
}
