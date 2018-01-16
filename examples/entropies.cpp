// Entropies
// Source: ./examples/entropies.cpp
#include <iostream>

#include "qpp.h"

int main() {
    using namespace qpp;
    cmat rho = st.pb00;
    cmat rhoA = ptrace(rho, {1});
    std::cout << ">> State:\n" << disp(rho) << '\n';
    std::cout << ">> Partial trace over B:\n" << disp(rhoA) << '\n';
    std::cout << ">> von-Neumann entropy: " << entropy(rhoA) << '\n';
    std::cout << ">> Renyi-0 (Hmax) entropy: " << renyi(rhoA, 0) << '\n';
    std::cout << ">> Renyi-1 entropy: " << renyi(rhoA, 1) << '\n';
    std::cout << ">> Renyi-2 entropy: " << renyi(rhoA, 2) << '\n';
    std::cout << ">> Renyi-inf (Hmin) entropy: " << renyi(rhoA, infty) << '\n';
    std::cout << ">> Tsallis-1 entropy: " << tsallis(rhoA, 1) << '\n';
    std::cout << ">> Tsallis-2 entropy: " << tsallis(rhoA, 2) << '\n';
    std::cout << ">> Quantum mutual information between A and B: "
              << qmutualinfo(rho, {0}, {1}) << '\n';
    std::cout << ">> Quantum mutual information between A and A: "
              << qmutualinfo(rho, {0}, {0}) << '\n';
    std::cout << ">> Quantum mutual information between B and B: "
              << qmutualinfo(rho, {1}, {1}) << '\n';
}
