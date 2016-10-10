// Entropies
// Source: ./examples/entropies.cpp
#include <iostream>
#include "qpp.h"

using namespace qpp;

int main()
{
    cmat rho = st.pb00;
    cmat rhoA = ptrace(rho, {1});
    std::cout << ">> State: " << std::endl << disp(rho) << std::endl;
    std::cout << ">> Partial trace over B: " << std::endl << disp(rhoA)
              << std::endl;
    std::cout << ">> von-Neumann entropy: " << entropy(rhoA) << std::endl;
    std::cout << ">> Renyi-0 (Hmax) entropy: " << renyi(rhoA, 0) << std::endl;
    std::cout << ">> Renyi-1 entropy: " << renyi(rhoA, 1) << std::endl;
    std::cout << ">> Renyi-2 entropy: " << renyi(rhoA, 2) << std::endl;
    std::cout << ">> Renyi-inf (Hmin) entropy: " << renyi(rhoA, infty)
              << std::endl;
    std::cout << ">> Tsallis-1 entropy: " << tsallis(rhoA, 1) << std::endl;
    std::cout << ">> Tsallis-2 entropy: " << tsallis(rhoA, 2) << std::endl;
    std::cout << ">> Quantum mutual information between A and B: "
              << qmutualinfo(rho, {0}, {1}) << std::endl;
    std::cout << ">> Quantum mutual information between A and A: "
              << qmutualinfo(rho, {0}, {0}) << std::endl;
    std::cout << ">> Quantum mutual information between B and B: "
              << qmutualinfo(rho, {1}, {1}) << std::endl;
}
