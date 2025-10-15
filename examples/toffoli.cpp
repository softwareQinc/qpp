// Source: ./examples/toffoli.cpp
//
// Toffoli gate simulation

#include <iostream>

#include "qpp/qpp.hpp"

int main() {
    using namespace qpp;

    std::cout << ">> Toffoli gate simulation\n\n";

    ket psi_in = randket(8);
    std::cout << ">> Input state:\n";
    std::cout << disp(dirac(psi_in)) << "\n\n";

    /**
     * Toffoli gate (control control not)
     *
     *  ---+---
     *     |
     *  ---+---
     *     |
     *  ---X---
     */
    ket result_qpp = apply(psi_in, gt.TOF, {0, 1, 2});
    std::cout << ">> Toffoli gate output state:\n";
    std::cout << disp(dirac(result_qpp)) << "\n\n";

    /**
     * Toffoli with T and CNOT
     *
     *  -------------+-------------+-----+---T---+--
     *               |             |     |       |
     *  -----+-------------+----------T--X--T_d--X--
     *       |       |     |       |
     *  --H--X--T_d--X--T--X--T_d--X--T--H----------
     */
    ket result = apply(psi_in, gt.H, {2});
    result = applyCTRL(result, gt.X, {1}, {2});
    result = apply(result, adjoint(gt.T), {2});
    result = applyCTRL(result, gt.X, {0}, {2});
    result = apply(result, gt.T, {2});
    result = applyCTRL(result, gt.X, {1}, {2});
    result = apply(result, adjoint(gt.T), {2});
    result = applyCTRL(result, gt.X, {0}, {2});
    result = apply(result, gt.T, {1});
    result = apply(result, gt.T, {2});
    result = applyCTRL(result, gt.X, {0}, {1});
    result = apply(result, gt.T, {0});
    result = apply(result, adjoint(gt.T), {1});
    result = apply(result, gt.H, {2});
    result = applyCTRL(result, gt.X, {0}, {1});
    std::cout << ">> Toffoli with T and CNOT output state:\n";
    std::cout << disp(dirac(result)) << '\n';
    std::cout << ">> Norm difference: " << norm(result - result_qpp) << "\n\n";

    /**
     * https://arxiv.org/abs/quant-ph/9503016 construction
     *
     * V * V = X
     *
     *  -----+-------+---+---
     *       |       |   |
     *  --+--X---+---X-------
     *    |      |       |
     *  --V-----V_d------V---
     */
    cmat sqrtx{cmat::Zero(2, 2)};
    sqrtx << 0.5 + 0.5 * 1_i, 0.5 - 0.5 * 1_i, 0.5 - 0.5 * 1_i, 0.5 + 0.5 * 1_i;
    result = applyCTRL(psi_in, sqrtx, {1}, {2});
    result = applyCTRL(result, gt.X, {0}, {1});
    result = applyCTRL(result, adjoint(sqrtx), {1}, {2});
    result = applyCTRL(result, gt.X, {0}, {1});
    result = applyCTRL(result, sqrtx, {0}, {2});
    std::cout
        << ">> Barenco et. al. [quant-ph/9503016] construction output state:\n";
    std::cout << disp(dirac(result)) << '\n';
    std::cout << ">> Norm difference: " << norm(result - result_qpp) << '\n';
}
