// Grover's searching
// Source: ./examples/grover.cpp
#include <cmath>
#include <iostream>
#include <numeric>
#include <tuple>
#include <vector>

#include "qpp.h"

int main() {
    using namespace qpp;
    idx n = 10; // number of qubits
    std::cout << ">> Grover on n = " << n << " qubits\n";

    std::vector<idx> dims(n, 2); // local dimensions
    std::vector<idx> subsys(n);  // ordered subsystems
    std::iota(std::begin(subsys), std::end(subsys), 0);

    // number of elements in the database
    idx N = static_cast<idx>(std::llround(std::pow(2, n)));
    std::cout << ">> Database size: " << N << '\n';

    // mark an element randomly
    idx marked = randidx(0, N - 1);
    std::cout << ">> Marked state: " << marked << " -> ";
    std::cout << disp(n2multiidx(marked, dims), " ") << '\n';

    ket psi = mket(n2multiidx(0, dims)); // computational |0>^\otimes n

    // apply H^\otimes n, no aliasing
    psi = (kronpow(gt.H, n) * psi).eval();

    cmat G = 2 * prj(psi) - gt.Id(N); // Diffusion operator

    Timer<> t;

    // number of queries
    idx nqueries = static_cast<idx>(std::ceil(pi / 4 * std::sqrt(N)));
    std::cout << ">> We run " << nqueries << " queries\n";
    for (idx i = 0; i < nqueries; ++i) {
        psi(marked) = -psi(marked); // apply the oracle first, no aliasing
        psi = (G * psi).eval();     // then the diffusion operator, no aliasing
    }

    // we now measure the state in the computational basis, destructively
    auto measured = measure_seq(psi, subsys, dims);
    std::cout << ">> Probability of the marked state: "
              << std::get<PROB>(measured) << '\n';

    // sample
    std::cout << ">> Let's sample...\n";
    auto result = std::get<RES>(measured);
    if (result == n2multiidx(marked, dims))
        std::cout << ">> Hooray, we obtained the correct result: ";
    else
        std::cout << ">> Not there yet... we obtained: ";
    std::cout << multiidx2n(result, dims) << " -> ";
    std::cout << disp(result, " ") << '\n';

    std::cout << ">> Run time: " << t.toc() << " seconds\n";
}
