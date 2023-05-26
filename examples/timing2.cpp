// Timing
// Source: ./examples/timing2.cpp
#include <cmath>
#include <iostream>
#include <numeric>
#include <vector>

#include "qpp.h"

int main() {
    using namespace qpp;

    idx n = 10;                                              // number of qubits
    auto D = static_cast<idx>(std::llround(std::pow(2, n))); // dimension 2^n
    std::cout << ">> n = " << n << " qubits, matrix size " << D << " x " << D
              << ".\n\n";
    cmat randcmat = cmat::Random(D, D);

    // qpp::ptrace()
    std::cout << "**** qpp::ptrace() timing ****\n";
    std::vector<idx> subsys_ptrace = {0};
    std::cout << ">> Subsytem(s): ";
    std::cout << disp(subsys_ptrace, ", ") << '\n';
    Timer<> t;
    ptrace(randcmat, subsys_ptrace);
    std::cout << ">> It took " << t.toc() << " seconds.\n\n";

    // qpp::ptranspose()
    std::cout << "**** qpp::ptranspose() timing ****\n";
    // partially transpose n-1 subsystems
    std::vector<idx> subsys_ptranspose(n - 1);
    std::iota(std::begin(subsys_ptranspose), std::end(subsys_ptranspose), 0);
    std::cout << ">> Subsytem(s): ";
    std::cout << disp(subsys_ptranspose, ", ") << '\n';
    t.tic();
    ptranspose(randcmat, subsys_ptranspose);
    std::cout << ">> It took " << t.toc() << " seconds.\n\n";

    // qpp::syspermute()
    std::cout << "**** qpp::syspermute() timing ****\n";
    std::vector<idx> perm(n); // left-shift all subsystems by 1
    for (idx i = 0; i < n; ++i)
        perm[i] = (i + 1) % n;
    std::cout << ">> Subsytem(s): ";
    std::cout << disp(perm, ", ") << '\n';
    t.tic();
    syspermute(randcmat, perm);
    std::cout << ">> It took " << t.toc() << " seconds.\n";
}
