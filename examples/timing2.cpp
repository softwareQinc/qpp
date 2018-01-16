// Timing, second take
// Source: ./examples/timing2.cpp
#include <cmath>
#include <iostream>
#include <vector>

#include "qpp.h"

int main() {
    using namespace qpp;
    idx n = 10; // number of qubits
    idx N = std::round(std::pow(2, n));
    std::cout << ">> n = " << n << " qubits, matrix size " << N << " x " << N
              << ".\n\n";
    cmat randcmat = cmat::Random(N, N);

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
    std::vector<idx> subsys_ptranspose;
    for (idx i = 0; i < n - 1; ++i)
        subsys_ptranspose.push_back(i);
    std::cout << ">> Subsytem(s): ";
    std::cout << disp(subsys_ptranspose, ", ") << '\n';
    t.tic();
    ptranspose(randcmat, subsys_ptranspose);
    std::cout << ">> It took " << t.toc() << " seconds.\n\n";

    // qpp::syspermute()
    std::cout << "**** qpp::syspermute() timing ****\n";
    std::vector<idx> perm; // left-shift all subsystems by 1
    for (idx i = 0; i < n; ++i)
        perm.push_back((i + 1) % n);
    std::cout << ">> Subsytem(s): ";
    std::cout << disp(perm, ", ") << '\n';
    t.tic();
    syspermute(randcmat, perm);
    std::cout << ">> It took " << t.toc() << " seconds.\n";
}
