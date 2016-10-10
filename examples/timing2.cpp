// Timing, second take
// Source: ./examples/timing2.cpp
#include <cmath>
#include <iostream>
#include <vector>
#include "qpp.h"

using namespace qpp;

int main()
{
    idx n = 10; // number of qubits
    idx N = std::round(std::pow(2, n));
    std::cout << ">> n = " << n << " qubits, matrix size " << N << " x " << N
              << "." << std::endl << std::endl;
    cmat randcmat = cmat::Random(N, N);

    // qpp::ptrace()
    std::cout << "**** qpp::ptrace() timing ****" << std::endl;
    std::vector<idx> subsys_ptrace = {0};
    std::cout << ">> Subsytem(s): ";
    std::cout << disp(subsys_ptrace, ", ") << std::endl;
    Timer<> t;
    ptrace(randcmat, subsys_ptrace);
    std::cout << ">> It took " << t.toc() << " seconds." << std::endl
              << std::endl;

    // qpp::ptranspose()
    std::cout << "**** qpp::ptranspose() timing ****" << std::endl;
    // partially transpose n-1 subsystems
    std::vector<idx> subsys_ptranspose;
    for (idx i = 0; i < n - 1; ++i)
        subsys_ptranspose.push_back(i);
    std::cout << ">> Subsytem(s): ";
    std::cout << disp(subsys_ptranspose, ", ") << std::endl;
    t.tic();
    ptranspose(randcmat, subsys_ptranspose);
    std::cout << ">> It took " << t.toc() << " seconds." << std::endl
              << std::endl;

    // qpp::syspermute()
    std::cout << "**** qpp::syspermute() timing ****" << std::endl;
    std::vector<idx> perm; // left-shift all subsystems by 1
    for (idx i = 0; i < n; ++i)
        perm.push_back((i + 1) % n);
    std::cout << ">> Subsytem(s): ";
    std::cout << disp(perm, ", ") << std::endl;
    t.tic();
    syspermute(randcmat, perm);
    std::cout << ">> It took " << t.toc() << " seconds." << std::endl;
}
