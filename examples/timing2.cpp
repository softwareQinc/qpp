// Timing, second take
// Source: ./examples/timing2.cpp
#include <qpp.h>

using namespace qpp;
using std::cout;
using std::endl;

int main()
{
    idx n = 10; // number of qubits
    idx N = std::round(std::pow(2, n));
    cout << ">> n = " << n << " qubits, matrix size " << N << " x " << N
            << "." << endl << endl;
    cmat randcmat = cmat::Random(N, N);

    // qpp::ptrace()
    cout << "**** qpp::ptrace() timing ****" << endl;
    std::vector<idx> subsys_ptrace = {0};
    cout << ">> Subsytem(s): ";
    cout << disp(subsys_ptrace, ", ") << endl;
    Timer<> t;
    ptrace(randcmat, subsys_ptrace);
    cout << ">> It took " << t.toc() << " seconds." << endl << endl;

    // qpp::ptranspose()
    cout << "**** qpp::ptranspose() timing ****" << endl;
    // partially transpose n-1 subsystems
    std::vector<idx> subsys_ptranspose;
    for (idx i = 0; i < n - 1; ++i)
        subsys_ptranspose.push_back(i);
    cout << ">> Subsytem(s): ";
    cout << disp(subsys_ptranspose, ", ") << endl;
    t.tic();
    ptranspose(randcmat, subsys_ptranspose);
    cout << ">> It took " << t.toc() << " seconds." << endl << endl;

    // qpp::syspermute()
    cout << "**** qpp::syspermute() timing ****" << endl;
    std::vector<idx> perm; // left-shift all subsystems by 1
    for (idx i = 0; i < n; ++i)
        perm.push_back((i + 1) % n);
    cout << ">> Subsytem(s): ";
    cout << disp(perm, ", ") << endl;
    t.tic();
    syspermute(randcmat, perm);
    cout << ">> It took " << t.toc() << " seconds." << endl;
}
