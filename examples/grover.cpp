// Grover's searching
// Source: ./examples/grover.cpp
#include <qpp.h>
using namespace qpp;
using std::cout;
using std::endl;

int main()
{
    idx n = 4; // number of qubits
    cout << ">> Grover on n = " << n << " qubits" << endl;

    std::vector<idx> dims(n, 2); // local dimensions
    // number of elements in the database
    idx N = std::round(std::pow(2, n));
    cout << ">> Database size: " << N << endl;

    // mark an element randomly
    idx marked = randidx(0, N - 1);
    cout << ">> Marked state: " << marked << " -> ";
    cout << disp(n2multiidx(marked, dims), " ") << endl;

    ket psi = mket(n2multiidx(0, dims)); // computational |0>^\otimes n

    // apply H^\otimes n, no aliasing
    psi = (kronpow(gt.H, n) * psi).eval();

    cmat G = 2 * prj(psi) - gt.Id(N); // Diffusion operator

    // number of queries
    idx nqueries = std::ceil(pi * std::sqrt((double) N) / 4.);
    cout << ">> We run " << nqueries << " queries" << endl;
    for (idx i = 0; i < nqueries; ++i)
    {
        psi(marked) = -psi(marked); // apply the oracle first, no aliasing
        psi = (G * psi).eval(); // then the diffusion operator, no aliasing
    }

    // we now measure the state in the computational basis
    auto measured = measure(psi, gt.Id(N));
    cout << ">> Probability of the marked state: "
    << std::get<1>(measured)[marked] << endl;
    cout << ">> Probability of all results: ";
    cout << disp(std::get<1>(measured), ", ") << endl;

    // sample
    cout << ">> Let's sample..." << endl;
    idx result = std::get<0>(measured);
    if (result == marked)
        cout << ">> Hooray, we obtained the correct result: ";
    else
        cout << ">> Not there yet... we obtained: ";
    cout << result << " -> ";
    cout << disp(n2multiidx(result, dims), " ") << endl;
}
