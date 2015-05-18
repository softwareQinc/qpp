// Entropies
// Source: ./examples/entropies.cpp
#include <qpp.h>
using namespace qpp;
using std::cout;
using std::endl;

int main()
{
    cmat rho = st.pb00;
    cmat rhoA = ptrace(rho, {1});
    cout << ">> State: " << endl << disp(rho) << endl;
    cout << ">> Partial trace over B: " << endl << disp(rhoA) << endl;
    cout << ">> von-Neumann entropy: " << entropy(rhoA) << endl;
    cout << ">> Renyi-0 (Hmax) entropy: " << renyi(rhoA, 0) << endl;
    cout << ">> Renyi-1 entropy: " << renyi(rhoA, 1) << endl;
    cout << ">> Renyi-2 entropy: " << renyi(rhoA, 2) << endl;
    cout << ">> Renyi-inf (Hmin) entropy: " << renyi(rhoA, infty) << endl;
    cout << ">> Tsallis-1 entropy: " << tsallis(rhoA, 1) << endl;
    cout << ">> Tsallis-2 entropy: " << tsallis(rhoA, 2) << endl;
    cout << ">> Quantum mutual information between A and B: "
            << qmutualinfo(rho, {0}, {1}) << endl;
    cout << ">> Quantum mutual information between A and A: "
            << qmutualinfo(rho, {0}, {0}) << endl;
    cout << ">> Quantum mutual information between B and B: "
            << qmutualinfo(rho, {1}, {1}) << endl;
}
