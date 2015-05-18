// Quantum channels
// Source: ./examples/channels.cpp
#include <qpp.h>
using namespace qpp;
using std::cout;
using std::endl;

int main()
{
    idx nk = 5;
    idx D = 3; // nk Kraus on d-dimensional system
    cout << ">> Generating a random channel with " << nk
    << " Kraus operators on a " << D << " dimensional space" << endl;
    std::vector<cmat> Ks = randkraus(nk, D);

    cmat rho_in = randrho(D); // random input state
    cmat rho_out = apply(rho_in, Ks); // output state

    cout << ">> Computing its Choi matrix..." << endl;
    cmat choim = kraus2choi(Ks);
    cout << ">> Choi matrix:" << endl << disp(choim) << endl;

    cout << ">> The eigenvalues of the Choi matrix are: "
    << endl << disp(transpose(hevals(choim))) << endl;

    cout << ">> Their sum is: " << sum(hevals(choim)) << endl;

    std::vector<cmat> Kperps = choi2kraus(choim);
    cout << ">> The Kraus rank of the channel is: " << Kperps.size() << endl;

    cmat rho_out1 = apply(rho_in, Kperps);
    // verification
    cout << ">> Norm difference on output states: "
    << norm(rho_out1 - rho_out) << endl;

    cout << ">> Superoperator matrix:" << endl;
    cmat smat = kraus2super(Ks);
    cout << disp(smat) << endl;

    cout << ">> The eigenvalues of the superoperator matrix are: " << endl;
    dyn_col_vect<cplx> evalsupop = evals(smat);
    cout << disp(transpose(evalsupop)) << endl;

    cout << ">> Their absolute values are: " << endl;
    for (idx i = 0; i < (idx) evalsupop.size(); ++i)
        cout << std::abs(evalsupop(i)) << " ";

    // verification
    cout << endl << ">> Norm difference for the superoperator action: ";
    cmat rho_out2 = transpose(
            reshape(smat * reshape(transpose(rho_in), D * D, 1), D, D));
    cout << norm(rho_out - rho_out2) << endl;
}
