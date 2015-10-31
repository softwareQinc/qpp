// Used for testing, do not use it as an example
#include <qpp.h>
// #include <experimental/experimental.h>
// #include <MATLAB/matlab.h>

using namespace qpp;

int main()
{
    std::vector<idx> dims{2, 2, 2}; // 3 qubits
    idx n = dims.size();            // total number of qudits
    idx D = prod(dims);             // total dimension

    std::vector<idx> ctrl{0, 2};    // where we apply the control
    std::vector<idx> target{1};     // target

    idx Dtarget = 1;                // dimension of the target subsystems
    for(idx i = 0; i < target.size(); ++i)
        Dtarget *= dims[target[i]]; // compute it here

    // some random n qudit pure state
    ket psi = randket(D);

    cmat rho = psi * adjoint(psi); // the corresponding density matrix
    cmat U = randU(Dtarget); // some random unitary on the target

    // applyCTRL on pure state
    ket A = applyCTRL(psi, U, ctrl, target, dims);

    // applyCTRL on density matrix
    cmat B = applyCTRL(rho, U, ctrl, target, dims);

    // result when using CTRL-U|psi><psi|CTRL-U^\dagger
    cmat result_psi = A * adjoint(A);
    // result when using CTRL-U(rho)CTRL-U^\dagger
    cmat result_rho = B;

    std::cout << "|psi><psi| = RHO =" << std::endl;
    std::cout << disp(rho) << std::endl << std::endl;

    std::cout << "U|psi><psi|U^\\dagger" << std::endl;
    std::cout << disp(result_psi) << std::endl << std::endl;

    std::cout << "U rho U^\\dagger" << std::endl;
    std::cout << disp(result_rho) << std::endl << std::endl;

    std::cout << "Norm difference (better be zero): ";
    std::cout << norm(result_psi - result_rho) << std::endl;
}