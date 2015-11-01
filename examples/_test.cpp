// Used for testing, do not use it as an example
#include <qpp.h>
#include <experimental/experimental.h>
// #include <MATLAB/matlab.h>

using namespace qpp;

// General TODOs:
// TODO: write a MatrixView class, for multi-index lazy evaluation view - wip
// TODO: use assert in performance-critical functions (qpp::internal)
// TODO: make sure we use column major order for performance-critical code
// TODO: remove "experimental/experimental.h" include directive at the end
// TODO: test applyCTRL with no controls

// In experimental.h
// TODO: exception checking
// TODO: check for noexcept
// TODO: write the documentation
// TODO: test

// Future work
// TODO: Implement a circuit description language, load/save from/to file

cmat gettmp()
{
    cmat tmp(2, 2);
    tmp << 10, 20, 30, 40;
    return tmp;
}

int main()
{
/*
    // testing qpp::experimental::MatrixView
    idx ROWS = 4, COLS = 4;
    cmat A = gt.CNOT;
    //A << 1, 2., 3., 4.;

    std::cout << "Initial:\n";
    std::cout << disp(A) << std::endl << std::endl;

    auto viewA = experimental::make_MatrixView(A, {1, 0});
    experimental::MatrixView<cmat> viewB = viewA;

    std::cout << "MatrixView: \n";
    std::cout << disp(viewA) << std::endl << std::endl;

    std::cout << "Copy via static_cast/get() of the MatrixView:\n";
    cmat result = static_cast<cmat>(viewA); // convert
    cmat result1 = viewA.get_copy(); // force evaluation
    std::cout << disp(result) << std::endl << std::endl;
    std::cout << disp(result1) << std::endl << std::endl;

    std::cout << "Copy MatrixView:\n";
    auto viewAcopy = viewA; // copy
    std::cout << disp(viewAcopy) << std::endl << std::endl;

    // rvalues should not bind

    // should not compile, and it doesn't :)
    // auto view_expression = experimental::make_MatrixView(A + A, {0, 1});
    // std::cout << disp(view_expression) << std::endl;

    // this line should not compile, and it doesn't :)
    // experimental::MatrixView<cmat> tmpview{gettmp(), {0, 1}};

    std::cout << "MatrixValue as a rvalue:\n";
    std::cout << disp(qpp::experimental::make_MatrixView(A, {1, 0}));
    std::cout << std::endl;

    std::cout << "Testing MatrixView via syspermute:\n";
    std::vector<idx> perm{1, 2, 3, 7, 4, 6, 5, 0};
    idx N = std::pow(2, perm.size());
    cmat rho = qpp::rand<cmat>(256, 256);
    cmat B = syspermute(rho, perm);
    cmat C = experimental::make_MatrixView(rho, perm).get_copy();
    std::cout << "Norm difference: " << norm(B - C) << std::endl;
*/
    std::vector<idx> dims{2, 2, 2, 2};  // 3 qubits
    idx n = dims.size();                // total number of qudits
    idx D = prod(dims);                 // total dimension

    std::vector<idx> ctrl{2};           // where we apply the control
    std::vector<idx> target{1, 0, 3};   // target

    idx Dtarget = 1;                    // dimension of the target subsystems
    for(idx i = 0; i < target.size(); ++i)
        Dtarget *= dims[target[i]];     // compute it here

    // some random n qudit pure state
    ket psi = randket(D);

    cmat rho = psi * adjoint(psi); // the corresponding density matrix
    cmat U = randU(Dtarget);       // some random unitary on the target

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
