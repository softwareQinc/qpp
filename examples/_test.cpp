// Used for testing, do not use it as an example
#include <qpp.h>
#include <experimental/experimental.h>
// #include <MATLAB/matlab.h>
using namespace qpp;

// Future work
// TODO: Implement a circuit description language, load/save from/to file
// TODO: Write a MatrixView class, for multi-index lazy evaluation view - wip

cmat gettmp()
{
    cmat tmp(2, 2);
    tmp << 10, 20, 30, 40;
    return tmp;
}

int main()
{
    // testing qpp::experimental::MatrixView
    idx ROWS = 4, COLS = 4;
    cmat A = gt.CNOT;
    //A << 1, 2., 3., 4.;

    std::cout << "Initial:\n";
    std::cout << disp(A) << std::endl << std::endl;

    auto viewA = experimental::make_MatrixView(A, {1, 0});
    experimental::MatrixView<cmat> viewB = viewA;

    // this should not compile, and it doesn't :)
    // experimental::MatrixView<cmat> tmpview{gettmp()};

    std::cout << "MatrixView: \n";
    std::cout << disp(viewA) << std::endl << std::endl;

    std::cout << "Copy via static_cast of the MatrixView:\n";
    cmat result = static_cast<cmat>(viewA); // convert
    std::cout << disp(result) << std::endl;
}