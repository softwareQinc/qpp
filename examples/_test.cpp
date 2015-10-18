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

    auto viewA = experimental::make_MatrixView(A, {1, 0});
    experimental::MatrixView<cmat> viewB = viewA;

    // this should not compile, and it doesn't :)
    // experimental::MatrixView<cmat> tmpview{gettmp()};

    for (idx i = 0; i < ROWS; ++i)
    {
        for (idx j = 0; j < COLS; ++j)
        {
            std::cout << disp(viewA(i, j)) << " ";
        }
        std::cout << std::endl;
    }

    cmat result = static_cast<cmat>(viewA); // convert
    std::cout << std::endl;
    std::cout << disp(result) << std::endl;
    // std::cout << disp(viewA) << std::endl;
}