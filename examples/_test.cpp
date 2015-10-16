// Used for testing, do not use it as an example
#include <qpp.h>
#include <experimental/experimental.h>
// #include <MATLAB/matlab.h>
using namespace qpp;

// Future work
// TODO: Implement a circuit description language, load/save from/to file
// TODO: Write a MatrixView class, for multi-index lazy evaluation view - wip

int main()
{
    // testing qpp::experimental::MatrixView
    idx ROWS = 2, COLS = 2;
    cmat A(ROWS, COLS);
    A << 1, 2., 3., 4.;

    auto viewA = experimental::make_MatrixView(A);
    experimental::MatrixView<cmat> viewB = viewA;

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
}