// Used for testing, do not use it as an example
#include <qpp.h>
#include <experimental/experimental.h>
// #include <MATLAB/matlab.h>
using namespace qpp;
using std::cout;
using std::endl;

// TODO: write a generalized inner product function (non-equal dimensions)
// TODO: measure pure state DONE, need testing
// TODO: modify measure_seq (if necessary)
// TODO: check that all functions are either templates or marked inline
// TODO: modify the pure-states examples accordingly
// TODO: add a duration-outputting member function to Timer

// Future work
// TODO: Implement a circuit description language, load/save from/to file
// TODO: Write a MatrixView class, for multi-index lazy evaluation view

int main()
{
    auto psi = st.GHZ;
    auto meas = experimental::_measure(psi, gt.Id(4), {0, 1});
    auto m = std::get<0>(meas);
    auto probs = std::get<1>(meas);
    auto outs = std::get<2>(meas);

    std::cout << m << std::endl;
    std::cout << disp(probs, " ") << std::endl << std::endl;
    std::cout << disp(outs[0]) << std::endl << std::endl << disp(outs[1]);

    std::cout << std::endl << "HERE\n";

    for(auto&& elem: outs)
        std::cout << disp(elem) << std::endl << std::endl;
}
