// Used for testing, do not use it as an example
#include <qpp.h>
#include <experimental/experimental.h>
// #include <MATLAB/matlab.h>
using namespace qpp;
using std::cout;
using std::endl;

// TODO: circuit description language, load/save from/to file
// TODO: measure pure state
// TODO: check that all functions are either templates or marked inline
// TODO: write a Matrix_view class, for multi-index lazy evaluation view

int main()
{
    ket psi = st.b00;
    auto meas = experimental::_measure(psi, gt.Id2, {0});
    auto m = std::get<0>(meas);
    auto probs = std::get<1>(meas);
    auto outs = std::get<2>(meas);

    std::cout << m << std::endl;
    std::cout << disp(probs, " ") << std::endl << std::endl;
    std::cout << disp(outs[0]) << std::endl << std::endl << disp(outs[1]);
}
