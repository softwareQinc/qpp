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
// TODO: check normalization in measure functions
// TODO: add a duration-outputting member function to Timer
// TODO: throw exception when d = 0

// Future work
// TODO: Implement a circuit description language, load/save from/to file
// TODO: Write a MatrixView class, for multi-index lazy evaluation view

int main()
{
    auto psi = st.GHZ;
    auto subsys{std::vector<idx>{0}};
    cmat Ybasis{2, 2};
    Ybasis << 1, 1, 1_i, -1_i;
    Ybasis /= std::sqrt(2);
    auto V = Ybasis;

    // Testing experimental::_measure()
    {
        std::cout << ">> experimental::_measure() ***\n";
        auto meas = experimental::_measure(psi, V, subsys);
        auto m = std::get<0>(meas);
        auto probs = std::get<1>(meas);
        auto outs = std::get<2>(meas);

        std::cout << "Initial state: \n";
        std::cout << disp(psi) << std::endl;
        std::cout << "Subsystems: " << disp(subsys, " ") << std::endl;
        std::cout << "Measurement matrix: \n";
        std::cout << disp(V) << std::endl;
        std::cout << "Probabilities: " << disp(probs, " ") << std::endl;
        std::cout << "Possible output states: \n";
        for (auto&& elem: outs)
            std::cout << disp(elem) << std::endl << std::endl;
        std::cout << "Result: " << m << std::endl;
        std::cout << "Output state: \n";
        std::cout << disp(outs[0]) << std::endl;
    }

    // Testing experimental::_measure_seq()
    {
        std::cout << std::endl << ">> experimental::_measure_seq()\n";
        auto meas = experimental::_measure_seq(psi, subsys);
        auto results = std::get<0>(meas);
        auto prob = std::get<1>(meas);
        auto out = std::get<2>(meas);

        std::cout << "Initial state: \n";
        std::cout << disp(psi) << std::endl;
        std::cout << "Subsystems: " << disp(subsys, " ") << std::endl;
        std::cout << "Measurement matrix: \n";
        std::cout << disp(V) << std::endl;
        std::cout << "Results: " << disp(results, " ") << std::endl;
        std::cout << "Probability: " << prob << std::endl;
        std::cout << "Output state: \n" << disp(out) << std::endl;
    }

    // Testing experimental::ip()
    {
        std::cout << std::endl << ">> experimental::ip() ***\n";
        psi = 0.8 * mket({0, 0}) + 0.6 * mket({1, 1});
        auto phi = st.z0;
        auto result = experimental::ip(phi, psi, {1});
        std::cout << "Initial state: \n";
        std::cout << disp(psi) << std::endl;
        std::cout << "Subsystems: " << disp(subsys, " ") << std::endl;
        std::cout << "Generalized inner product:\n";
        std::cout << disp(result);
    }
}
