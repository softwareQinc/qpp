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
// TODO: throw exceptions when d = 0

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
        cout << ">> experimental::_measure()\n";
        auto meas = experimental::_measure(psi, V, subsys);
        auto m = std::get<0>(meas);
        auto probs = std::get<1>(meas);
        auto outs = std::get<2>(meas);

        cout << "Initial state: \n";
        cout << disp(psi) << endl;
        cout << "Subsystems: " << disp(subsys, " ") << endl;
        cout << "Measurement matrix: \n";
        cout << disp(V) << endl;
        cout << "Probabilities: " << disp(probs, " ") << endl;
        cout << "Possible output states: \n";
        for (auto&& elem: outs)
            cout << disp(elem) << endl << endl;
        cout << "Result: " << m << endl;
        cout << "Output state: \n";
        cout << disp(outs[0]) << endl;
    }

    // Testing experimental::_measure_seq()
    {
        cout << endl << ">> experimental::_measure_seq()\n";
        auto meas = experimental::_measure_seq(psi, subsys);
        auto results = std::get<0>(meas);
        auto prob = std::get<1>(meas);
        auto out = std::get<2>(meas);

        cout << "Initial state: \n";
        cout << disp(psi) << endl;
        cout << "Subsystems: " << disp(subsys, " ") << endl;
        cout << "Measurement matrix: \n";
        cout << disp(V) << endl;
        cout << "Results: " << disp(results, " ") << endl;
        cout << "Probability: " << prob << endl;
        cout << "Output state: \n" << disp(out) << endl;
    }

    // Testing experimental::ip()
    {
        cout << endl << ">> experimental::ip()\n";
        psi = 0.8 * mket({0, 0}) + 0.6 * mket({1, 1});
        auto phi = st.y1;
        auto result = experimental::ip(phi, psi, {1});
        cout << "Initial state: \n";
        cout << disp(psi) << endl;
        cout << "Subsystems: " << disp(subsys, " ") << endl;
        cout << "Generalized inner product:\n";
        cout << disp(result) << endl;
    }
}
