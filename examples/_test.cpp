// Used for testing, do not use it as an example
#include <qpp.h>
#include <experimental/experimental.h>
// #include <MATLAB/matlab.h>
using namespace qpp;
using std::cout;
using std::endl;

// TODO: write ptrace in terms of ip(), test it, time it
// TODO: modify the pure-states examples accordingly
// TODO: add a duration-outputting member function to Timer
// TODO: Extensive testing of qpp::ip(), qpp::measure() (rank-one POVMs)

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

    // Testing measure()
    {
        cout << ">> measure()\n";
        auto meas = measure(psi, V, subsys);
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

    // Testing measure_seq()
    {
        cout << endl << ">> measure_seq()\n";
        auto meas = measure_seq(psi, subsys);
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

    // Testing ip()
    {
        cout << endl << ">> ip()\n";
        psi = 0.8 * mket({0, 0}) + 0.6 * mket({1, 1});
        auto phi = st.y0;
        auto result = ip(phi, psi, {0});
        cout << "Initial state: \n";
        cout << disp(psi) << endl;
        cout << "Subsystems: " << disp(subsys, " ") << endl;
        cout << "Generalized inner product:\n";
        cout << disp(result) << endl;

        cout << "Additional testing, dim(phi) == dim(psi)\n";
        cout << disp(ip(st.z0, st.y0, {0})) << endl;
        cout << disp(adjoint(st.z0) * st.y0) << endl;
    }
}
