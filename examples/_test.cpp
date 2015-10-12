// Used for testing, do not use it as an example
#include <qpp.h>
#include <experimental/experimental.h>
// #include <MATLAB/matlab.h>
using namespace qpp;
using std::cout;
using std::endl;

// TODO: Extensive testing of qpp::ip(), qpp::measure() (rank-one POVMs)
// TODO: modify the pure-states examples that use qpp::measure()
// TODO: add a duration-outputting member function to qpp::Timer

// Future work
// TODO: Implement a circuit description language, load/save from/to file
// TODO: Write a MatrixView class, for multi-index lazy evaluation view

int main()
{
    ket psi = st.GHZ;
    std::vector<idx> subsys{0};
    cmat Ybasis{2, 2};
    Ybasis << 1, 1, 1_i, -1_i;
    Ybasis /= std::sqrt(2);
    cmat V = Ybasis;

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
        ket phi = st.y0;
        ket result = ip(phi, psi, {0});
        cout << "Initial state: \n";
        cout << disp(psi) << endl;
        cout << "Subsystems: " << disp(subsys, " ") << endl;
        cout << "Generalized inner product:\n";
        cout << disp(result) << endl;

        cout << "Additional testing, dim(phi) == dim(psi)\n";
        cout << disp(ip(st.z0, st.y0, {0})) << endl;
        cout << disp(adjoint(st.z0) * st.y0) << endl;
    }

    // Testing ptrace()
    {
        cout << endl << ">> ptrace()\n";
        idx N = 4;
        idx D = std::round(std::pow(2, N));
        std::vector<idx> dims(2, N); // N qubits
        std::vector<idx> subsys(N / 2);
        std::iota(std::begin(subsys), std::end(subsys), 0);
        std::cout << "Generating a random matrix on N = ";
        std::cout << N << " qubits..\n";
        cmat rho = rand<cmat>(D, D);
        cout << "Taking the partial trace over: " << disp(subsys, " ") << endl;
        Timer t;

        // qpp::experimental::ptrace()
        t.tic();
        cmat result2 = experimental::ptrace(rho, subsys);
        t.toc();
        cout << "qpp::experimental::ptrace() took: " << t << " seconds\n";

        // qpp::ptrace
        t.tic();
        cmat result1 = ptrace(rho, subsys);
        t.toc();
        cout << "qpp::ptrace() took: " << t << " seconds\n";

//        cout << disp(result1) << endl << endl;
//        cout << disp(result2) << endl << endl;

        cout << "Norm difference: " << norm(result2 - result1) << endl;
    }
}
