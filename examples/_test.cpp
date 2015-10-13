// Used for testing, do not use it as an example
#include <qpp.h>
#include <experimental/experimental.h>
// #include <MATLAB/matlab.h>
using namespace qpp;
using std::cout;
using std::endl;

#include <thread>

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
        auto measured = measure(psi, V, subsys);
        idx m = std::get<0>(measured);
        std::vector<double> probs = std::get<1>(measured);
        std::vector<cmat> outs = std::get<2>(measured);

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
        auto measured = measure_seq(psi, subsys);
        std::vector<idx> results = std::get<0>(measured);
        double prob = std::get<1>(measured);
        ket out = std::get<2>(measured);

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
        cout << "ip(): " << disp(ip(st.z0, st.y0, {0})) << endl;
        cout << "Inner product: " << disp(adjoint(st.z0) * st.y0) << endl;
    }

    // Testing Timer<>
    {
        cout << endl << ">> Testing Timer<>\n";
        Timer<std::chrono::duration<float>,
                std::chrono::high_resolution_clock> t;
        std::size_t sec = 2;
        cout << "Waiting " << sec << " seconds...\n";
        std::this_thread::sleep_for(std::chrono::seconds{sec});
        cout << "Done.\n";
        cout << "Time waited in seconds: " << t.toc() << endl;
        cout << "Time waited in milliseconds: ";
        cout << t.get_duration<std::chrono::milliseconds>().count() << endl;
    }
}