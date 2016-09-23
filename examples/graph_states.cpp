// Graph states
// Source: ./examples/graph_states.cpp
#include <qpp.h>

using namespace qpp;
using std::cout;
using std::endl;

int main()
{
    // adjacency matrix, triangle graph (LU equivalent to a GHZ state)
    idx Gamma[3][3] = {{0, 1, 1},
                       {1, 0, 1},
                       {1, 1, 0}};

    // start with 2 states in |000>
    ket G0 = mket({0, 0, 0});
    ket G1 = mket({0, 0, 0});

    // and their density matrices
    cmat rhoG0 = prj(G0);
    cmat rhoG1 = prj(G1);

    // then construct the graph state via 2 methods:
    // qpp::apply() and qpp::applyCTRL()
    // result should be the same, we check later
    cmat H3 = kronpow(gt.H, 3); // all |+>
    G0 = (H3 * G0).eval();
    G1 = G0;
    rhoG0 = (H3 * rhoG0 * adjoint(H3)).eval();
    rhoG1 = rhoG0;
    // apply pairwise Control-Phases
    for ( idx i = 0; i < 3; ++i )
        for ( idx j = i + 1; j < 3; ++j )
        {
            if ( Gamma[i][j] )
            {
                G0 = apply(G0, gt.CZ, {i, j});
                G1 = applyCTRL(G1, gt.Z, {i}, {j});
                rhoG0 = apply(rhoG0, gt.CZ, {i, j});
                rhoG1 = applyCTRL(rhoG1, gt.Z, {i}, {j});
            }
        }
    // end construction

    cout << ">> Resulting graph states: " << endl;
    cout << disp(G0) << endl << endl;
    cout << disp(G1) << endl;
    // verification
    cout << ">> Norm difference: " << norm(G0 - G1) << endl;

    // check the corresponding density matrices
    cout << ">> Resulting density matrices: " << endl;
    cout << disp(rhoG0) << endl << endl;
    cout << disp(rhoG1) << endl;
    cout << ">> Norm difference: " << norm(rhoG0 - rhoG1) << endl;

    // check the X-Z rule
    // applying X to a vertex is equivalent to applying Z to its neighbors
    ket G0X0 = apply(G0, gt.X, {0});
    cmat rhoG0X0 = apply(rhoG0, gt.X, {0});
    ket G0Z1Z2 = apply(G0, kron(gt.Z, gt.Z), {1, 2});
    cmat rhoG0Z1Z2 = apply(rhoG0, kron(gt.Z, gt.Z), {1, 2});

    // verification
    cout << ">> Checking the X-Z rule" << endl;
    cout << ">> X-Z rule. Norm difference for the kets: ";
    cout << norm(G0X0 - G0Z1Z2) << endl;
    cout << ">> X-Z rule. Norm difference for the corresponding "
            "density matrices: ";
    cout << norm(rhoG0X0 - rhoG0Z1Z2) << endl;
}
