/*
 * Quantum++
 *
 * Copyright (c) 2013 - 2014 Vlad Gheorghiu (vgheorgh@gmail.com)
 *
 * This file is part of Quantum++.
 *
 * Quantum++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Quantum++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Quantum++.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <qpp.h>

// #include <MATLAB/matlab.h> // support for MATLAB

using namespace qpp;

using std::cout;
using std::endl;

// We define each example as an independent function
void MEASUREMENTS()
{
    cout << "**** Measurements ****" << endl;

    ket psi = st.b00;
    std::vector<idx> subsys = {0};
    idx result;
    std::vector<double> probs;
    std::vector<cmat> states;

    // measures the first subsystem of the Bell state (|00> + |11>)/sqrt(2)
    // in the X basis
    std::tie(result, probs, states) = measure(psi, gt.H, {0});
    cout << ">> Measuring part " << disp(subsys, " ")
            << " of the state: " << endl;
    cout << disp(psi) << endl;
    cout << ">> Measurement result: " << result << endl;
    cout << ">> Probabilities: " << disp(probs, ", ") << endl;
    cout << ">> Resulting normalized post-measurement states: " << endl;

    for (auto&& it: states)
        cout << disp(it) << endl << endl;

    // measure 2 subsystems out of a 4-qubit random density matrix
    cmat rho = randrho(16);
    subsys = {1, 2};
    cmat U = randU(4); // random basis on 2 qubits

    cout << ">> Measuring qubits " << disp(subsys, " ")
            << " of a 4-qubit random state in the random basis:" << endl;
    cout << disp(U) << endl;

    std::tie(result, probs, states) = measure(rho, U, {1, 2});
    cout << ">> Measurement result: " << result << endl;
    cout << ">> Probabilities: " << disp(probs, ", ") << endl;
    cout << ">> Sum of the probabilities: "
            << sum(probs.begin(), probs.end()) << endl;
    cout << ">> Resulting normalized post-measurement states: " << endl;

    for (auto&& it: states)
        cout << disp(it) << endl << endl;

    // check now how the state after the measurement "looks"
    // on the left over subsystems {0, 3}

    // it should be the same as the partial trace over {1, 2}
    // of the initial state (before the measurement), as local CPTP maps
    // do not influence the complementary subsystems

    cmat rho_bar = ptrace(rho, subsys);
    cmat rho_out_bar = cmat::Zero(4, 4);

    // compute the resulting mixed state after the measurement
    for (idx i = 0; i < probs.size(); ++i)
        rho_out_bar += probs[i] * states[i];

    // verification
    cout << ">> Norm difference: " << norm(rho_bar - rho_out_bar) << endl;

    // random Kraus
    cout << ">> Random channel on part of the state " << endl;
    rho_bar = ptrace(rho, {1, 3});
    rho_out_bar = ptrace(apply(rho, randkraus(3, 4), {1, 3}), {1, 3});

    // verification
    cout << ">> Norm difference: " << norm(rho_bar - rho_out_bar)
            << endl << endl;
}

void TELEPORTATION()
{
    idx D = 3; // size of the system
    cout << "**** Qudit teleportation, D = " << D << " ****" << endl;

    ket mes_AB = ket::Zero(D * D); // maximally entangled state resource
    for (idx i = 0; i < D; ++i)
        mes_AB += mket({i, i}, D);
    mes_AB /= std::sqrt((double) D);

    // circuit that measures in the qudit Bell basis
    cmat Bell_aA = adjoint(gt.CTRL(gt.Xd(D), {0}, {1}, 2, D)
            * kron(gt.Fd(D), gt.Id(D)));

    ket psi_a = randket(D); // random state as input on a
    cout << ">> Initial state:" << endl;
    cout << disp(psi_a) << endl;

    ket input_aAB = kron(psi_a, mes_AB); // joint input state aAB
    // output before measurement
    ket output_aAB = apply(input_aAB, Bell_aA, {0, 1}, D);

    // measure on aA
    auto measured_aA = measure(output_aAB, gt.Id(D * D), {0, 1}, D);
    idx m = std::get<0>(measured_aA); // measurement result

    auto midx = n2multiidx(m, {D, D});
    cout << ">> Alice's measurement result: ";
    cout << m << " -> " << disp(midx, " ") << endl;
    cout << ">> Alice's measurement probabilities: ";
    cout << disp(std::get<1>(measured_aA), ", ") << endl;

    // conditional result on B before correction
    cmat output_m_B = std::get<2>(measured_aA)[m];
    cmat correction_B = powm(gt.Zd(D), midx[0]) *
            powm(adjoint(gt.Xd(D)), midx[1]); // correction operator
    // apply correction on B
    cout << ">> Bob must apply the correction operator Z^" << midx[0]
            << " X^" << D - midx[1] << endl;
    cmat rho_B = correction_B * output_m_B * adjoint(correction_B);

    cout << ">> Bob's final state (after correction): " << endl;
    cout << disp(rho_B) << endl;

    // verification
    cout << ">> Norm difference: " << norm(rho_B - prj(psi_a)) << endl << endl;
}

void DENSE_CODING()
{
    idx D = 3; // size of the system
    cout << "**** Qudit dense coding, D = " << D << " ****" << endl;

    ket mes_AB = ket::Zero(D * D); // maximally entangled state resource
    for (idx i = 0; i < D; ++i)
        mes_AB += mket({i, i}, D);
    mes_AB /= std::sqrt((double) D);

    // circuit that measures in the qudit Bell basis
    cmat Bell_AB = adjoint(gt.CTRL(gt.Xd(D), {0}, {1}, 2, D)
            * kron(gt.Fd(D), gt.Id(D)));

    // equal probabilities of choosing a message
    idx m_A = randidx(0, D * D - 1);
    auto midx = n2multiidx(m_A, {D, D});
    cout << ">> Alice sent: " << m_A << " -> ";
    cout << disp(midx, " ") << endl;

    // Alice's operation
    cmat U_A = powm(gt.Zd(D), midx[0]) * powm(adjoint(gt.Xd(D)), midx[1]);

    // Alice encodes the message
    ket psi_AB = apply(mes_AB, U_A, {0}, D);

    // Bob measures the joint system in the qudit Bell basis
    // use swap to avoid temporaries until Eigen supports move semantics,
    // although the compiler should use copy elision and
    // remove the temporary even without using the swap trick
    psi_AB.swap(apply(psi_AB, Bell_AB, {0, 1}, D));

    auto measured = measure(psi_AB, gt.Id(D * D));
    cout << ">> Bob's measurement probabilities: ";
    cout << disp(std::get<1>(measured), ", ") << endl;

    // Bob samples according to the measurement probabilities
    idx m_B = std::get<0>(measured);
    cout << ">> Bob received: ";
    cout << m_B << " -> " << disp(n2multiidx(m_B, {D, D}), " ")
            << endl << endl;
}

void GROVER()
{
    Timer t; // set a timer

    idx n = 4; // number of qubits
    cout << "**** Grover on n = " << n << " qubits ****" << endl;

    std::vector<idx> dims(n, 2); // local dimensions
    // number of elements in the database
    idx N = std::round(std::pow(2, n));
    cout << ">> Database size: " << N << endl;

    // mark an element randomly
    idx marked = randidx(0, N - 1);
    cout << ">> Marked state: " << marked << " -> ";
    cout << disp(n2multiidx(marked, dims), " ") << endl;

    ket psi = mket(n2multiidx(0, dims)); // computational |0>^\otimes n

    // apply H^\otimes n, no aliasing
    psi = (kronpow(gt.H, n) * psi).eval();

    cmat G = 2 * prj(psi) - gt.Id(N); // Diffusion operator

    // number of queries
    idx nqueries = std::ceil(pi * std::sqrt((double) N) / 4.);
    cout << ">> We run " << nqueries << " queries" << endl;
    for (idx i = 0; i < nqueries; ++i)
    {
        psi(marked) = -psi(marked); // apply the oracle first, no aliasing
        psi = (G * psi).eval(); // then the diffusion operator, no aliasing
    }

    // we now measure the state in the computational basis
    auto measured = measure(psi, gt.Id(N));
    cout << ">> Probability of the marked state: "
            << std::get<1>(measured)[marked] << endl;
    cout << ">> Probability of all results: ";
    cout << disp(std::get<1>(measured), ", ") << endl;

    // sample
    cout << ">> Let's sample..." << endl;
    idx result = std::get<0>(measured);
    if (result == marked)
        cout << ">> Hooray, we obtained the correct result: ";
    else
        cout << ">> Not there yet... we obtained: ";
    cout << result << " -> ";
    cout << disp(n2multiidx(result, dims), " ") << endl;

    // stop the timer and display it
    cout << ">> It took " << t.toc() << " seconds to simulate Grover on "
            << n << " qubits." << endl << endl;
}

void ENTANGLEMENT()
{
    cout << "**** Entanglement ****" << endl;

    cmat rho = 0.2 * st.pb00 + 0.8 * st.pb11;
    cout << ">> rho: " << endl;
    cout << disp(rho) << endl;

    cout << ">> Concurrence of rho: " << concurrence(rho) << endl;

    cout << ">> Negativity of rho: " << negativity(rho, {2, 2}) << endl;

    cout << ">> Logarithimc negativity of rho: "
            << lognegativity(rho, {2, 2}) << endl;

    ket psi = 0.8 * mket({0, 0}) + 0.6 * mket({1, 1});

    // apply some local random unitaries
    psi = kron(randU(2), randU(2)) * psi;

    cout << ">> psi: " << endl;
    cout << disp(psi) << endl;

    cout << ">> Entanglement of psi: " << entanglement(psi, {2, 2}) << endl;

    cout << ">> Concurrence of psi: " << concurrence(prj(psi)) << endl;

    cout << ">> G-Concurrence of psi: " << gconcurrence(psi) << endl;

    cout << ">> Schmidt coefficients of psi: " << endl;
    cout << disp(schmidtcoeff(psi, {2, 2})) << endl;

    cout << ">> Schmidt probabilities of psi: " << endl;
    cout << disp(schmidtprob(psi, {2, 2})) << endl;

    cmat UA = schmidtA(psi, {2, 2});
    cmat UB = schmidtB(psi, {2, 2});

    cout << ">> Schmidt vectors on Alice's side: " << endl;
    cout << disp(UA) << endl;

    cout << ">> Schmidt vectors on Bob's side: " << endl;
    cout << disp(UB) << endl;

    cout << ">> State psi in the Schmidt basis: " << endl;
    cout << disp(adjoint(kron(UA, UB)) * psi) << endl;

    // reconstructed state
    ket psi_from_schmidt =
            schmidtcoeff(psi, {2, 2})(0) * kron(UA.col(0), UB.col(0))
                    + schmidtcoeff(psi, {2, 2})(1) * kron(UA.col(1), UB.col(1));
    cout << ">> State psi reconstructed from the Schmidt decomposition: "
            << endl;
    cout << disp(psi_from_schmidt) << endl;

    // verification
    cout << ">> Norm difference: " << norm(psi - psi_from_schmidt)
            << endl << endl;
}

void QECC()
{
    cout << "**** Quantum error correcting codes ****" << endl;

    ket a0 = codes.codeword(Codes::Type::FIVE_QUBIT, 0);
    ket a1 = codes.codeword(Codes::Type::FIVE_QUBIT, 1);

    ket b0 = codes.codeword(Codes::Type::SEVEN_QUBIT_STEANE, 0);
    ket b1 = codes.codeword(Codes::Type::SEVEN_QUBIT_STEANE, 1);

    ket c0 = codes.codeword(Codes::Type::NINE_QUBIT_SHOR, 0);
    ket c1 = codes.codeword(Codes::Type::NINE_QUBIT_SHOR, 1);

    cout << ">> [[5, 1, 3]] Five qubit code. ";
    cout << "Checking codeword orthogonality." << endl;
    cout << ">> |<0L|1L>| = ";
    cout << disp(adjoint(a0) * a1) << endl;

    cout << ">> [[7, 1, 3]] Seven qubit Steane code. ";
    cout << "Checking codeword orthogonality." << endl;
    cout << ">> |<0L|1L>| = ";
    cout << disp(adjoint(b0) * b1) << endl;

    cout << ">> [[9, 1, 3]] Nine qubit Shor code. ";
    cout << "Checking codeword orthogonality." << endl;
    cout << ">> |<0L|1L>| = ";
    cout << disp(adjoint(c0) * c1) << endl << endl;
}

void CHANNEL()
{
    cout << "**** Channel tests ****" << endl;

    idx nk = 5;
    idx D = 3; // nk Kraus on d-dimensional system
    cout << ">> Generating a random channel with " << nk
            << " Kraus operators on a " << D << " dimensional space" << endl;
    std::vector<cmat> Ks = randkraus(nk, D);

    cmat rho_in = randrho(D); // random input state
    cmat rho_out = apply(rho_in, Ks); // output state

    cout << ">> Computing its Choi matrix..." << endl;
    cmat choim = choi(Ks);
    cout << ">> Choi matrix:" << endl << disp(choim) << endl;

    cout << ">> The eigenvalues of the Choi matrix are: "
            << endl << disp(transpose(hevals(choim))) << endl;

    cout << ">> Their sum is: " << sum(hevals(choim)) << endl;

    std::vector<cmat> Kperps = choi2kraus(choim);
    cout << ">> The Kraus rank of the channel is: " << Kperps.size() << endl;

    cmat rho_out1 = apply(rho_in, Kperps);
    // verification
    cout << ">> Norm difference on output states: "
            << norm(rho_out1 - rho_out) << endl;

    cout << ">> Superoperator matrix:" << endl;
    cmat smat = super(Ks);
    cout << disp(smat) << endl;

    cout << ">> The eigenvalues of the superoperator matrix are: " << endl;
    cmat evalsupop = evals(smat);
    cout << disp(transpose(evalsupop)) << endl;

    cout << ">> Their absolute values are: " << endl;
    for (idx i = 0; i < (idx) evalsupop.size(); ++i)
        cout << std::abs(evalsupop(i)) << " ";

    // verification
    cout << endl << ">> Norm difference for the superoperator action: ";
    cmat rho_out2 = transpose(
            reshape(smat * reshape(transpose(rho_in), D * D, 1), D, D));
    cout << norm(rho_out - rho_out2) << endl << endl;
}

// test function used by qpp::cwise() in FUNCTOR()
cplx pow3(const cplx& z)
{
    return std::pow(z, 3);
}

void FUNCTOR()
{
    cout << "**** Functor ****" << endl;

    // functor test
    cout << ">> Functor z^3 acting component-wise on:" << endl;
    cmat A(2, 2);
    A << 1, 2, 3, 4;
    cout << disp(A) << endl;

    cout << ">> Result (with lambda):" << endl;
    // functor z^3 componentwise, specify OutputScalar and Derived for lambdas
    cout << disp(cwise<cplx, cmat>(A, [](const cplx& z) -> cplx
    {
        return z * z * z;
    })) << endl;

    cout << ">> Result (with proper function):" << endl;
    // automatic type deduction for proper functions
    cout << disp(cwise(A, &pow3)) << endl << endl;
}

void GRAMSCHMIDT()
{
    cout << "**** Gram-Schmidt ****" << endl;

    cmat A(3, 3);
    A << 1, 1, 0, 0, 2, 0, 0, 0, 0;
    cout << ">> Input matrix:" << endl << disp(A) << endl;

    cmat Ags = grams(A);
    cout << ">> Result:" << endl << disp(Ags) << endl;

    cout << ">> Projector onto G.S. vectors:" << endl;
    cout << disp(Ags * adjoint(Ags)) << endl << endl;
}

void SPECTRAL()
{
    cout << "**** Spectral decomposition tests ****" << endl;
    idx D = 4;
    cmat rH = randH(D);
    cout << ">> Original matrix: " << endl << disp(rH) << endl;

    // spectral decomposition here
    dyn_col_vect<double> evalsH = hevals(rH);
    cmat evectsH = hevects(rH);
    cmat spec = cmat::Zero(D, D);
    // reconstruct the matrix
    for (idx i = 0; i < D; ++i)
        spec += evalsH(i) * prj(evectsH.col(i));

    cout << ">> Reconstructed from spectral decomposition: " << endl;
    cout << disp(spec) << endl;

    // verification
    cout << ">> Norm difference: " << norm(spec - rH) << endl << endl;
}

void RANDOM()
{
    cout << "**** Randomness ****" << endl;

    cout << ">> Generating a random ket on D = 5" << endl;
    ket rket = randket(5);
    cout << disp(rket) << endl;

    std::vector<double> probs = abssq(rket);
    cout << "Probabilities: " << disp(probs, ", ") << endl;

    cout << "Sum of the probabilities: ";
    cout << sum(probs.begin(), probs.end()) << endl << endl;
}

void ENTROPIES()
{
    cout << "*** Entropies ****" << endl;

    cmat rho = st.pb00;
    cmat rhoA = ptrace(rho, {1});
    cout << ">> State: " << endl << disp(rho) << endl;
    cout << ">> Partial trace over B: " << endl << disp(rhoA) << endl;
    cout << ">> Shannon entropy: " << shannon(rhoA) << endl;
    cout << ">> Renyi-0 (Hmax) entropy: " << renyi(rhoA, 0) << endl;
    cout << ">> Renyi-1 entropy: " << renyi(rhoA, 1) << endl;
    cout << ">> Renyi-2 entropy: " << renyi(rhoA, 2) << endl;
    cout << ">> Renyi-inf (Hmin) entropy: " << renyi(rhoA, infty) << endl;
    cout << ">> Tsallis-1 entropy: " << tsallis(rhoA, 1) << endl;
    cout << ">> Tsallis-2 entropy: " << tsallis(rhoA, 2) << endl;
    cout << ">> Quantum mutual information between A and B: "
            << qmutualinfo(rho, {0}, {1}) << endl;
    cout << ">> Quantum mutual information between A and A: "
            << qmutualinfo(rho, {0}, {0}) << endl;
    cout << ">> Quantum mutual information between B and B: "
            << qmutualinfo(rho, {1}, {1}) << endl << endl;
}

void GRAPHSTATES()
{
    cout << "**** Graph states ****" << endl;

    // adjacency matrix, triangle graph (LU equivalent to a GHZ state)
    idx Gamma[3][3] = {{0, 1, 1}, {1, 0, 1}, {1, 1, 0}};

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
    for (idx i = 0; i < 3; ++i)
        for (idx j = i + 1; j < 3; ++j)
        {
            if (Gamma[i][j])
            {
                // use swap to avoid temporaries
                // until Eigen supports move semantics,
                // although the compiler should use copy elision and
                // remove the temporary even without using the swap trick
                G0.swap(apply(G0, gt.CZ, {i, j}));
                G1.swap(applyCTRL(G1, gt.Z, {i}, {j}));
                rhoG0.swap(apply(rhoG0, gt.CZ, {i, j}));
                rhoG1.swap(applyCTRL(rhoG1, gt.Z, {i}, {j}));
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
    cout << norm(rhoG0X0 - rhoG0Z1Z2) << endl << endl;
}

void TIMING()
{
    cout << "**** Timing tests ****" << endl;

    idx n = 12; // number of qubits
    idx N = std::round(std::pow(2, n));
    cout << ">> n = " << n << " qubits, matrix size " << N << " x " << N
            << "." << endl << endl;
    cmat randcmat = cmat::Random(N, N);

    // qpp::ptrace()
    cout << "**** qpp::ptrace() timing ****" << endl;
    std::vector<idx> subsys_ptrace = {0};
    cout << ">> Subsytem(s): ";
    cout << disp(subsys_ptrace, ", ") << endl;
    Timer t;
    ptrace(randcmat, subsys_ptrace);
    cout << ">> It took " << t.toc() << " seconds." << endl << endl;

    // qpp::ptranspose()
    cout << "**** qpp::ptranspose() timing ****" << endl;
    // partially transpose n-1 subsystems
    std::vector<idx> subsys_ptranspose;
    for (idx i = 0; i < n - 1; ++i)
        subsys_ptranspose.push_back(i);
    cout << ">> Subsytem(s): ";
    cout << disp(subsys_ptranspose, ", ") << endl;
    t.tic();
    ptranspose(randcmat, subsys_ptranspose);
    cout << ">> It took " << t.toc() << " seconds." << endl << endl;

    // qpp::syspermute()
    cout << "**** qpp::syspermute() timing ****" << endl;
    std::vector<idx> perm; // left-shift all subsystems by 1
    for (idx i = 0; i < n; ++i)
        perm.push_back((i + 1) % n);
    cout << ">> Subsytem(s): ";
    cout << disp(perm, ", ") << endl;
    t.tic();
    syspermute(randcmat, perm);
    cout << ">> It took " << t.toc() << " seconds." << endl << endl;
}

int main()
{
    // Examples
    MEASUREMENTS();
    TELEPORTATION();
    DENSE_CODING();
    GROVER();
    ENTANGLEMENT();
    QECC();
    CHANNEL();
    FUNCTOR();
    GRAMSCHMIDT();
    SPECTRAL();
    RANDOM();
    ENTROPIES();
    GRAPHSTATES();
    TIMING();
}
