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

// TODO: check for roundoff errors when using double->std::size_t conversions
// TODO: check the entropy functions, use svlas instead of hevals
// TODO: use .rows() instead of .cols() whenever possible

#include "qpp.h"

// #include "MATLAB/matlab.h" // support for MATLAB

using namespace qpp;

// We define each example as an independent lambda function

auto MEASUREMENTS = []
{
    std::cout << "**** Measurements ****" << std::endl;

    ket psi = st.b00;
    std::vector<std::size_t> subsys = {0};
    std::size_t result;
    std::vector<double> probs;
    std::vector<cmat> states;

    // measures the first subsystem of the Bell state (|00> + |11>)/sqrt(2)
    // in the X basis
    std::tie(result, probs, states) = measure(psi, gt.H, {0});
    std::cout << ">> Measuring part " << disp(subsys, " ")
            << " of the state: " << std::endl;
    std::cout << disp(psi) << std::endl;
    std::cout << ">> Measurement result: " << result << std::endl;
    std::cout << ">> Probabilities: " << disp(probs, ", ") << std::endl;
    std::cout << ">> Resulting normalized post-measurement states: "
            << std::endl;

    for (auto&& it: states)
        std::cout << disp(it) << std::endl << std::endl;

    // measure 2 subsystems out of a 4-qubit random density matrix
    cmat rho = randrho(16);
    subsys = {1, 2};
    cmat U = randU(4); // random basis on 2 qubits

    std::cout << ">> Measuring qubits " << disp(subsys, " ")
            << " of a 4-qubit random state in the random basis:" << std::endl;
    std::cout << disp(U) << std::endl;

    std::tie(result, probs, states) = measure(rho, U, {1, 2});
    std::cout << ">> Measurement result: " << result << std::endl;
    std::cout << ">> Probabilities: " << disp(probs, ", ") << std::endl;
    std::cout << ">> Sum of the probabilities: "
            << sum(probs.begin(), probs.end()) << std::endl;
    std::cout << ">> Resulting normalized post-measurement states: "
            << std::endl;

    for (auto&& it: states)
        std::cout << disp(it) << std::endl << std::endl;

    // check now how the state after the measurement "looks"
    // on the left over subsystems {0, 3}

    // it should be the same as the partial trace over {1, 2}
    // of the initial state (before the measurement), as local CPTP maps
    // do not influence the complementary subsystems

    cmat rho_bar = ptrace(rho, subsys);
    cmat rho_out_bar = cmat::Zero(4, 4);

    // compute the resulting mixed state after the measurement
    for (std::size_t i = 0; i < probs.size(); ++i)
        rho_out_bar += probs[i] * states[i];

    // verification
    std::cout << ">> Norm difference: " << norm(rho_bar - rho_out_bar)
            << std::endl;

    // random Kraus
    std::cout << ">> Random channel on part of the state " << std::endl;
    rho_bar = ptrace(rho, {1, 3});
    rho_out_bar = ptrace(channel(rho, randkraus(3, 4), {1, 3}), {1, 3});

    // verification
    std::cout << ">> Norm difference: " << norm(rho_bar - rho_out_bar)
            << std::endl << std::endl;
};

auto TELEPORTATION = []
{
    std::size_t D = 3; // size of the system
    std::cout << "**** Qudit teleportation, D = " << D << " ****" << std::endl;

    ket mes_AB = ket::Zero(D * D); // maximally entangled state resource
    for (std::size_t i = 0; i < D; ++i)
        mes_AB += mket({i, i}, D);
    mes_AB /= std::sqrt((double) D);

    // circuit that measures in the qudit Bell basis
    cmat Bell_aA = adjoint(gt.CTRL(gt.Xd(D), {0}, {1}, 2, D)
            * kron(gt.Fd(D), gt.Id(D)));

    ket psi_a = randket(D); // random state as input on a
    std::cout << ">> Initial state:" << std::endl;
    std::cout << disp(psi_a) << std::endl;

    ket input_aAB = kron(psi_a, mes_AB); // joint input state aAB
    // output before measurement
    ket output_aAB = apply(input_aAB, Bell_aA, {0, 1}, D);

    // measure on aA
    auto measured_aA = measure(output_aAB, gt.Id(D * D), {0, 1}, D);
    std::size_t m = std::get<0>(measured_aA); // measurement result

    auto midx = n2multiidx(m, {D, D});
    std::cout << ">> Alice's measurement result: ";
    std::cout << m << " -> " << disp(midx, " ") << std::endl;
    std::cout << ">> Alice's measurement probabilities: ";
    std::cout << disp(std::get<1>(measured_aA), ", ") << std::endl;

    // conditional result on B before correction
    cmat output_m_B = std::get<2>(measured_aA)[m];
    cmat correction_B = powm(gt.Zd(D), midx[0]) *
            powm(adjoint(gt.Xd(D)), midx[1]); // correction operator
    // apply correction on B
    std::cout << ">> Bob must apply the correction operator Z^" << midx[0]
            << " X^" << D - midx[1] << std::endl;
    cmat rho_B = correction_B * output_m_B * adjoint(correction_B);

    std::cout << ">> Bob's final state (after correction): " << std::endl;
    std::cout << disp(rho_B) << std::endl;

    // verification
    std::cout << ">> Norm difference: " << norm(rho_B - prj(psi_a))
            << std::endl << std::endl;
};

auto DENSE_CODING = []
{
    std::size_t D = 3; // size of the system
    std::cout << "**** Qudit dense coding, D = " << D << " ****" << std::endl;

    ket mes_AB = ket::Zero(D * D); // maximally entangled state resource
    for (std::size_t i = 0; i < D; ++i)
        mes_AB += mket({i, i}, D);
    mes_AB /= std::sqrt((double) D);

    // circuit that measures in the qudit Bell basis
    cmat Bell_AB = adjoint(gt.CTRL(gt.Xd(D), {0}, {1}, 2, D)
            * kron(gt.Fd(D), gt.Id(D)));

    // equal probabilities of choosing a message
    std::size_t m_A = randint(0, D * D - 1);
    auto midx = n2multiidx(m_A, {D, D});
    std::cout << ">> Alice sent: " << m_A << " -> ";
    std::cout << disp(midx, " ") << std::endl;

    // Alice's operation
    cmat U_A = powm(gt.Zd(D), midx[0]) * powm(adjoint(gt.Xd(D)), midx[1]);

    // Alice encodes the message
    ket psi_AB = apply(mes_AB, U_A, {0}, D);

    // Bob measures the joint system in the qudit Bell basis
    psi_AB = apply(psi_AB, Bell_AB, {0, 1}, D);
    auto measured = measure(psi_AB, gt.Id(D * D));

    std::cout << ">> Bob's measurement probabilities: ";
    std::cout << disp(std::get<1>(measured), ", ") << std::endl;

    // Bob samples according to the measurement probabilities
    std::size_t m_B = std::get<0>(measured);
    std::cout << ">> Bob received: ";
    std::cout << m_B << " -> " << disp(n2multiidx(m_B, {D, D}), " ")
            << std::endl << std::endl;
};

auto GROVER = []
{
    Timer t; // set a timer

    std::size_t n = 4; // number of qubits
    std::cout << "**** Grover on n = " << n << " qubits ****" << std::endl;

    std::vector<std::size_t> dims(n, 2); // local dimensions
    std::size_t N = std::pow(2, n); // number of elements in the database
    std::cout << ">> Database size: " << N << std::endl;

    // mark an element randomly
    std::size_t marked = randint(0, N - 1);
    std::cout << ">> Marked state: " << marked << " -> ";
    std::cout << disp(n2multiidx(marked, dims), " ") << std::endl;

    ket psi = mket(n2multiidx(0, dims)); // computational |0>^\otimes n

    // apply H^\otimes n, no aliasing
    psi = (kronpow(gt.H, n) * psi).eval();

    cmat G = 2 * prj(psi) - gt.Id(N); // Diffusion operator

    // number of queries
    std::size_t nqueries = std::ceil(pi * std::sqrt((double) N) / 4.);
    std::cout << ">> We run " << nqueries << " queries" << std::endl;
    for (std::size_t i = 0; i < nqueries; ++i)
    {
        psi(marked) = -psi(marked); // apply the oracle first, no aliasing
        psi = (G * psi).eval(); // then the diffusion operator, no aliasing
    }

    // we now measure the state in the computational basis
    auto measured = measure(psi, gt.Id(N));
    std::cout << ">> Probability of the marked state: "
            << std::get<1>(measured)[marked] << std::endl;
    std::cout << ">> Probability of all results: ";
    std::cout << disp(std::get<1>(measured), ", ") << std::endl;

    // sample
    std::cout << ">> Let's sample..." << std::endl;
    std::size_t result = std::get<0>(measured);
    if (result == marked)
        std::cout << ">> Hooray, we obtained the correct result: ";
    else
        std::cout << ">> Not there yet... we obtained: ";
    std::cout << result << " -> ";
    std::cout << disp(n2multiidx(result, dims), " ") << std::endl;

    // stop the timer and display it
    std::cout << ">> It took " << t.toc() << " seconds to simulate Grover on "
            << n << " qubits." << std::endl << std::endl;
};

auto ENTANGLEMENT = []
{
    std::cout << "**** Entanglement ****" << std::endl;

    cmat rho = 0.2 * st.pb00 + 0.8 * st.pb11;
    std::cout << ">> rho: " << std::endl;
    std::cout << disp(rho) << std::endl;

    std::cout << ">> Concurrence of rho: "
            << concurrence(rho) << std::endl;

    std::cout << ">> Negativity of rho: "
            << negativity(rho, {2, 2}) << std::endl;

    std::cout << ">> Logarithimc negativity of rho: "
            << lognegativity(rho, {2, 2}) << std::endl;

    ket psi = 0.8 * mket({0, 0}) + 0.6 * mket({1, 1});

    // apply some local random unitaries
    psi = kron(randU(2), randU(2)) * psi;

    std::cout << ">> psi: " << std::endl;
    std::cout << disp(psi) << std::endl;

    std::cout << ">> Entanglement of psi: "
            << entanglement(psi, {2, 2}) << std::endl;

    std::cout << ">> Concurrence of psi: "
            << concurrence(prj(psi)) << std::endl;

    std::cout << ">> G-Concurrence of psi: "
            << gconcurrence(psi) << std::endl;

    std::cout << ">> Schmidt coefficients of psi: " << std::endl;
    std::cout << disp(schmidtcoeff(psi, {2, 2})) << std::endl;

    std::cout << ">> Schmidt probabilities of psi: " << std::endl;
    std::cout << disp(schmidtprob(psi, {2, 2})) << std::endl;

    cmat U = schmidtU(psi, {2, 2});
    cmat V = schmidtV(psi, {2, 2});

    std::cout << ">> Schmidt vectors on Alice's side: " << std::endl;
    std::cout << disp(U) << std::endl;

    std::cout << ">> Schmidt vectors on Bob's side: " << std::endl;
    std::cout << disp(V) << std::endl;

    std::cout << ">> State psi in the Schmidt basis: " << std::endl;
    std::cout << disp(adjoint(kron(U, V)) * psi) << std::endl;

    // reconstructed state
    ket psi_from_schmidt =
            schmidtcoeff(psi, {2, 2})(0) * kron(U.col(0), V.col(0))
                    + schmidtcoeff(psi, {2, 2})(1) * kron(U.col(1), V.col(1));
    std::cout << ">> State psi reconstructed from the Schmidt decomposition: "
            << std::endl;
    std::cout << disp(psi_from_schmidt) << std::endl;

    // verification
    std::cout << ">> Norm difference: " << norm(psi - psi_from_schmidt)
            << std::endl << std::endl;
};

auto QECC = []
{
    std::cout << "**** Quantum error correcting codes ****" << std::endl;

    ket a0 = codes.codeword(Codes::Type::FIVE_QUBIT, 0);
    ket a1 = codes.codeword(Codes::Type::FIVE_QUBIT, 1);

    ket b0 = codes.codeword(Codes::Type::SEVEN_QUBIT_STEANE, 0);
    ket b1 = codes.codeword(Codes::Type::SEVEN_QUBIT_STEANE, 1);

    ket c0 = codes.codeword(Codes::Type::NINE_QUBIT_SHOR, 0);
    ket c1 = codes.codeword(Codes::Type::NINE_QUBIT_SHOR, 1);

    std::cout << ">> [[5, 1, 3]] Five qubit code. ";
    std::cout << "Checking codeword orthogonality." << std::endl;
    std::cout << ">> |<0L|1L>| = ";
    std::cout << disp(adjoint(a0) * a1) << std::endl;

    std::cout << ">> [[7, 1, 3]] Seven qubit Steane code. ";
    std::cout << "Checking codeword orthogonality." << std::endl;
    std::cout << ">> |<0L|1L>| = ";
    std::cout << disp(adjoint(b0) * b1) << std::endl;

    std::cout << ">> [[9, 1, 3]] Nine qubit Shor code. ";
    std::cout << "Checking codeword orthogonality." << std::endl;
    std::cout << ">> |<0L|1L>| = ";
    std::cout << disp(adjoint(c0) * c1) << std::endl << std::endl;
};

auto CHANNEL = []
{
    std::cout << "**** Channel tests ****" << std::endl;

    std::size_t nk = 5;
    std::size_t D = 3; // nk Kraus on d-dimensional system
    std::cout << ">> Generating a random channel with " << nk
            << " Kraus operators on a " << D << " dimensional space"
            << std::endl;
    std::vector<cmat> Ks = randkraus(nk, D);

    cmat rho_in = randrho(D); // random input state
    cmat rho_out = channel(rho_in, Ks); // output state

    std::cout << ">> Computing its Choi matrix..." << std::endl;
    cmat choim = choi(Ks);
    std::cout << ">> Choi matrix:" << std::endl << disp(choim) << std::endl;

    std::cout << ">> The eigenvalues of the Choi matrix are: "
            << std::endl << disp(transpose(hevals(choim))) << std::endl;

    std::cout << ">> Their sum is: " << sum(hevals(choim))
            << std::endl;

    std::vector<cmat> Kperps = choi2kraus(choim);
    std::cout << ">> The Kraus rank of the channel is: "
            << Kperps.size() << std::endl;

    cmat rho_out1 = channel(rho_in, Kperps);
    // verification
    std::cout << ">> Norm difference on output states: "
            << norm(rho_out1 - rho_out) << std::endl;

    std::cout << ">> Superoperator matrix:" << std::endl;
    cmat smat = super(Ks);
    std::cout << disp(smat) << std::endl;

    std::cout << ">> The eigenvalues of the superoperator matrix are: "
            << std::endl;
    cmat evalsupop = evals(smat);
    std::cout << disp(transpose(evalsupop)) << std::endl;

    std::cout << ">> Their absolute values are: " << std::endl;
    for (std::size_t i = 0; i < (std::size_t) evalsupop.size(); i++)
        std::cout << std::abs(evalsupop(i)) << " ";

    // verification
    std::cout << std::endl
            << ">> Diference in norm for superoperator action: ";
    cmat rho_out2 = transpose(
            reshape(smat * reshape(transpose(rho_in), D * D, 1), D, D));
    std::cout << norm(rho_out - rho_out2) << std::endl << std::endl;
};

cplx pow3(const cplx& z) // test function used by qpp::cwise() in FUNCTOR()
{
    return std::pow(z, 3);
}

auto FUNCTOR = []
{
    std::cout << "**** Functor ****" << std::endl;

    // functor test
    std::cout << ">> Functor z^3 acting component-wise on:" << std::endl;
    cmat A(2, 2);
    A << 1, 2, 3, 4;
    std::cout << disp(A) << std::endl;

    std::cout << ">> Result (with lambda):" << std::endl;
    // functor z^3 componentwise, specify OutputScalar and Derived for lambdas
    std::cout << disp(cwise<cplx, cmat>(A, [](const cplx& z) -> cplx
    {
        return z * z * z;
    })) << std::endl;

    std::cout << ">> Result (with proper function):" << std::endl;
    // automatic type deduction for proper functions
    std::cout << disp(cwise(A, &pow3)) << std::endl << std::endl;
};

auto GRAMSCHMIDT = []
{
    std::cout << "**** Gram-Schmidt ****" << std::endl;

    cmat A(3, 3);
    A << 1, 1, 0, 0, 2, 0, 0, 0, 0;
    std::cout << ">> Input matrix:" << std::endl << disp(A) << std::endl;

    cmat Ags = grams(A);
    std::cout << ">> Result:" << std::endl << disp(Ags) << std::endl;

    std::cout << ">> Projector onto G.S. vectors:" << std::endl;
    std::cout << disp(Ags * adjoint(Ags)) << std::endl << std::endl;
};

auto SPECTRAL = []
{
    std::cout << "**** Spectral decomposition tests ****" << std::endl;
    std::size_t D = 4;
    cmat rH = randH(D);
    std::cout << ">> Original matrix: " << std::endl << disp(rH) << std::endl;

    // spectral decomposition here
    DynColVect<double> evalsH = hevals(rH);
    cmat evectsH = hevects(rH);
    cmat spec = cmat::Zero(D, D);
    // reconstruct the matrix
    for (std::size_t i = 0; i < D; i++)
        spec += evalsH(i) * prj(evectsH.col(i));

    std::cout << ">> Reconstructed from spectral decomposition: " << std::endl;
    std::cout << disp(spec) << std::endl;

    // verification
    std::cout << ">> Norm difference: " << norm(spec - rH)
            << std::endl << std::endl;
};

auto RANDOM = []
{
    std::cout << "**** Randomness ****" << std::endl;

    std::cout << ">> Generating a random ket on D = 5" << std::endl;
    ket rket = randket(5);
    std::cout << disp(rket) << std::endl;

    std::vector<double> probs = abssq(rket);
    std::cout << "Probabilities: " << disp(probs, ", ") << std::endl;

    std::cout << "Sum of the probabilities: ";
    std::cout << sum(probs.begin(), probs.end()) << std::endl << std::endl;
};

auto TIMING = []
{
    std::cout << "**** Timing tests ****" << std::endl;

    std::size_t n = 12; // number of qubits
    std::size_t N = std::pow(2, n);
    std::cout << ">> n = " << n << " qubits, matrix size "
            << N << " x " << N << "." << std::endl << std::endl;
    cmat randcmat = cmat::Random(N, N);

    // qpp::ptrace()
    std::cout << "**** qpp::ptrace() timing ****" << std::endl;
    std::vector<std::size_t> subsys_ptrace = {0};
    std::cout << ">> Subsytem(s): ";
    std::cout << disp(subsys_ptrace, ", ") << std::endl;
    Timer t;
    ptrace(randcmat, subsys_ptrace);
    std::cout << ">> It took " << t.toc() << " seconds."
            << std::endl << std::endl;

    // qpp::ptranspose()
    std::cout << "**** qpp::ptranspose() timing ****" << std::endl;
    // partially transpose n-1 subsystems
    std::vector<std::size_t> subsys_ptranspose;
    for (std::size_t i = 0; i < n - 1; ++i)
        subsys_ptranspose.push_back(i);
    std::cout << ">> Subsytem(s): ";
    std::cout << disp(subsys_ptranspose, ", ") << std::endl;
    t.tic();
    ptranspose(randcmat, subsys_ptranspose);
    std::cout << ">> It took " << t.toc() << " seconds."
            << std::endl << std::endl;

    // qpp::syspermute()
    std::cout << "**** qpp::syspermute() timing ****" << std::endl;
    std::vector<std::size_t> perm; // left-shift all subsystems by 1
    for (std::size_t i = 0; i < n; ++i)
        perm.push_back((i + 1) % n);
    std::cout << ">> Subsytem(s): ";
    std::cout << disp(perm, ", ") << std::endl;
    t.tic();
    syspermute(randcmat, perm);
    std::cout << ">> It took " << t.toc() << " seconds."
            << std::endl << std::endl;
};

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
    TIMING();
}
