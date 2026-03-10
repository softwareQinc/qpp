// Source: ./examples/circuits/hhl_circuit.cpp
//
// HHL algorithm: solving a 2x2 linear system A x = b
//
// Matrix:
// A =  |1  0|
//      |0  2|
//
// Right-hand side (encoded as a quantum state):
//
// |b> = |+> = (|0> + |1>)/\sqrt{2}
//
// The classical solution is:
//
// x = A^{-1} b
//
// Since A^{-1} = diag(1, 1/2),
//
// A^{-1}|b> = (1|0> + (1/2)|1>)/\sqrt{2}
//
// After normalization the expected quantum state is:
//
// |x> = (|0> + 0.5|1>) / \sqrt{1 + 0.25}
//     ~= 0.894427 |0> + 0.447214 |1>
//
// The HHL algorithm prepares a quantum state proportional to A^{-1}|b>
// (up to normalization) using:
//
//   1. Quantum Phase Estimation on U = exp(iAt)
//   2. Controlled rotation proportional to 1/\lambda
//   3. Uncomputation of the phase register
//   4. Postselection on an ancilla qubit
//
// The resulting system qubit should approximate |x>.

#include <cmath>
#include <complex>
#include <iostream>
#include <vector>

#include <qpp/qpp.hpp>

int main() {
    using namespace qpp;

    std::cout << ">> HHL algorithm (QCircuit version)\n";
    std::cout << ">> Solving A|x> = |b>\n";
    std::cout << ">> A = diag(1, 2) and |b> = (1, 1)^T/\\sqrt{2}\n";
    std::cout << ">> Expected x = (0.894427, 0.447214)^T\n";

    // registers
    const idx n_phase = 2;         // phase qubits: 0, 1
    const idx q_sys = n_phase;     // system qubit: 2
    const idx q_anc = n_phase + 1; // ancilla qubit: 3
    const idx n_total = 4;         // total number of qubits

    const realT t = pi / 2.0; // QPE time parameter
    const realT C = 0.5;      // Scaling constant for rotation

    // base U = exp(i A t) for A = diag(1,2) => diag(e^{it}, e^{i2t})
    cmat Ubase = cmat::Zero(2, 2);
    Ubase(0, 0) = std::exp(cplx(0, 1) * (1.0 * t));
    Ubase(1, 1) = std::exp(cplx(0, 1) * (2.0 * t));

    // QCircuit construction
    QCircuit qc{n_total};

    // prepare |b> on system qubit: |+> = H|0>
    qc.gate(gt.H, q_sys);

    // quantum Phase Estimation (QPE)
    for (idx i = 0; i < n_phase; ++i) {
        qc.gate(gt.H, i);
    }
    cmat U = Ubase;
    for (idx k = 0; k < n_phase; ++k) {
        const idx ctrl = n_phase - k - 1;
        qc.CTRL(U, ctrl, q_sys);
        U = powm(U, 2);
    }

    // inverse QFT on phase register
    qc.TFQ({0, 1});

    // controlled rotation on ancilla
    auto Ry = [&](realT theta) -> cmat {
        cmat R(2, 2);
        const realT c = std::cos(theta / 2.0);
        const realT s = std::sin(theta / 2.0);
        R << c, -s, s, c;
        return R;
    };

    // phase value -> lambda (for this choice)
    auto lambda_of_v = [&](idx v) -> realT {
        if (v == 1) {
            return 1.0; // Phase "01"
        }
        if (v == 2) {
            return 2.0; // Phase "10"
        }
        return 0.0;
    };

    // we build a multi-qubit gate G that acts on {phase0, phase1, ancilla}
    cmat G = cmat::Identity(8, 8);
    for (idx v = 0; v < 4; ++v) {
        const realT lam = lambda_of_v(v);
        if (lam > 0.0) {
            const realT theta = 2.0 * std::asin(C / lam);
            const idx base = 2 * v;
            G.block(base, base, 2, 2) = Ry(theta);
        }
    }
    qc.gate(G, {0, 1, q_anc});

    // uncompute QPE
    qc.QFT({0, 1});
    U = Ubase;
    std::vector<cmat> powers;
    for (idx k = 0; k < n_phase; ++k) {
        powers.push_back(U);
        U = powm(U, 2);
    }
    for (idx k = n_phase; k-- > 0;) {
        const idx ctrl = n_phase - k - 1;
        qc.CTRL(adjoint(powers[k]), ctrl, q_sys);
    }
    for (idx i = 0; i < n_phase; ++i) {
        qc.gate(gt.H, i);
    }

    // execution
    QEngine engine{qc};
    engine.execute();
    ket psi = engine.get_state();

    // postselect ancilla = |1>
    const cmat P1 = 1_prj;
    psi = apply(psi, P1, {q_anc});

    const realT p_success = psi.squaredNorm();
    if (p_success < 1e-15) {
        std::cout
            << ">> Ancilla success probability ~ 0; something went wrong.\n";
        return 1;
    }
    psi /= std::sqrt(p_success);

    // reduce to system density matrix (trace out {0, 1, 3})
    const cmat rho_sys = ptrace(psi, {0, 1, q_anc});

    // classical normalized solution: A^{-1}|b> ∝ (1*|0> + 1/2*|1>)
    ket x_class(2);
    x_class << 1.0, 0.5;
    x_class /= x_class.norm();

    // fidelity = <x_class|rho_sys|x_class>
    realT fidelity = (adjoint(x_class) * rho_sys * x_class)(0, 0).real();

    // optionally extract a representative |x> from rho_sys (dominant
    // eigenvector)
    auto [evals, evects] = heig(rho_sys);
    Eigen::SelfAdjointEigenSolver<cmat> es(rho_sys);
    idx imax = (evals(1) > evals(0)) ? 1 : 0;
    ket x_quant = evects.col(imax);
    x_quant /= x_quant.norm();

    auto dirac_disp_opt = IOManipDiracOpts{}
                              .set_amplitudes_after(false)
                              .set_plus_op(" + ")
                              .set_mul_op("");
    std::cout << ">> Postselection success probability p = " << p_success
              << "\n";
    std::cout << ">> rho_sys (system reduced density matrix):\n"
              << disp(rho_sys) << "\n";
    std::cout << ">> Quantum |x> (dominant eigenvector of rho_sys):\n"
              << disp(dirac(x_quant), dirac_disp_opt) << "\n";
    std::cout << ">> Classical x (normalized):\n" << disp(x_class) << "\n";
    std::cout << ">> Fidelity = " << fidelity << "\n";
}
