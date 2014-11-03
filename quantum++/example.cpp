/*
 * File:   example.cpp
 * Author: vlad
 *
 * Created on December 12, 2013, 10:38 PM
 */

#include "qpp.h"

// #include "MATLAB/matlab.h" // support for MATLAB

// TODO: testing (IN PROGRESS)
// TODO: parallelize everything (Gates::CTRL etc) - IMPORTANT! (IN PROGRESS)
// TODO: change loops to column major order (IN PROGRESS)
// TODO: inline small functions (IN PROGRESS)

// TODO: write introductory page (documentation)
// TODO: more robust exception parameter checking

using namespace qpp;

cplx pow3(const cplx& z) // a test function
{
	return std::pow(z, 3);
}

int main()
{
	std::discrete_distribution<std::size_t> dd{0.1,0.9};
	std::cout << dd(rdevs._rng); // sample
	std::cout << dd(rdevs._rng); // sample again



	/*
	// TESTING

	// testing applyCTRL
	size_t d = 3;
	size_t n = 3;
	std::vector<std::size_t> ctrl = { };
	std::vector<std::size_t> subsys = { 2 };
	ket psi = mket( { 1, 0, 1 }, d);
	cmat rho = mprj( { 1, 0, 1 }, d);
	cmat gate = gt.Xd(3);
	displn(psi);
	std::cout << std::endl;
	auto result_ket = applyCTRL(psi, gate, ctrl, subsys, n, d);
	displn(result_ket);
	std::cout << std::endl;
	std::cout << std::endl;
	auto result_rho = applyCTRL(rho, gate, ctrl, subsys, n, d);
	displn(result_rho);
	std::cout << std::endl << "Norm difference: ";
	displn(norm(result_rho - result_ket * adjoint(result_ket)));

	// testing channel and apply
	std::cout << std::endl << "Testing channel(...) and apply(...)."
			<< std::endl;
	rho = randrho(16);
	cmat K = kron(gt.Id2, gt.X, gt.Y, gt.Z);
	std::vector<std::size_t> p = randperm(4); // permutation
	std::cout << "Permutation: ";
	displn(p, ", ");
	std::vector<std::size_t> invp = invperm(p); // inverse permutation
	std::cout << "Inverse permutation: ";
	displn(invp, ", ");
	cmat r1 = channel(rho, { K }, p, 4, 2);
	cmat r2 = syspermute(channel(syspermute(rho, p, { 2, 2, 2, 2 }), { K }, { 0,
			1, 2, 3 }, 4, 2), invp, { 2, 2, 2, 2 });
	std::cout << norm(r1 - r2) << std::endl << std::endl;

	r1 = apply(rho, K, p, 4, 2);
	r2 = syspermute(
			apply(syspermute(rho, p, { 2, 2, 2, 2 }), K, { 0, 1, 2, 3 }, 4, 2),
			invp, { 2, 2, 2, 2 });
	std::cout << norm(r1 - r2) << std::endl << std::endl;

	displn(channel(prj(mket( { 0, 1 })), { gt.CNOTab }, { 1, 0 }, 2, 2));
	std::cout << std::endl;
	displn(apply(mket( { 0, 0 }), gt.CNOTab, { 0, 1 }, 2, 2));

	// quantum teleportation
	std::cout << std::endl << "Qudit teleportation." << std::endl;
	psi = randket(2); // a random state;
	std::cout << "|psi><psi|:" << std::endl;
	displn(prj(psi));
	cmat telecircuit = gt.expandout(gt.H, { 0 }, { 2, 2, 2 }) * gt.CTRL(gt.X, {
			0 }, { 1 }, 3);
	ket psiin = kron(psi, st.b00); // input state
	ket psiout = telecircuit * psiin; // output state before measurement
	// measure Alice's qubits, measurement results are 1 0
	psiout = kron(prj(st.z1), prj(st.z0), gt.Id2) * psiout;
	// apply correction
	psiout = gt.expandout(powm(gt.Z, 1) * powm(gt.X, 0), { 2 }, { 2, 2, 2 })
			* psiout;
	// not necessary to normalize, prj() takes care of it below
	cmat rhoout = ptrace(prj(psiout), { 0, 1 }, { 2, 2, 2 });
	std::cout << std::endl << "Teleported state:" << std::endl;
	displn(rhoout);
	std::cout << "Difference in norm: " << norm(prj(psi) - rhoout) << std::endl;

	// qudit measurements
	std::cout << std::endl << "Qudit measurements." << std::endl;
	std::cout << "Initially in state |0><0|." << std::endl;
	ket zd0(3);
	zd0 << 1, 0, 0;
	experimental::Qudit q(prj(zd0));
	std::cout << "Measuring Z operator non-destructively. Results:"
			<< std::endl;
	std::cout << q.measure() << std::endl;
	std::cout << q.measure() << std::endl;
	std::cout << q.measure() << std::endl;
	std::cout << "Measuring X operator non-destructively. Results:"
			<< std::endl;
	std::cout << q.measure(gt.Xd(3)) << std::endl;
	std::cout << q.measure(gt.Xd(3)) << std::endl;
	std::cout << q.measure(gt.Xd(3)) << std::endl;
	// von Neumann projective measurement
	std::cout << "Measuring X operator destructively (collapse). Results:"
			<< std::endl;
	std::cout << q.measure(gt.Xd(3), true) << std::endl;
	std::cout << q.measure(gt.Xd(3)) << std::endl;
	std::cout << q.measure(gt.Xd(3)) << std::endl;
	std::cout << "Finally measuring Z operator destructively. Results:"
			<< std::endl;
	std::cout << q.measure(true) << std::endl;
	std::cout << q.measure() << std::endl;
	std::cout << q.measure() << std::endl;
	std::cout << "Final state of qudit:" << std::endl;
	displn(q.getRho());

	// Bell state generator
	std::cout << std::endl << "Bell state generator: " << std::endl;
	cmat circuit;
	circuit = gt.CTRL(gt.X, { 0 }, { 1 }, 2) * gt.expandout(gt.H, 0, { 2, 2 });
	cmat input = kron(st.z0, st.z0);
	cmat output = circuit * input;
	std::cout << "Circuit matrix representation: " << std::endl;
	displn(circuit);
	std::cout << std::endl << "Output (|Bell_0> state) of the circuit on |00>: "
			<< std::endl;
	displn(output);

	// 3-qubit repetion code
	std::cout << std::endl << "3-qubit repetition code: " << std::endl;
	cmat rep;
	rep = gt.CTRL(gt.X, { 0 }, { 2 }, 3) * gt.CTRL(gt.X, { 0 }, { 1 }, 3);
	input = kron(st.z1, st.z0, st.z0);
	output = rep * input;
	std::cout << "Circuit acting on |000> produces |111>. Check: " << std::endl;
	displn(output);

	// functor test
	std::cout << std::endl << "Functor z^3 acting on:" << std::endl;
	cmat a(2, 2);
	a << 1, 2, 3, 4;
	displn(a);
	std::cout << "Result (with lambda):" << std::endl;
	// functor z^3 componentwise, specify OutputScalar and Derived for lambdas
	displn(cwise<cplx, cmat>(a, [](const cplx& z)->cplx
	{	return z*z*z;}));
	std::cout << "Result (with proper function):" << std::endl;
	// automatic type deduction for proper functions
	displn(cwise(a, &pow3));

	// Gram-Schmidt
	std::cout << std::endl << "Gram-Schmidt on matrix:" << std::endl;
	cmat A(3, 3);
	A << 1, 1, 0, 0, 2, 0, 0, 0, 0;
	displn(A);
	cmat Ags = grams(A);
	std::cout << std::endl << "Result:" << std::endl;
	displn(Ags);
	std::cout << std::endl << "Projector is:" << std::endl;
	displn(Ags * adjoint(Ags));

	// spectral decomposition test
	std::cout << std::endl << "Spectral decomposition tests." << std::endl;
	std::size_t D = 4;
	cmat rH = randH(D);
	DynColVect<double> evalsH = hevals(rH);
	cmat evectsH = hevects(rH);
	cmat spec = cmat::Zero(D, D);
	for (std::size_t i = 0; i < D; i++)
		spec += evalsH(i) * prj(evectsH.col(i));
	std::cout << "Original matrix: " << std::endl;
	displn(rH);
	std::cout << std::endl << "Reconstructed from spectral decomposition: "
			<< std::endl;
	displn(spec);
	std::cout << "Difference in norm: " << norm(spec - rH) << std::endl;

	// channel tests
	std::cout << std::endl << "Channel tests." << std::endl;
	std::size_t nk = 5;
	d = 3; // nk Kraus on d-dimensional system
	std::cout << "Generating a random channel with " << nk
			<< " Kraus operators on a " << d << " dimensional space..."
			<< std::endl;
	std::vector<cmat> Ks = randkraus(nk, d);

	cmat rho_in = randrho(d); // input state
	cmat rho_out = channel(rho_in, Ks); // output state

	std::cout << "Computing its Choi matrix..." << std::endl;
	cmat choim = choi(Ks);
	std::cout << "Choi matrix:" << std::endl;
	displn(choim);
	std::cout << std::endl << "The eigenvalues of the Choi matrix are: "
			<< std::endl;
	displn(transpose(hevals(choim)));
	std::cout << std::endl << "Their sum is: " << sum(hevals(choim))
			<< std::endl;
	std::vector<cmat> Kperps = choi2kraus(choim);
	std::cout << std::endl << "The Kraus rank of the channel is: "
			<< Kperps.size() << std::endl;
	cmat rho_out1 = channel(rho_in, Kperps);
	std::cout << std::endl << "Difference in norm on output states: "
			<< norm(rho_out1 - rho_out) << std::endl;
	std::cout << std::endl << "Superoperator matrix:" << std::endl;

	Timer tmp;

	cmat smat1 = super(Ks);
	tmp.toc();
	std::clog << tmp << std::endl;

	tmp.tic();
	cmat smat = experimental::super(Ks);
	tmp.toc();
	std::clog << tmp << std::endl;

	std::clog << norm(smat - smat1) << std::endl;

	displn(smat);
	std::cout << std::endl
			<< "The eigenvalues of the superoperator matrix are: " << std::endl;
	cmat evalsupop = evals(smat);
	displn(transpose(evalsupop));
	std::cout << std::endl << "Their absolute values are: " << std::endl;
	for (std::size_t i = 0; i < (std::size_t) evalsupop.size(); i++)
		std::cout << std::abs((cplx) evalsupop(i)) << " ";
	std::cout << std::endl << std::endl
			<< "Diference in norm for superoperator action: ";
	cmat rho_out2 = transpose(
			reshape(smat * reshape(transpose(rho_in), d * d, 1), d, d));
	std::cout << norm(rho_out - rho_out2) << std::endl;

	// statistics tests
	std::cout << std::endl << "Statistics tests." << std::endl;
	std::vector<cplx> ampl = { 1. + 1_i, 1. - 1_i };
	ket va(4);
	va << 0.1, 1, 1. + 1_i, 1. + 2_i;
	std::vector<double> weights = amplitudes(std::begin(ampl), std::end(ampl));
	std::discrete_distribution<std::size_t> dc(std::begin(weights),
			std::end(weights));

	std::cout << "The probabilities are: ";
	std::vector<double> probs = dc.probabilities();
	displn(probs, ", ", "{", "}");
	std::cout << "Their sum is: " << sum(std::begin(probs), std::end(probs))
			<< std::endl;

	// 	// TIMING tests
	std::cout << std::endl << "Timing tests..." << std::endl;
	n = 12; // number of qubits
	std::size_t N = std::pow(2, n);
	std::vector<std::size_t> dims(n, 2); // local dimensions
	std::cout << "n = " << n << " qubits, matrix size " << N << " x " << N
			<< "." << std::endl;

	// matrix initialization
	std::cout << std::endl << "Matrix initialization timing." << std::endl;
	// start the timer, automatic tic() in the constructor
	Timer t, total;
	cmat randcmat = cmat::Random(N, N);
	t.toc(); // read the time
	std::cout << "Took " << t << " seconds." << std::endl;

	// lazy matrix product
	std::cout << std::endl << "Lazy matrix product timing." << std::endl;
	t.tic();
	auto lazyprod = randcmat * randcmat; // lazyprod has type GenMatProduct
	t.toc(); // read the time
	std::cout << "Took " << t << " seconds." << std::endl;

	// ptrace1 timing
	std::cout << std::endl << "ptrace1 timing." << std::endl;
	t.tic(); // reset the chronometer
	// trace away half of the qubits
	ptrace1(randcmat,
			{ (std::size_t) std::sqrt(N), (std::size_t) std::sqrt(N) });
	t.toc(); // read the time
	std::cout << "Took " << t << " seconds." << std::endl;

	// ptrace2 timing
	std::cout << std::endl << "ptrace2 timing." << std::endl;
	t.tic(); // reset the chronometer
	// trace away half of the qubits
	ptrace2(randcmat,
			{ (std::size_t) std::sqrt(N), (std::size_t) std::sqrt(N) });
	t.toc(); // read the time
	std::cout << "Took " << t << " seconds." << std::endl;

	// ptrace
	std::cout << std::endl << "ptrace timing." << std::endl;
	std::vector<std::size_t> subsys_ptrace = { 0 };
	std::cout << "Subsytem(s): ";
	displn(subsys_ptrace, ", ");
	t.tic();
	ptrace(randcmat, subsys_ptrace, dims);
	t.toc();
	std::cout << "Took " << t << " seconds." << std::endl;

	// ptranspose
	std::cout << std::endl << "ptranspose timing." << std::endl;
	std::vector<std::size_t> subsys_ptranspose; // partially transpose n-1 subsystems
	for (std::size_t i = 0; i < n - 1; i++)
		subsys_ptranspose.push_back(i);
	std::cout << "Subsytem(s): ";
	displn(subsys_ptranspose, ", ");
	t.tic();
	ptranspose(randcmat, subsys_ptranspose, dims);
	t.toc();
	std::cout << "Took " << t << " seconds." << std::endl;

	// syspermute
	std::cout << std::endl << "syspermute timing." << std::endl;
	std::vector<std::size_t> perm; // left-shift all subsystems by 1
	for (std::size_t i = 0; i < n; i++)
		perm.push_back((i + 1) % n);
	std::cout << "Subsytem(s): ";
	displn(perm, ", ");
	t.tic();
	syspermute(randcmat, perm, dims);
	t.toc();
	std::cout << "Took " << t << " seconds." << std::endl;

	//	 // matrix product
	//	 std::cout << std::endl << "Matrix product timing." << std::endl;
	//	 t.tic(); // reset the chronometer
	//	 cmat prodmat = randcmat * randcmat; // explicit cmat now
	//	 t.toc(); // read the time
	//	 std::cout << "Took " << t << " seconds." << std::endl;

	// END TIMING
	total.toc(); // read the total running time
	std::cout << std::endl << "Total time: " << total.seconds() << " seconds.";
	//	 */
}
