/*
 * gates.h
 *
 *  Created on: Apr 7, 2014
 *      Author: vlad
 */

#ifndef GATES_H_
#define GATES_H_

namespace qpp
{

class Gates: public Singleton<const Gates> // const Singleton
{
	friend class Singleton<const Gates> ;
public:
	// one qubit gates
	types::cmat Id2 { types::cmat::Identity(2, 2) }; // Identity matrix
	types::cmat H { types::cmat::Zero(2, 2) }; // Hadamard matrix
	types::cmat X { types::cmat::Zero(2, 2) }; // X matrix
	types::cmat Y { types::cmat::Zero(2, 2) }; // Y matrix
	types::cmat Z { types::cmat::Zero(2, 2) }; // Z matrix
	types::cmat S { types::cmat::Zero(2, 2) }; // S gate
	types::cmat T { types::cmat::Zero(2, 2) }; // T gate

	// two qubit gates
	types::cmat CNOTab { types::cmat::Identity(4, 4) }; // CNOT ctrl1 target2
	types::cmat CZ { types::cmat::Identity(4, 4) }; // Controlled-Phase (Controlled-Z)
	types::cmat CNOTba { types::cmat::Zero(4, 4) }; // CNOT ctrl2 target1
	types::cmat SWAP { types::cmat::Identity(4, 4) }; // SWAP gate

	// three qubit gates
	types::cmat TOF { types::cmat::Identity(8, 8) }; // Toffoli
	types::cmat FRED { types::cmat::Identity(8, 8) }; // Fredkin
private:
	Gates()
	{
		// initialize the constants and gates
		H << 1 / std::sqrt(2.), 1 / std::sqrt(2.), 1 / std::sqrt(2.), -1
				/ std::sqrt(2.);
		X << 0, 1, 1, 0;
		Z << 1, 0, 0, -1;
		Y << 0, -1_i, 1_i, 0;
		S << 1, 0, 0, 1_i;
		T << 1, 0, 0, std::exp(1_i * ct::pi / 4.0);
		CNOTab.block(2, 2, 2, 2) = X;
		CNOTba(0, 0) = 1;
		CNOTba(1, 3) = 1;
		CNOTba(2, 2) = 1;
		CNOTba(3, 1) = 1;
		CZ(3, 3) = -1;

		SWAP.block(1, 1, 2, 2) = X;
		TOF.block(6, 6, 2, 2) = X;
		FRED.block(4, 4, 4, 4) = SWAP;
	}
public:
	// gates with variable dimension

	// one qubit gates

	// Rotation of theta about n (a unit vector {nx, ny, nz})
	types::cmat Rn(double theta, std::vector<double> n) const
	{
		if (n.size() != 3) // not a 3-D vector
			throw Exception("Gates::Rn", "n is not a 3-D vector!");

		types::cmat result(2, 2);
		result = std::cos(theta / 2) * Id2
				- 1_i * std::sin(theta / 2) * (n[0] * X + n[1] * Y + n[2] * Z);
		return result;
	}

	// one quDit gates

	types::cmat Zd(std::size_t D) const
	{
		if (D == 0)
			throw Exception("Gates::Zd", Exception::Type::DIMS_INVALID);

		types::cmat result(D, D);
		result = types::cmat::Zero(D, D);
		for (std::size_t i = 0; i < D; i++)
			result(i, i) = std::pow(ct::omega(D), i);
		return result;
	}

	types::cmat Fd(std::size_t D) const
	{
		if (D == 0)
			throw Exception("Gates::Fd", Exception::Type::DIMS_INVALID);

		types::cmat result(D, D);
		result = types::cmat::Zero(D, D);
		for (std::size_t j = 0; j < D; j++)
			for (std::size_t i = 0; i < D; i++)
				result(i, j) = 1 / std::sqrt((double)D) * std::pow(ct::omega(D), i * j);
		return result;
	}

	types::cmat Xd(std::size_t D) const // X|k>=|k+1>
	{
		if (D == 0)
			throw Exception("Gates::Xd", Exception::Type::DIMS_INVALID);

		return static_cast<types::cmat>(Fd(D).inverse() * Zd(D) * Fd(D));
	}

	template<typename Derived = Eigen::MatrixXcd>
	Derived Id(std::size_t D) const
	{
		if (D == 0)
			throw Exception("Gates::Id", Exception::Type::DIMS_INVALID);
		return Derived::Identity(D, D);
	}

	// applies controlled-gate A to part of state vector
	// or density matrix specified by subsys
	template<typename Derived1, typename Derived2>
	types::DynMat<typename Derived1::Scalar> applyCTRL(
			const Eigen::MatrixBase<Derived1>& state,
			const Eigen::MatrixBase<Derived2>& A,
			const std::vector<std::size_t>& ctrl,
			const std::vector<std::size_t>& subsys, std::size_t n,
			std::size_t d = 2) const
	{
	}

	// applies gate A to part of state vector
	// or density matrix specified by subsys
	template<typename Derived1, typename Derived2>
	types::DynMat<typename Derived1::Scalar> apply(
			const Eigen::MatrixBase<Derived1>& state,
			const Eigen::MatrixBase<Derived2>& A,
			const std::vector<std::size_t>& subsys,
			const std::vector<std::size_t>& dims) const
	{
		const types::DynMat<typename Derived1::Scalar> & rstate = state;
		const types::DynMat<typename Derived2::Scalar> & rA = A;

		// EXCEPTION CHECKS

		// check types
		if (!std::is_same<typename Derived1::Scalar, typename Derived2::Scalar>::value)
			throw Exception("Gates::apply", Exception::Type::TYPE_MISMATCH);

		// check zero sizes
		if (!internal::_check_nonzero_size(rA))
			throw Exception("Gates::apply", Exception::Type::ZERO_SIZE);

		// check zero sizes
		if (!internal::_check_nonzero_size(rstate))
			throw Exception("Gates::apply", Exception::Type::ZERO_SIZE);

		// check square matrix for the subsys
		if (!internal::_check_square_mat(rA))
			throw Exception("Gates::apply", Exception::Type::MATRIX_NOT_SQUARE);

		// check that dims is a valid dimension vector
		if (!internal::_check_dims(dims))
			throw Exception("Gates::apply", Exception::Type::DIMS_INVALID);

		// check subsys is valid w.r.t. dims
		if (!internal::_check_subsys_match_dims(subsys, dims))
			throw Exception("Gates::apply",
					Exception::Type::SUBSYS_MISMATCH_DIMS);

		// Use static allocation for speed!
		std::size_t Cdims[ct::maxn];
		std::size_t midx_row[ct::maxn];
		std::size_t midx_rho_row[ct::maxn];

		std::size_t CdimsA[ct::maxn];
		std::size_t CsubsysA[ct::maxn];
		std::size_t midxA_row[ct::maxn];
		std::size_t midxA_rho_row[ct::maxn];

		std::size_t CdimsA_bar[ct::maxn];
		std::size_t CsubsysA_bar[ct::maxn];
		std::size_t midxA_bar_row[ct::maxn];

		std::size_t n = dims.size();
		std::size_t nA = subsys.size();
		std::size_t nA_bar = n - nA;

		std::size_t D = 1;
		std::size_t DA_bar = 1;

		for (std::size_t k = 0, cnt = 0; k < n; k++)
		{
			midx_row[k] = midx_rho_row[k] = 0;
			Cdims[k] = dims[k];
			D *= dims[k];

			// compute the complement of subsys w.r.t. dims
			if (std::find(std::begin(subsys), std::end(subsys), k)
					== std::end(subsys))
			{
				CsubsysA_bar[cnt] = k;
				CdimsA_bar[cnt] = dims[k];
				midxA_bar_row[cnt] = 0;
				DA_bar *= dims[k];
				cnt++;
			}
		}

		std::size_t DA = 1;
		for (std::size_t k = 0; k < nA; k++)
		{
			midxA_row[k] = midxA_rho_row[k] = 0;
			CdimsA[k] = dims[subsys[k]];
			CsubsysA[k] = subsys[k];
			DA *= dims[subsys[k]];
		}

		if (internal::_check_col_vector(rstate)) // we have a ket
		{
			// check that dims match state vector
			if (!internal::_check_dims_match_cvect(dims, rstate))
				throw Exception("Gates::apply",
						Exception::Type::DIMS_MISMATCH_CVECTOR);

			// check that state vector matches the dimensions of the subsys
			if (static_cast<std::size_t>(rA.cols()) != DA)
				throw Exception("Gates::apply",
						Exception::Type::DIMS_MISMATCH_CVECTOR);

			types::DynMat<typename Derived1::Scalar> result(D, 1);

			// run over the subsys's row multi-index
			for (std::size_t a = 0; a < DA; a++)
			{
				// get the subsys's row multi-index
				internal::_n2multiidx(a, nA, CdimsA, midxA_row);
				// compute subsys part of the result's row multi-index
				for (std::size_t k = 0; k < nA; k++)
					midx_row[CsubsysA[k]] = midxA_row[k];

				// run over the complement's row multi-index
				for (std::size_t i = 0; i < DA_bar; i++)
				{
					// get the complement's row multi-index
					internal::_n2multiidx(i, nA_bar, CdimsA_bar, midxA_bar_row);
					// now compute the complement part of the
					// result's row multi-index
					// and the complement part
					// of the state's total row multi-index
					for (std::size_t k = 0; k < nA_bar; k++)
						midx_row[CsubsysA_bar[k]] =
								midx_rho_row[CsubsysA_bar[k]] =
										midxA_bar_row[k];
					// compute the results's row index
					std::size_t result_row_idx = internal::_multiidx2n(midx_row,
							n, Cdims);

					// compute the coefficient
					typename Derived1::Scalar coeff = 0;
					for (std::size_t c = 0; c < DA; c++)
					{
						// compute the subsys part state's row multi-index
						internal::_n2multiidx(c, nA, CdimsA, midxA_rho_row);
						// now we have the total state's row multi-index
						for (std::size_t k = 0; k < nA; k++)
							midx_rho_row[CsubsysA[k]] = midxA_rho_row[k];

						coeff += rA(a, c)
								* rstate(
										internal::_multiidx2n(midx_rho_row, n,
												Cdims));
					}
					// write down the result
					result(result_row_idx) = coeff;
				}
			}
			return result;
		}
		else if (internal::_check_square_mat(rstate)) // we have a matrix
		{

			// check that dims match state matrix
			if (!internal::_check_dims_match_mat(dims, rstate))
				throw Exception("Gates::apply",
						Exception::Type::DIMS_MISMATCH_MATRIX);

			// check that state matrix matches the dimensions of the subsys
			if (static_cast<std::size_t>(rA.cols()) != DA)
				throw Exception("Gates::apply",
						Exception::Type::DIMS_MISMATCH_MATRIX);

			types::DynMat<typename Derived1::Scalar> result(D, D);

			// run over the subsys's row multi-index
			for (std::size_t a = 0; a < DA; a++)
			{
				// get the subsys's row multi-index
				internal::_n2multiidx(a, nA, CdimsA, midxA_row);
				// compute subsys part of the result's row multi-index
				for (std::size_t k = 0; k < nA; k++)
					midx_row[CsubsysA[k]] = midxA_row[k];

				// run over the complement's row multi-index
				for (std::size_t i = 0; i < DA_bar; i++)
				{
					// get the complement's row multi-index
					internal::_n2multiidx(i, nA_bar, CdimsA_bar, midxA_bar_row);
					// now compute the complement part
					// of the result's row multi-index
					// and the complement part of the
					// state's total row multi-index
					for (std::size_t k = 0; k < nA_bar; k++)
						midx_row[CsubsysA_bar[k]] =
								midx_rho_row[CsubsysA_bar[k]] =
										midxA_bar_row[k];
					// compute the results's row index
					std::size_t result_row_idx = internal::_multiidx2n(midx_row,
							n, Cdims);

					// run over the col index
					for (std::size_t j = 0; j < D; j++)
					{
						// compute the coefficient
						typename Derived1::Scalar coeff = 0;
						for (std::size_t c = 0; c < DA; c++)
						{
							// compute the subsys part state's row multi-index
							internal::_n2multiidx(c, nA, CdimsA, midxA_rho_row);
							// now we have the total state's row multi-index
							for (std::size_t k = 0; k < nA; k++)
								midx_rho_row[CsubsysA[k]] = midxA_rho_row[k];

							coeff += rA(a, c)
									* rstate(
											internal::_multiidx2n(midx_rho_row,
													n, Cdims), j);

						}
						// write down the result
						result(result_row_idx, j) = coeff;
					}
				}
			}
			return result;
		}
		else
			throw Exception("Gates::apply",
					Exception::Type::MATRIX_NOT_SQUARE_OR_CVECTOR);
	}

	// returns a multi-quDit multi-controlled-gate in matrix form
	template<typename Derived>
	types::DynMat<typename Derived::Scalar> CTRL(
			const Eigen::MatrixBase<Derived>& A,
			const std::vector<std::size_t>& ctrl,
			const std::vector<std::size_t>& subsys, std::size_t n,
			std::size_t d = 2) const
	{
		const types::DynMat<typename Derived::Scalar> & rA = A;

		// EXCEPTION CHECKS
		// check matrix zero size
		if (!internal::_check_nonzero_size(rA))
			throw Exception("Gates::CTRL", Exception::Type::ZERO_SIZE);

		// check square matrix
		if (!internal::_check_square_mat(rA))
			throw Exception("Gates::CTRL", Exception::Type::MATRIX_NOT_SQUARE);

		// check lists zero size
		if (ctrl.size() == 0)
			throw Exception("Gates::CTRL", Exception::Type::ZERO_SIZE);
		if (subsys.size() == 0)
			throw Exception("Gates::CTRL", Exception::Type::ZERO_SIZE);

		// check out of range
		if (n == 0)
			throw Exception("Gates::CTRL", Exception::Type::OUT_OF_RANGE);

		// check valid local dimension
		if (d == 0)
			throw Exception("Gates::CTRL", Exception::Type::DIMS_INVALID);

		std::vector<std::size_t> ctrlgate = ctrl; // ctrl + gate subsystem vector
		ctrlgate.insert(std::end(ctrlgate), std::begin(subsys),
				std::end(subsys));

		std::vector<std::size_t> dims; // local dimensions vector
		dims.insert(std::begin(dims), n, d);

		// check that ctrl + gate subsystem is valid
		// with respect to local dimensions
		if (!internal::_check_subsys_match_dims(ctrlgate, dims))
			throw Exception("Gates::CTRL",
					Exception::Type::SUBSYS_MISMATCH_DIMS);

		// check that subsys list match the dimension of the matrix
		if (rA.cols() != std::pow(d, subsys.size()))
			throw Exception("Gates::CTRL",
					Exception::Type::DIMS_MISMATCH_MATRIX);
		// END EXCEPTION CHECKS

		// Use static allocation for speed!
		std::size_t Cdims[ct::maxn];
		std::size_t midx_row[ct::maxn];
		std::size_t midx_col[ct::maxn];

		std::size_t CdimsA[ct::maxn];
		std::size_t midxA_row[ct::maxn];
		std::size_t midxA_col[ct::maxn];

		std::size_t Cdims_bar[ct::maxn];
		std::size_t Csubsys_bar[ct::maxn];
		std::size_t midx_bar[ct::maxn];

		std::size_t ngate = subsys.size();
		std::size_t nctrl = ctrl.size();
		std::size_t nsubsys_bar = n - ctrlgate.size();
		std::size_t D = std::pow(d, n);
		std::size_t DA = static_cast<std::size_t>(rA.cols());
		std::size_t Dsubsys_bar = std::pow(d, nsubsys_bar);

		for (std::size_t k = 0, cnt = 0; k < n; k++)
		{
			midx_row[k] = midx_col[k] = 0;
			Cdims[k] = d;

			// compute the complementary subsystem of ctrlgate w.r.t. dims
			if (std::find(std::begin(ctrlgate), std::end(ctrlgate), k)
					== std::end(ctrlgate))
			{
				Csubsys_bar[cnt] = k;
				Cdims_bar[cnt] = d;
				midx_bar[cnt] = 0;
				cnt++;
			}
		}

		for (std::size_t k = 0; k < ngate; k++)
		{
			midxA_row[k] = midxA_col[k] = 0;
			CdimsA[k] = d;
		}

		types::DynMat<typename Derived::Scalar> result = types::DynMat<
				typename Derived::Scalar>::Identity(D, D);
		types::DynMat<typename Derived::Scalar> Ak;

		// run over the complement indexes
		for (std::size_t i = 0; i < Dsubsys_bar; i++)
		{
			// get the complement's row multi-index
			internal::_n2multiidx(i, nsubsys_bar, Cdims_bar, midx_bar);
			for (std::size_t k = 0; k < d; k++)
			{
				Ak = powm(rA, k); // compute rA^k
				// run over the subsys's row multi-index
				for (std::size_t a = 0; a < DA; a++)
				{
					// get the subsys's row multi-index
					internal::_n2multiidx(a, ngate, CdimsA, midxA_row);

					// construct the result's row multi-index

					// first the ctrl part (equal for both row and column)
					for (std::size_t c = 0; c < nctrl; c++)
						midx_row[ctrl[c]] = midx_col[ctrl[c]] = k;

					// then the complement part (equal for column)
					for (std::size_t c = 0; c < nsubsys_bar; c++)
						midx_row[Csubsys_bar[c]] = midx_col[Csubsys_bar[c]] =
								midx_bar[c];

					// then the subsys part
					for (std::size_t c = 0; c < ngate; c++)
						midx_row[subsys[c]] = midxA_row[c];

					// run over the subsys's column multi-index
					for (std::size_t b = 0; b < DA; b++)
					{
						// get the subsys's column multi-index
						internal::_n2multiidx(b, ngate, CdimsA, midxA_col);

						// construct the result's column multi-index
						for (std::size_t c = 0; c < ngate; c++)
							midx_col[subsys[c]] = midxA_col[c];

						// finally write the values
						result(internal::_multiidx2n(midx_row, n, Cdims),
								internal::_multiidx2n(midx_col, n, Cdims)) = Ak(
								a, b);
					}
				}

			}
		}
		return result;
	}

};
/* class Gates */

} /* namespace qpp */

#endif /* GATES_H_ */
