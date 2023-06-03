/*
 * This file is part of Quantum++.
 *
 * Copyright (c) 2013 - 2023 softwareQ Inc. All rights reserved.
 *
 * MIT License
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

/**
 * \file classes/gates.hpp
 * \brief Quantum gates
 */

#ifndef QPP_CLASSES_GATES_HPP_
#define QPP_CLASSES_GATES_HPP_

namespace qpp {
/**
 * \class qpp::Gates
 * \brief const Singleton class that implements most commonly used gates
 */
class Gates final : public internal::Singleton<const Gates> // const Singleton
{
    friend class internal::Singleton<const Gates>;

  public:
    // One qubit gates
    cmat Id2{cmat::Identity(2, 2)}; ///< Identity gate
    cmat H{cmat::Zero(2, 2)};       ///< Hadamard gate
    cmat X{cmat::Zero(2, 2)};       ///< Pauli Sigma-X gate
    cmat Y{cmat::Zero(2, 2)};       ///< Pauli Sigma-Y gate
    cmat Z{cmat::Zero(2, 2)};       ///< Pauli Sigma-Z gate
    cmat S{cmat::Zero(2, 2)};       ///< S gate
    cmat T{cmat::Zero(2, 2)};       ///< T gate

    // two qubit gates
    cmat CNOT{cmat::Identity(4, 4)}; ///< Controlled-NOT control target gate
    cmat CZ{cmat::Identity(4, 4)};   ///< Controlled-Phase gate
    cmat CNOTba{cmat::Zero(4, 4)};   ///< Controlled-NOT target->control gate
    cmat SWAP{cmat::Identity(4, 4)}; ///< SWAP gate

    // three qubit gates
    cmat TOF{cmat::Identity(8, 8)};  ///< Toffoli gate
    cmat FRED{cmat::Identity(8, 8)}; ///< Fredkin gate
  private:
    /**
     * \brief Initializes the gates
     */
    Gates() {
        H << 1 / std::sqrt(2.), 1 / std::sqrt(2.), 1 / std::sqrt(2.),
            -1 / std::sqrt(2.);
        X << 0, 1, 1, 0;
        Z << 1, 0, 0, -1;
        Y << 0, -1_i, 1_i, 0;
        S << 1, 0, 0, 1_i;
        T << 1, 0, 0, std::exp(1_i * static_cast<cplx::value_type>(pi / 4.0));
        CNOT.block(2, 2, 2, 2) = X;
        CNOTba(0, 0) = 1;
        CNOTba(1, 3) = 1;
        CNOTba(2, 2) = 1;
        CNOTba(3, 1) = 1;
        CZ(3, 3) = -1;

        SWAP.block(1, 1, 2, 2) = X;
        TOF.block(6, 6, 2, 2) = X;
        FRED.block(4, 4, 4, 4) = SWAP;
    }

    /**
     * \brief Default destructor
     */
    ~Gates() override = default;

  public:
    // variable gates

    // one qubit gates

    /**
     * \brief Qubit rotation of \a theta about the 3-dimensional real (unit)
     * vector \a n
     *
     * \param theta Rotation angle
     * \param n 3-dimensional real (unit) vector
     * \return Rotation gate
     */
    cmat Rn(realT theta, const std::array<realT, 3>& n) const {
        cmat result(2, 2);
        result = std::cos(theta / 2) * Id2 -
                 1_i * static_cast<cplx::value_type>(std::sin(theta / 2)) *
                     (n[0] * X + n[1] * Y + n[2] * Z);

        return result;
    }

    /**
     * \brief Qubit rotation of \a theta about the X axis
     *
     * \param theta Rotation angle
     * \return Rotation gate
     */
    cmat RX(realT theta) const {
        // EXCEPTION CHECKS

        // END EXCEPTION CHECKS

        return Rn(theta, {1, 0, 0});
    }

    /**
     * \brief Qubit rotation of \a theta about the Y axis
     *
     * \param theta Rotation angle
     * \return Rotation gate
     */
    cmat RY(realT theta) const {
        // EXCEPTION CHECKS

        // END EXCEPTION CHECKS

        return Rn(theta, {0, 1, 0});
    }

    /**
     * \brief Qubit rotation of \a theta about the Z axis
     *
     * \param theta Rotation angle
     * \return Rotation gate
     */
    cmat RZ(realT theta) const {
        // EXCEPTION CHECKS

        // END EXCEPTION CHECKS

        return Rn(theta, {0, 0, 1});
    }

    // one quDit gates

    /**
     * \brief Generalized Z gate for qudits
     *
     * \note Defined as \f$Z = \sum_{j=0}^{D-1} \exp(2\pi \mathrm{i} j/D)
     * |j\rangle\langle j|\f$
     *
     * \param D Dimension of the Hilbert space
     * \return Generalized Z gate for qudits
     */
    cmat Zd(idx D = 2) const {
        // EXCEPTION CHECKS

        // check valid dimension
        if (D == 0)
            throw exception::DimsInvalid("qpp::Gates::Zd()", "D");
        // END EXCEPTION CHECKS

        if (D == 2)
            return Z;

        cmat result = cmat::Zero(D, D);
        for (idx i = 0; i < D; ++i)
            result(i, i) = std::pow(omega(D), static_cast<realT>(i));

        return result;
    }

    /**
     * \brief SWAP gate for qudits
     *
     * \param D Dimension of the Hilbert space
     * \return SWAP gate for qudits
     */
    [[qpp::parallel]] cmat SWAPd(idx D = 2) const {
        // EXCEPTION CHECKS

        // check valid dimension
        if (D == 0)
            throw exception::DimsInvalid("qpp::Gates::SWAPd()", "D");
        // END EXCEPTION CHECKS

        if (D == 2)
            return SWAP;

        cmat result = cmat::Zero(D * D, D * D);

#ifdef QPP_OPENMP
// NOLINTNEXTLINE
#pragma omp parallel for collapse(2)
#endif // QPP_OPENMP
       // column major order for speed
        for (idx j = 0; j < D; ++j)
            for (idx i = 0; i < D; ++i)
                result(D * i + j, i + D * j) = 1;

        return result;
    }

    /**
     * \brief Quantum Fourier transform gate for qudits
     *
     * \note Defined as
     * \f$F = \sum_{j,k=0}^{D-1} \exp(2\pi \mathrm{i} jk/D) |j\rangle\langle k|
     * \f$
     *
     * \param D Dimension of the Hilbert space
     * \return Fourier transform gate for qudits
     */
    [[qpp::parallel]] cmat Fd(idx D = 2) const {
        // EXCEPTION CHECKS

        // check valid dimension
        if (D == 0)
            throw exception::DimsInvalid("qpp::Gates::Fd()", "D");
        // END EXCEPTION CHECKS

        if (D == 2)
            return H;

        cmat result(D, D);

#ifdef QPP_OPENMP
// NOLINTNEXTLINE
#pragma omp parallel for collapse(2)
#endif // QPP_OPENMP
       // column major order for speed
        for (idx j = 0; j < D; ++j)
            for (idx i = 0; i < D; ++i)
                result(i, j) = 1 / static_cast<cplx::value_type>(std::sqrt(D)) *
                               std::pow(omega(D), static_cast<realT>(i * j));

        return result;
    }

    /**
     * \brief Modular multiplication gate for qubits
     * Implements \f$|x\rangle  \longrightarrow |ax \mathrm{ mod } N\rangle\f$
     *
     * \note For the gate to be unitary, \a a and \a N should be co-prime. The
     * function does not check co-primality in release versions!
     *
     * \note The number of qubits required to implement the gate should satisfy
     * \f$n \geq \lceil\log_2(N)\rceil\f$
     *
     * \param a Positive integer less than \a N
     * \param N Positive integer
     * \param n Number of qubits required for implementing the gate
     * \return Modular multiplication gate
     */
    [[qpp::parallel]] cmat MODMUL(idx a, idx N, idx n) const {
        // check co-primality (unitarity) only in DEBUG version
        assert(gcd(a, N) == 1);
        // EXCEPTION CHECKS

        // check valid arguments
        if (N < 3 || a >= N) {
            throw exception::OutOfRange("qpp::Gates::MODMUL()", "a/N");
        }

        // check enough qubits
        if (n < static_cast<idx>(std::ceil(std::log2(N)))) {
            throw exception::OutOfRange("qpp::Gates::MODMUL()", "n/N");
        }
        // END EXCEPTION CHECKS

        // minimum number of qubits required to implement the gate
        idx D = static_cast<idx>(std::llround(std::pow(2, n)));

        cmat result = cmat::Zero(D, D);

#ifdef QPP_OPENMP
// NOLINTNEXTLINE
#pragma omp parallel for collapse(2)
#endif // QPP_OPENMP
       // column major order for speed
        for (idx j = 0; j < N; ++j)
            for (idx i = 0; i < N; ++i)
                if (static_cast<idx>(modmul(static_cast<bigint>(j),
                                            static_cast<bigint>(a),
                                            static_cast<bigint>(N))) == i)
                    result(i, j) = 1;

#ifdef QPP_OPENMP
// NOLINTNEXTLINE
#pragma omp parallel for
#endif // QPP_OPENMP
       // complete the matrix
        for (idx i = N; i < D; ++i)
            result(i, i) = 1;

        return result;
    }

    /**
     * \brief Generalized X gate for qudits
     *
     * \note Defined as \f$X = \sum_{j=0}^{D-1} |j\oplus 1\rangle\langle j|\f$,
     * i.e., raising operator \f$X|j\rangle = |j\oplus 1\rangle\f$
     *
     * \param D Dimension of the Hilbert space
     * \return Generalized X gate for qudits
     */
    cmat Xd(idx D = 2) const {
        // EXCEPTION CHECKS

        // check valid dimension
        if (D == 0)
            throw exception::DimsInvalid("qpp::Gates::Xd()", "D");
        // END EXCEPTION CHECKS

        if (D == 2)
            return X;

        return Fd(D).inverse() * Zd(D) * Fd(D);
    }

    /**
     * \brief Identity gate
     *
     * \note Can change the return type from complex matrix (default) by
     * explicitly specifying the template parameter
     *
     * \param D Dimension of the Hilbert space
     * \return Identity gate on a Hilbert space of dimension \a D
     */
    template <typename Derived = cmat>
    Derived Id(idx D = 2) const {
        // EXCEPTION CHECKS

        // check valid dimension
        if (D == 0)
            throw exception::DimsInvalid("qpp::Gates::Id()", "D");
        // END EXCEPTION CHECKS

        return Derived::Identity(D, D);
    }

    /**
     * \brief Generates the multi-partite \a A gate in matrix form
     * \see qpp::apply()
     *
     * \note The dimension of the gate \a A must match the dimension of
     * \a target
     *
     * \param A Eigen expression
     * \param target Target subsystem indexes where the gate \a A is applied
     * \param dims Dimensions of the multi-partite system
     * \return \a A gate, as a matrix over the same scalar field as \a A
     */
    template <typename Derived>
    dyn_mat<typename Derived::Scalar> GATE(const Eigen::MatrixBase<Derived>& A,
                                           const std::vector<idx>& target,
                                           const std::vector<idx>& dims) const {
        const dyn_mat<typename Derived::Scalar>& rA = A.derived();

        // EXCEPTION CHECKS

        // check matrix zero-size
        if (!internal::check_nonzero_size(rA))
            throw exception::ZeroSize("qpp::Gates::GATE()", "A");

        // check square matrix
        if (!internal::check_square_mat(rA))
            throw exception::MatrixNotSquare("qpp::Gates::GATE()", "A");

        // check zero-size
        if (target.empty())
            throw exception::ZeroSize("qpp::Gates::GATE()", "target");

        // check that dims is a valid dimension vector
        if (!internal::check_dims(dims))
            throw exception::DimsInvalid("qpp::Gates::GATE()", "dims");

        // check that target is valid w.r.t. dims
        if (!internal::check_subsys_match_dims(target, dims))
            throw exception::SubsysMismatchDims("qpp::Gates::GATE()",
                                                "dims/target");

        // check that target list match the dimension of the matrix
        using Index = typename dyn_mat<typename Derived::Scalar>::Index;

        idx DA = 1;
        for (idx elem : target)
            DA *= dims[elem];

        if (rA.rows() != static_cast<Index>(DA))
            throw exception::MatrixMismatchSubsys("qpp::Gates::GATE()",
                                                  "A/dims/target");

        // END EXCEPTION CHECKS

        // Use static allocation for speed!
        idx Cdims[internal::maxn];
        idx midx_row[internal::maxn];
        idx midx_col[internal::maxn];

        idx CdimsA[internal::maxn];
        idx midxA_row[internal::maxn];
        idx midxA_col[internal::maxn];

        idx CdimsA_bar[internal::maxn];
        idx Csubsys_bar[internal::maxn];
        idx midx_bar[internal::maxn];

        idx n = dims.size();
        idx n_gate = target.size();
        idx n_subsys_bar = n - target.size();

        // compute the complementary subsystem of ctrlgate w.r.t. dims
        std::vector<idx> subsys_bar = complement(target, n);

        idx D = prod(dims);
        idx Dsubsys_bar = 1;
        for (idx elem : subsys_bar)
            Dsubsys_bar *= dims[elem];

        std::copy(subsys_bar.begin(), subsys_bar.end(),
                  std::begin(Csubsys_bar));

        for (idx k = 0; k < n; ++k) {
            midx_row[k] = midx_col[k] = 0;
            Cdims[k] = dims[k];
        }

        for (idx k = 0; k < n_subsys_bar; ++k) {
            CdimsA_bar[k] = dims[subsys_bar[k]];
            midx_bar[k] = 0;
        }

        for (idx k = 0; k < n_gate; ++k) {
            midxA_row[k] = midxA_col[k] = 0;
            CdimsA[k] = dims[target[k]];
        }

        dyn_mat<typename Derived::Scalar> result =
            dyn_mat<typename Derived::Scalar>::Identity(D, D);

        // run over the complement indexes
        for (idx i = 0; i < Dsubsys_bar; ++i) {
            // get the complement row multi-index
            internal::n2multiidx(i, n_subsys_bar, CdimsA_bar, midx_bar);

            // run over the target row multi-index
            for (idx a = 0; a < DA; ++a) {
                // get the target row multi-index
                internal::n2multiidx(a, n_gate, CdimsA, midxA_row);

                // construct the result row multi-index

                // first the target part
                for (idx k = 0; k < n_gate; ++k)
                    midx_row[target[k]] = midxA_row[k];

                // then the complement part (equal for column)
                for (idx k = 0; k < n_subsys_bar; ++k)
                    midx_row[Csubsys_bar[k]] = midx_col[Csubsys_bar[k]] =
                        midx_bar[k];

                // run over the target column multi-index
                for (idx b = 0; b < DA; ++b) {
                    // get the target column multi-index
                    internal::n2multiidx(b, n_gate, CdimsA, midxA_col);

                    // construct the result column multi-index
                    for (idx k = 0; k < n_gate; ++k)
                        midx_col[target[k]] = midxA_col[k];

                    // finally write the values
                    result(internal::multiidx2n(midx_row, n, Cdims),
                           internal::multiidx2n(midx_col, n, Cdims)) = rA(a, b);
                }
            }
        }

        return result;
    }

    /**
     * \brief Generates the multi-partite \a A gate in matrix form
     * \see qpp::apply()
     *
     * \note The dimension of the gate \a A must match the dimension of
     * \a target
     *
     * \param A Eigen expression
     * \param target Target subsystem indexes where the gate \a A is applied
     * \param n Total number of subsystems
     * \param d Subsystem dimensions
     * \return \a A gate, as a matrix over the same scalar field as \a A
     */
    template <typename Derived>
    dyn_mat<typename Derived::Scalar> GATE(const Eigen::MatrixBase<Derived>& A,
                                           const std::vector<idx>& target,
                                           idx n, idx d = 2) const {
        // EXCEPTION CHECKS

        // check valid local dimension
        if (d == 0)
            throw exception::DimsInvalid("qpp::Gates::GATE()", "d");
        // END EXCEPTION CHECKS

        return GATE(A, target, std::vector<idx>(n, d));
    }

    /**
     * \brief Generates the multi-partite multiple-controlled-\a A gate in
     * matrix form
     * \see qpp::applyCTRL()
     *
     * \note The dimension of the gate \a A must match the dimension of
     * \a target. All subsystems must have the same dimension.
     *
     * \param A Eigen expression
     * \param ctrl Control subsystem indexes
     * \param target Target subsystem indexes where the gate \a A is applied
     * \param n Total number of subsystems
     * \param d Subsystem dimensions
     * \param shift Optional, performs the control as if the \a ctrl qudit
     * states were \f$X\f$-incremented component-wise by \a shift. If present,
     * the size of \a shift must be the same as the size of \a ctrl.
     * \return CTRL-A gate, as a matrix over the same scalar field as \a A
     */
    template <typename Derived>
    dyn_mat<typename Derived::Scalar>
    CTRL(const Eigen::MatrixBase<Derived>& A, const std::vector<idx>& ctrl,
         const std::vector<idx>& target, idx n, idx d = 2,
         std::optional<std::vector<idx>> shift = std::nullopt) const {
        const dyn_mat<typename Derived::Scalar>& rA = A.derived();

        // EXCEPTION CHECKS

        // check matrix zero-size
        if (!internal::check_nonzero_size(rA))
            throw exception::ZeroSize("qpp::Gates::CTRL()", "A");

        // check square matrix
        if (!internal::check_square_mat(rA))
            throw exception::MatrixNotSquare("qpp::Gates::CTRL()", "A");

        // check lists zero-size
        if (ctrl.empty())
            throw exception::ZeroSize("qpp::Gates::CTRL()", "ctrl");
        if (target.empty())
            throw exception::ZeroSize("qpp::Gates::CTRL()", "target");

        // check out of range
        if (n == 0)
            throw exception::OutOfRange("qpp::Gates::CTRL()", "n");

        // check valid local dimension
        if (d == 0)
            throw exception::DimsInvalid("qpp::Gates::CTRL()", "d");

        // ctrl + gate subsystem vector
        std::vector<idx> ctrlgate = ctrl;
        std::sort(ctrlgate.begin(), ctrlgate.end());
        ctrlgate.insert(ctrlgate.end(), target.begin(), target.end());
        // FIXME if needed
        // std::sort(ctrlgate.begin(), ctrlgate.end());

        // check ctrl and target don't share common elements
        for (idx elem_ctrl : ctrl)
            for (idx elem_target : target)
                if (elem_ctrl == elem_target)
                    throw exception::OutOfRange("qpp::Gates::CTRL()",
                                                "ctrl/target");

        std::vector<idx> dims(n, d); // local dimensions vector
        // check that ctrl + gate subsystem is valid
        // with respect to local dimensions
        if (!internal::check_subsys_match_dims(ctrlgate, dims))
            throw exception::SubsysMismatchDims("qpp::Gates::CTRL()",
                                                "ctrl/dims");

        // check that target list match the dimension of the matrix
        idx DA = rA.rows();
        if (DA != static_cast<idx>(std::llround(std::pow(d, target.size()))))
            throw exception::MatrixMismatchSubsys("qpp::Gates::CTRL()",
                                                  "A/d/target");

        // check shift
        if (shift.has_value() && (shift.value().size() != ctrl.size()))
            throw exception::SizeMismatch("qpp::Gates::CTRL()", "ctrl/shift");
        if (shift.has_value())
            for (idx& elem : shift.value()) {
                if (elem >= d)
                    throw exception::OutOfRange("qpp::Gates::CTRL()", "shift");
                elem = d - elem;
            }
        // END EXCEPTION CHECKS

        if (!shift.has_value())
            shift = std::vector<idx>(ctrl.size(), 0);

        idx D = prod(dims);
        idx Dctrl = static_cast<idx>(std::llround(std::pow(d, ctrl.size())));
        idx ctrl_size = ctrl.size();

        dyn_mat<typename Derived::Scalar> result =
            dyn_mat<typename Derived::Scalar>::Zero(D, D);

        std::vector<idx> ctrl_bar = complement(ctrlgate, n);
        std::vector<idx> ctrlgate_bar = complement(ctrlgate, n);
        idx Dctrlgate_bar = 1;
        for (idx elem : ctrlgate_bar)
            Dctrlgate_bar *= dims[elem];

        dyn_mat<typename Derived::Scalar> Id_ctrlgate_bar =
            dyn_mat<typename Derived::Scalar>::Identity(Dctrlgate_bar,
                                                        Dctrlgate_bar);

        dyn_mat<typename Derived::Scalar> prj_bar =
            dyn_mat<typename Derived::Scalar>::Identity(Dctrl, Dctrl);
        for (idx k = 0; k < d; ++k) {
            // copy shift
            std::vector<idx> ctrl_shift = shift.value();
            // increment ctrl_shift by k (mod d)
            std::transform(ctrl_shift.begin(), ctrl_shift.end(),
                           ctrl_shift.begin(),
                           [k, d](idx elem) { return (elem + k) % d; });
            // compute the projector multi-index
            idx pos = multiidx2n(ctrl_shift, std::vector<idx>(ctrl_size, d));

            // compute the projection matrix
            dyn_mat<typename Derived::Scalar> prj_mat =
                dyn_mat<typename Derived::Scalar>::Zero(Dctrl, Dctrl);
            prj_mat(pos, pos) = 1;

            // subtract from identity on the ctrl side
            prj_bar -= prj_mat;

            // now compute [prj] x Ak and add to result
            dyn_mat<typename Derived::Scalar> Ak = powm(rA, k);
            dyn_mat<typename Derived::Scalar> gate = kron(prj_mat, Ak);

            result += GATE(gate, ctrlgate, n, d);
        }

        // finally, add projector's complement
        result += GATE(
            kron(prj_bar, dyn_mat<typename Derived::Scalar>::Identity(DA, DA)),
            ctrlgate, n, d);

        return result;
    }

    /**
     * \brief Expands out
     * \see qpp::kron()
     *
     * Expands out \a A as a matrix in a multi-partite system. Faster than using
     * qpp::kron(I, I, ..., I, A, I, ..., I).
     *
     * \param A Eigen expression
     * \param pos Position
     * \param dims Dimensions of the multi-partite system
     * \return Tensor product
     * \f$I\otimes\cdots\otimes I\otimes A \otimes I \otimes\cdots\otimes I\f$,
     * with \a A on position \a pos, as a dynamic matrix over the same scalar
     * field as \a A
     */
    template <typename Derived>
    dyn_mat<typename Derived::Scalar>
    expandout(const Eigen::MatrixBase<Derived>& A, idx pos,
              const std::vector<idx>& dims) const {
        const dyn_mat<typename Derived::Scalar>& rA = A.derived();

        // EXCEPTION CHECKS

        // check zero-size
        if (!internal::check_nonzero_size(rA))
            throw exception::ZeroSize("qpp::Gates::expandout()", "A");

        // check that dims is a valid dimension vector
        if (!internal::check_dims(dims))
            throw exception::DimsInvalid("qpp::Gates::expandout()", "dims");

        // check square matrix
        if (!internal::check_square_mat(rA))
            throw exception::MatrixNotSquare("qpp::Gates::expandout()", "A");

        // check that position is valid
        if (pos + 1 > static_cast<idx>(dims.size()))
            throw exception::OutOfRange("qpp::Gates::expandout()", "dims/pos");

        // check that dims[pos] match the dimension of A
        if (static_cast<idx>(rA.rows()) != dims[pos])
            throw exception::DimsMismatchMatrix("qpp::Gates::expandout()",
                                                "A/dims");
        // END EXCEPTION CHECKS

        return GATE(A, {pos}, dims);
    }

    /**
     * \brief Expands out
     * \see qpp::kron()
     *
     * Expands out \a A as a matrix in a multi-partite system. Faster than using
     * qpp::kron(I, I, ..., I, A, I, ..., I).
     *
     * \note The std::initializer_list overload exists because otherwise, in the
     * degenerate case when \a dims has only one element, the one element list
     * is implicitly converted to the element's underlying type, i.e., qpp::idx,
     * which has the net effect of picking the wrong (non-vector)
     * qpp::expandout() overload
     *
     * \param A Eigen expression
     * \param pos Position
     * \param dims Dimensions of the multi-partite system
     * \return Tensor product
     * \f$I\otimes\cdots\otimes I\otimes A \otimes I \otimes\cdots\otimes I\f$,
     * with \a A on position \a pos, as a dynamic matrix over the same scalar
     * field as \a A
     */
    template <typename Derived>
    dyn_mat<typename Derived::Scalar>
    expandout(const Eigen::MatrixBase<Derived>& A, idx pos,
              const std::initializer_list<idx>& dims) const {
        return expandout(A, pos, std::vector<idx>(dims));
    }

    /**
     * \brief Expands out
     * \see qpp::kron()
     *
     * Expands out \a A as a matrix in a multi-partite system. Faster than using
     * qpp::kron(I, I, ..., I, A, I, ..., I).
     *
     * \param A Eigen expression
     * \param pos Position
     * \param n Number of subsystems
     * \param d Subsystem dimensions
     * \return Tensor product
     * \f$I\otimes\cdots\otimes I\otimes A \otimes I \otimes\cdots\otimes I\f$,
     * with \a A on position \a pos, as a dynamic matrix over the same scalar
     * field as \a A
     */
    template <typename Derived>
    dyn_mat<typename Derived::Scalar>
    expandout(const Eigen::MatrixBase<Derived>& A, idx pos, idx n,
              idx d = 2) const {
        // EXCEPTION CHECKS

        // check zero size
        if (!internal::check_nonzero_size(A))
            throw exception::ZeroSize("qpp::Gates::expandout()", "A");

        // check valid dims
        if (d == 0)
            throw exception::DimsInvalid("qpp::Gates::expandout()", "d");
        // END EXCEPTION CHECKS

        return expandout(A, pos, std::vector<idx>(n, d));
    }

    // getters

    /**
     * \brief Get the name of the most common qubit gates
     *
     * \note Assumes that the gate \a U is represented by a square matrix. If
     * not, returns the empty string
     *
     * \param U Complex matrix representing the quantum gate
     * \return Optional name of the gate
     */
    std::optional<std::string> get_name(const cmat& U) const {
        // EXCEPTION CHECKS

        // check zero size
        if (!internal::check_nonzero_size(U))
            throw exception::ZeroSize("qpp::Gates::get_name()", "U");

        // check square matrix
        if (!internal::check_square_mat(U))
            return {};

        // END EXCEPTION CHECKS

        const idx D = static_cast<idx>(U.rows());

        switch (D) {
                // 1 qubit gates
            case 2:
                if (U == Id2)
                    return "Id2";
                else if (U == H)
                    return "H";
                else if (U == X)
                    return "X";
                else if (U == Y)
                    return "Y";
                else if (U == Z)
                    return "Z";
                else if (U == S)
                    return "S";
                else if (U == adjoint(S))
                    return "S+";
                else if (U == T)
                    return "T";
                else if (U == adjoint(T))
                    return "T+";
                else
                    return {};
                // 2 qubit gates
            case 4:
                if (U == CNOT)
                    return "CNOT";
                else if (U == CZ)
                    return "CZ";
                else if (U == CNOTba)
                    return "CNOTba";
                else if (U == SWAP)
                    return "SWAP";
                else
                    return {};
                // 3 qubit gates
            case 8:
                if (U == TOF)
                    return "TOF";
                else if (U == FRED)
                    return "FRED";
                else
                    return {};

            default:
                return {};
        }
    }
    // end getters
}; /* class Gates */

} /* namespace qpp */

#endif /* QPP_CLASSES_GATES_HPP_ */
