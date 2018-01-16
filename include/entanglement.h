/*
 * This file is part of Quantum++.
 *
 * MIT License
 *
 * Copyright (c) 2013 - 2018 Vlad Gheorghiu (vgheorgh@gmail.com)
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
* \file entanglement.h
* \brief Entanglement functions
*/

#ifndef ENTANGLEMENT_H_
#define ENTANGLEMENT_H_

namespace qpp {
/**
* \brief Schmidt coefficients of the bi-partite pure state \a A
*
* \note The sum of the squares of the Schmidt coefficients equals 1
* \see qpp::schmidtprobs()
*
* \param A Eigen expression
* \param dims Dimensions of the bi-partite system
* \return Schmidt coefficients of \a A, ordered in decreasing order, as a
* real dynamic column vector
*/
template <typename Derived>
dyn_col_vect<double> schmidtcoeffs(const Eigen::MatrixBase<Derived>& A,
                                   const std::vector<idx>& dims) {
    const dyn_mat<typename Derived::Scalar>& rA = A.derived();

    // EXCEPTION CHECKS

    // check zero-size
    if (!internal::check_nonzero_size(rA))
        throw exception::ZeroSize("qpp::schmidtcoeffs()");
    // check bi-partite
    if (dims.size() != 2)
        throw exception::NotBipartite("qpp::schmidtcoeffs()");
    // check column vector
    if (!internal::check_cvector(rA))
        throw exception::MatrixNotCvector("qpp::schmidtcoeffs()");
    // check matching dimensions
    if (!internal::check_dims_match_cvect(dims, rA))
        throw exception::DimsMismatchCvector("qpp::schmidtcoeffs()");
    // END EXCEPTION CHECKS

    return svals(transpose(reshape(rA, dims[1], dims[0])));
}

/**
* \brief Schmidt coefficients of the bi-partite pure state \a A
*
* \note The sum of the squares of the Schmidt coefficients equals 1
* \see qpp::schmidtprobs()
*
* \param A Eigen expression
* \param d Subsystem dimensions
* \return Schmidt coefficients of \a A, ordered in decreasing order, as a
* real dynamic column vector
*/
template <typename Derived>
dyn_col_vect<double> schmidtcoeffs(const Eigen::MatrixBase<Derived>& A,
                                   idx d = 2) {
    const dyn_mat<typename Derived::Scalar>& rA = A.derived();

    // EXCEPTION CHECKS

    // check zero size
    if (!internal::check_nonzero_size(A))
        throw exception::ZeroSize("qpp::schmidtcoeffs()");

    // check valid dims
    if (d < 2)
        throw exception::DimsInvalid("qpp::schmidtcoeffs()");
    // END EXCEPTION CHECKS

    idx N = internal::get_num_subsys(static_cast<idx>(rA.rows()), d);
    std::vector<idx> dims(N, d); // local dimensions vector

    return schmidtcoeffs(A, dims);
}

/**
* \brief Schmidt basis on Alice side
*
* \param A Eigen expression
* \param dims Dimensions of the bi-partite system
* \return Unitary matrix \f$ U \f$ whose columns represent
* the Schmidt basis vectors on Alice side.
*/
template <typename Derived>
cmat schmidtA(const Eigen::MatrixBase<Derived>& A,
              const std::vector<idx>& dims) {
    const dyn_mat<typename Derived::Scalar>& rA = A.derived();

    // EXCEPTION CHECKS

    // check zero-size
    if (!internal::check_nonzero_size(rA))
        throw exception::ZeroSize("qpp::schmidtA()");
    // check bi-partite
    if (dims.size() != 2)
        throw exception::NotBipartite("qpp::schmidtA()");
    // check column vector
    if (!internal::check_cvector(rA))
        throw exception::MatrixNotCvector("qpp::schmidtA()");
    // check matching dimensions
    if (!internal::check_dims_match_cvect(dims, rA))
        throw exception::DimsMismatchCvector("qpp::schmidtA()");
    // END EXCEPTION CHECKS

    return svdU(transpose(reshape(rA, dims[1], dims[0])));
}

/**
* \brief Schmidt basis on Alice side
*
* \param A Eigen expression
* \param d Subsystem dimensions
* \return Unitary matrix \f$ U \f$ whose columns represent
* the Schmidt basis vectors on Alice side.
*/
template <typename Derived>
cmat schmidtA(const Eigen::MatrixBase<Derived>& A, idx d = 2) {
    const dyn_mat<typename Derived::Scalar>& rA = A.derived();

    // EXCEPTION CHECKS

    // check zero size
    if (!internal::check_nonzero_size(A))
        throw exception::ZeroSize("qpp::schmidtA()");

    // check valid dims
    if (d < 2)
        throw exception::DimsInvalid("qpp::schmidtA()");
    // END EXCEPTION CHECKS

    idx N = internal::get_num_subsys(static_cast<idx>(rA.rows()), d);
    std::vector<idx> dims(N, d); // local dimensions vector

    return schmidtA(A, dims);
}

/**
* \brief Schmidt basis on Bob side
*
* \param A Eigen expression
* \param dims Dimensions of the bi-partite system
* \return Unitary matrix \f$ V \f$ whose columns represent
* the Schmidt basis vectors on Bob side.
*/
template <typename Derived>
cmat schmidtB(const Eigen::MatrixBase<Derived>& A,
              const std::vector<idx>& dims) {
    const dyn_mat<typename Derived::Scalar>& rA = A.derived();

    // EXCEPTION CHECKS

    // check zero-size
    if (!internal::check_nonzero_size(rA))
        throw exception::ZeroSize("qpp::schmidtB()");
    // check bi-partite
    if (dims.size() != 2)
        throw exception::NotBipartite("qpp::schmidtB()");
    // check column vector
    if (!internal::check_cvector(rA))
        throw exception::MatrixNotCvector("qpp::schmidtB()");
    // check matching dimensions
    if (!internal::check_dims_match_cvect(dims, rA))
        throw exception::DimsMismatchCvector("qpp::schmidtB()");
    // END EXCEPTION CHECKS

    // by default returns U_B^*, we need U_B, i.e. the complex conjugate,
    // i.e. adjoint(transpose(U_B))
    return svdV(transpose(reshape(conjugate(rA), dims[1], dims[0])));
}

/**
* \brief Schmidt basis on Bob side
*
* \param A Eigen expression
* \param d Subsystem dimensions
* \return Unitary matrix \f$ V \f$ whose columns represent
* the Schmidt basis vectors on Bob side.
*/
template <typename Derived>
cmat schmidtB(const Eigen::MatrixBase<Derived>& A, idx d = 2) {
    const dyn_mat<typename Derived::Scalar>& rA = A.derived();

    // EXCEPTION CHECKS

    // check zero size
    if (!internal::check_nonzero_size(A))
        throw exception::ZeroSize("qpp::schmidtB()");

    // check valid dims
    if (d < 2)
        throw exception::DimsInvalid("qpp::schmidtB()");
    // END EXCEPTION CHECKS

    idx N = internal::get_num_subsys(static_cast<idx>(rA.rows()), d);
    std::vector<idx> dims(N, d); // local dimensions vector

    return schmidtB(A, dims);
}

/**
* \brief Schmidt probabilities of the bi-partite pure state \a A
*
* Defined as the squares of the Schmidt coefficients.
* The sum of the Schmidt probabilities equals 1.
* \see qpp::schmidtcoeffs()
*
* \param A Eigen expression
* \param dims Dimensions of the bi-partite system
* \return Real vector consisting of the Schmidt probabilites of \a A,
* ordered in decreasing order
*/
template <typename Derived>
std::vector<double> schmidtprobs(const Eigen::MatrixBase<Derived>& A,
                                 const std::vector<idx>& dims) {
    const dyn_mat<typename Derived::Scalar>& rA = A.derived();

    // EXCEPTION CHECKS

    // check zero-size
    if (!internal::check_nonzero_size(rA))
        throw exception::ZeroSize("qpp::schmidtprobs()");
    // check bi-partite
    if (dims.size() != 2)
        throw exception::NotBipartite("qpp::schmidtprobs()");
    // check column vector
    if (!internal::check_cvector(rA))
        throw exception::MatrixNotCvector("qpp::schmidtprobs()");
    // check matching dimensions
    if (!internal::check_dims_match_cvect(dims, rA))
        throw exception::DimsMismatchCvector("qpp::schmidtprobs()");
    // END EXCEPTION CHECKS

    std::vector<double> result;
    dyn_col_vect<double> scf = schmidtcoeffs(rA, dims);
    for (idx i = 0; i < static_cast<idx>(scf.rows()); ++i)
        result.push_back(std::pow(scf(i), 2));

    return result;
}

/**
* \brief Schmidt probabilities of the bi-partite pure state \a A
*
* Defined as the squares of the Schmidt coefficients.
* The sum of the Schmidt probabilities equals 1.
* \see qpp::schmidtcoeffs()
*
* \param A Eigen expression
* \param d Subsystem dimensions
* \return Real vector consisting of the Schmidt probabilites of \a A,
* ordered in decreasing order
*/
template <typename Derived>
std::vector<double> schmidtprobs(const Eigen::MatrixBase<Derived>& A,
                                 idx d = 2) {
    const dyn_mat<typename Derived::Scalar>& rA = A.derived();

    // EXCEPTION CHECKS

    // check zero size
    if (!internal::check_nonzero_size(A))
        throw exception::ZeroSize("qpp::schmidtprobs()");

    // check valid dims
    if (d < 2)
        throw exception::DimsInvalid("qpp::schmidtprobs()");
    // END EXCEPTION CHECKS

    idx N = internal::get_num_subsys(static_cast<idx>(rA.rows()), d);
    std::vector<idx> dims(N, d); // local dimensions vector

    return schmidtprobs(A, dims);
}

/**
* \brief Entanglement of the bi-partite pure state \a A
*
* Defined as the von-Neumann entropy of the reduced density matrix
* of one of the subsystems
* \see qpp::entropy()
*
* \param A Eigen expression
* \param dims Dimensions of the bi-partite system
* \return Entanglement, with the logarithm in base 2
*/
template <typename Derived>
double entanglement(const Eigen::MatrixBase<Derived>& A,
                    const std::vector<idx>& dims) {
    const dyn_mat<typename Derived::Scalar>& rA = A.derived();

    // EXCEPTION CHECKS

    // check zero-size
    if (!internal::check_nonzero_size(rA))
        throw exception::ZeroSize("qpp::entanglement()");
    // check bi-partite
    if (dims.size() != 2)
        throw exception::NotBipartite("qpp::entanglement()");
    // check column vector
    if (!internal::check_cvector(rA))
        throw exception::MatrixNotCvector("qpp::entanglement()");
    // check matching dimensions
    if (!internal::check_dims_match_cvect(dims, rA))
        throw exception::DimsMismatchCvector("qpp::entanglement()");
    // END EXCEPTION CHECKS

    return entropy(schmidtprobs(rA, dims));
}

/**
* \brief Entanglement of the bi-partite pure state \a A
*
* Defined as the von-Neumann entropy of the reduced density matrix
* of one of the subsystems
* \see qpp::entropy()
*
* \param A Eigen expression
* \param d Subsystem dimensions
* \return Entanglement, with the logarithm in base 2
*/
template <typename Derived>
double entanglement(const Eigen::MatrixBase<Derived>& A, idx d = 2) {
    const dyn_mat<typename Derived::Scalar>& rA = A.derived();

    // EXCEPTION CHECKS

    // check zero size
    if (!internal::check_nonzero_size(A))
        throw exception::ZeroSize("qpp::entanglement()");

    // check valid dims
    if (d < 2)
        throw exception::DimsInvalid("qpp::entanglement()");
    // END EXCEPTION CHECKS

    idx N = internal::get_num_subsys(static_cast<idx>(rA.rows()), d);
    std::vector<idx> dims(N, d); // local dimensions vector

    return entanglement(A, dims);
}

/**
* \brief G-concurrence of the bi-partite pure state \a A
*
* \note Both local dimensions must be equal
*
* Uses qpp::logdet() to avoid overflows
* \see qpp::logdet()
*
* \param A Eigen expression
* \return G-concurrence
*/
// the G-concurrence
template <typename Derived>
double gconcurrence(const Eigen::MatrixBase<Derived>& A) {
    const dyn_mat<typename Derived::Scalar>& rA = A.derived();

    // EXCEPTION CHECKS

    // check zero-size
    if (!internal::check_nonzero_size(rA))
        throw exception::ZeroSize("qpp::gconcurrence()");
    // check column vector
    if (!internal::check_cvector(rA))
        throw exception::MatrixNotCvector("qpp::gconcurrence()");

    idx d = internal::get_dim_subsys(static_cast<idx>(rA.rows()), 2);

    // check equal local dimensions
    if (d * d != static_cast<idx>(rA.rows()))
        throw exception::DimsNotEqual("qpp::gconcurrence()");
    // END EXCEPTION CHECKS

    // we compute exp(logdet()) to avoid underflow
    return d * std::abs(std::exp(2. / d * logdet(reshape(rA, d, d))));
}

/**
* \brief Negativity of the bi-partite mixed state \a A
*
* \param A Eigen expression
* \param dims Dimensions of the bi-partite system
* \return Negativity
*/
template <typename Derived>
double negativity(const Eigen::MatrixBase<Derived>& A,
                  const std::vector<idx>& dims) {
    const dyn_mat<typename Derived::Scalar>& rA = A.derived();

    // EXCEPTION CHECKS

    // check zero-size
    if (!internal::check_nonzero_size(rA))
        throw exception::ZeroSize("qpp::negativity()");
    // check bi-partite
    if (dims.size() != 2)
        throw exception::NotBipartite("qpp::negativity()");
    // check square matrix vector
    if (!internal::check_square_mat(rA))
        throw exception::MatrixNotSquare("qpp::negativity()");
    // check matching dimensions
    if (!internal::check_dims_match_mat(dims, rA))
        throw exception::DimsMismatchMatrix("qpp::negativity()");
    // END EXCEPTION CHECKS

    return (schatten(ptranspose(rA, {0}, dims), 1) - 1.) / 2.;
}

/**
* \brief Negativity of the bi-partite mixed state \a A
*
* \param A Eigen expression
* \param d Subsystem dimensions
* \return Negativity
*/
template <typename Derived>
double negativity(const Eigen::MatrixBase<Derived>& A, idx d = 2) {
    const dyn_mat<typename Derived::Scalar>& rA = A.derived();

    // EXCEPTION CHECKS

    // check zero size
    if (!internal::check_nonzero_size(A))
        throw exception::ZeroSize("qpp::negativity()");

    // check valid dims
    if (d < 2)
        throw exception::DimsInvalid("qpp::negativity()");
    // END EXCEPTION CHECKS

    idx N = internal::get_num_subsys(static_cast<idx>(rA.rows()), d);
    std::vector<idx> dims(N, d); // local dimensions vector

    return negativity(A, dims);
}

/**
* \brief Logarithmic negativity of the bi-partite mixed state \a A
*
* \param A Eigen expression
* \param dims Dimensions of the bi-partite system
* \return Logarithmic negativity, with the logarithm in base 2
*/
template <typename Derived>
double lognegativity(const Eigen::MatrixBase<Derived>& A,
                     const std::vector<idx>& dims) {
    const dyn_mat<typename Derived::Scalar>& rA = A.derived();

    // EXCEPTION CHECKS

    // check zero-size
    if (!internal::check_nonzero_size(rA))
        throw exception::ZeroSize("qpp::lognegativity()");
    // check bi-partite
    if (dims.size() != 2)
        throw exception::NotBipartite("qpp::lognegativity()");
    // check square matrix vector
    if (!internal::check_square_mat(rA))
        throw exception::MatrixNotSquare("qpp::lognegativity()");
    // check matching dimensions
    if (!internal::check_dims_match_mat(dims, rA))
        throw exception::DimsMismatchMatrix("qpp::lognegativity()");
    // END EXCEPTION CHECKS

    return std::log2(2 * negativity(rA, dims) + 1);
}

/**
* \brief Logarithmic negativity of the bi-partite mixed state \a A
*
* \param A Eigen expression
* \param d Subsystem dimensions
* \return Logarithmic negativity, with the logarithm in base 2
*/
template <typename Derived>
double lognegativity(const Eigen::MatrixBase<Derived>& A, idx d = 2) {
    const dyn_mat<typename Derived::Scalar>& rA = A.derived();

    // EXCEPTION CHECKS

    // check zero size
    if (!internal::check_nonzero_size(A))
        throw exception::ZeroSize("qpp::lognegativity()");

    // check valid dims
    if (d < 2)
        throw exception::DimsInvalid("qpp::lognegativity()");
    // END EXCEPTION CHECKS

    idx N = internal::get_num_subsys(static_cast<idx>(rA.rows()), d);
    std::vector<idx> dims(N, d); // local dimensions vector

    return lognegativity(A, dims);
}

/**
* \brief Wootters concurrence of the bi-partite qubit mixed state \a A
*
* \param A Eigen expression
* \return Wootters concurrence
*/
template <typename Derived>
double concurrence(const Eigen::MatrixBase<Derived>& A) {
    const dyn_mat<typename Derived::Scalar>& rA = A.derived();

    // EXCEPTION CHECKS

    // check zero-size
    if (!internal::check_nonzero_size(rA))
        throw exception::ZeroSize("qpp::concurrence()");
    // check square matrix vector
    if (!internal::check_square_mat(rA))
        throw exception::MatrixNotSquare("qpp::concurrence()");
    // check that the state is a 2-qubit state
    if (rA.rows() != 4)
        throw exception::NotQubitSubsys("qpp::concurrence()");
    // END EXCEPTION CHECKS

    cmat sigmaY = Gates::get_instance().Y;
    dyn_col_vect<double> lambdas =
        evals(rA * kron(sigmaY, sigmaY) * conjugate(rA) * kron(sigmaY, sigmaY))
            .real();

    std::vector<double> lambdas_sorted(lambdas.data(),
                                       lambdas.data() + lambdas.size());

    std::sort(std::begin(lambdas_sorted), std::end(lambdas_sorted),
              std::greater<double>());
    std::transform(std::begin(lambdas_sorted), std::end(lambdas_sorted),
                   std::begin(lambdas_sorted), [](double elem) {
                       return std::sqrt(std::abs(elem));
                   }); // chop tiny negatives

    return std::max(0., lambdas_sorted[0] - lambdas_sorted[1] -
                            lambdas_sorted[2] - lambdas_sorted[3]);
}

} /* namespace qpp */

#endif /* ENTANGLEMENT_H_ */
