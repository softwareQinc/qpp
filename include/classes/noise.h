/*
 * This file is part of Quantum++.
 *
 * MIT License
 *
 * Copyright (c) 2013 - 2019 Vlad Gheorghiu (vgheorgh@gmail.com)
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
 * \file classes/noise.h
 * \brief Noise models
 */

#ifndef CLASSES_NOISE_H_
#define CLASSES_NOISE_H_

namespace qpp {
/**
 * \class qpp::NoiseType
 * \brief Contains template tags used to specify the noise type
 */
struct NoiseType {
    /**
     * \class qpp::NoiseType::StateDependent
     * \brief Template tag, used whenever the noise is state-dependent
     */
    struct StateDependent;

    /**
     * \class qpp::NoiseType::StateIndependent
     * \brief Template tag, used whenever the noise is state-independent
     */
    struct StateIndependent;
};
/**
 * \class qpp::NoiseBase
 * \brief Base class for all noise models, derive your particular noise model
 */
template <class T>
class NoiseBase {
  public:
    using noise_type = T;
    static_assert(
        std::is_same<NoiseType::StateDependent, noise_type>::value ||
            std::is_same<NoiseType::StateIndependent, noise_type>::value,
        "");

  protected:
    const std::vector<cmat> Ks_;        ///< Kraus operators
    mutable std::vector<double> probs_; ///< probabilities
    mutable idx d_{};                   ///< qudit dimension

    mutable idx i_{}; ///< index of the last occurring noise element
    mutable bool generated_{false}; ///< set to true after compute_state_() is
    ///< invoked, or if the noise is state-independent

    /**
     * \brief Compute probability outcomes for StateDependent noise type,
     * otherwise returns without performing any operation (no-op)
     *
     * \param state State vector or density matrix
     * \param target Qudit indexes where the noise is applied
     */
    void compute_probs_(const cmat& state,
                        const std::vector<idx>& target) const {
        if (!std::is_same<NoiseType::StateDependent, noise_type>::value)
            return; // no-op

        // minimal EXCEPTION CHECKS

        if (!internal::check_nonzero_size(state))
            throw exception::ZeroSize("qpp::Noise::compute_probs_()");
        // END EXCEPTION CHECKS

        cmat rho_i;
        idx n = internal::get_num_subsys(state.rows(), d_);

        for (idx i = 0; i < Ks_.size(); ++i) {
            rho_i = ptrace(state, complement(target, n), d_);
            probs_[i] = trace(Ks_[i] * rho_i * adjoint(Ks_[i])).real();
        }
    } /* compute_probs_() */

    /**
     * \brief Compute the resulting state after the noise was applied
     *
     * \param state State vector or density matrix
     * \param target Qudit indexes where the noise is applied
     * \return Resulting state after the noise was applied
     */
    cmat compute_state_(const cmat& state,
                        const std::vector<idx>& target) const {
        cmat result;
        idx D = static_cast<idx>(state.rows());

        //************ ket ************//
        if (internal::check_cvector(state)) {
            result.resize(D, 1);
        }
        //************ density matrix ************//
        else if (internal::check_square_mat(state)) {
            result.resize(D, D);
        }
        //************ Exception: not ket nor density matrix ************//
        else
            throw exception::MatrixNotSquareNorCvector(
                "qpp::Noise::compute_state_()");

        // now do the actual noise generation
        std::discrete_distribution<idx> dd{std::begin(probs_),
                                           std::end(probs_)};
        auto& gen =
#ifdef NO_THREAD_LOCAL_
            RandomDevices::get_instance().get_prng();
#else
            RandomDevices::get_thread_local_instance().get_prng();
#endif
        i_ = dd(gen);
        result = apply(state, Ks_[i_], target, d_);
        generated_ = true;

        return normalize(result);
    } /* compute_state_() */

  public:
    /**
     * \brief Constructs a noise instance for StateDependent noise type
     *
     * \note SFINAEd-out for StateIndependent noise
     *
     * \param A Eigen expression (state vector or density matrix)
     * \param Ks Vector of noise (Kraus) operators that specify the noise
     * \param d Subsystem dimension
     */
    template <typename U = noise_type>
    explicit NoiseBase(
        const std::vector<cmat>& Ks,
        typename std::enable_if<
            std::is_same<NoiseType::StateDependent, U>::value>::type* = nullptr)
        : Ks_{Ks}, probs_(Ks.size()) {
        // EXCEPTION CHECKS

        if (Ks.size() == 0)
            throw exception::ZeroSize("qpp::Noise::Noise()");
        if (!internal::check_nonzero_size(Ks[0]))
            throw exception::ZeroSize("qpp::Noise::Noise()");
        if (!internal::check_square_mat(Ks[0]))
            throw exception::MatrixNotSquare("qpp::Noise::Noise()");
        for (auto&& elem : Ks)
            if (elem.rows() != Ks[0].rows() || elem.cols() != Ks[0].rows())
                throw exception::DimsNotEqual("qpp::Noise::Noise()");
        // END EXCEPTION CHECKS

        d_ = Ks[0].rows(); // set the local dimension
    }

    /**
     * \brief Constructs a noise instance for StateIndependent noise type
     *
     * \note SFINAEd-out for StateDependent noise
     *
     * \param A Eigen expression (state vector or density matrix)
     * \param Ks Vector of noise (Kraus) operators that specify the noise
     * \param d Subsystem dimension
     */
    template <typename U = noise_type>
    explicit NoiseBase(
        const std::vector<cmat>& Ks, const std::vector<double>& probs,
        typename std::enable_if<std::is_same<NoiseType::StateIndependent,
                                             U>::value>::type* = nullptr)
        : Ks_{Ks}, probs_(probs) {
        // EXCEPTION CHECKS

        if (Ks.size() == 0)
            throw exception::ZeroSize("qpp::Noise::Noise()");
        if (Ks.size() != probs.size())
            throw exception::SizeMismatch("qpp::Noise::Noise");
        if (!internal::check_nonzero_size(Ks[0]))
            throw exception::ZeroSize("qpp::Noise::Noise()");
        if (!internal::check_square_mat(Ks[0]))
            throw exception::MatrixNotSquare("qpp::Noise::Noise()");
        for (auto&& elem : Ks)
            if (elem.rows() != Ks[0].rows() || elem.cols() != Ks[0].rows())
                throw exception::DimsNotEqual("qpp::Noise::Noise()");
        for (auto&& elem : probs)
            if (elem < 0 || elem > 1)
                throw exception::OutOfRange("qpp::Noise::Noise");
        // END EXCEPTION CHECKS

        d_ = Ks[0].rows(); // set the local dimension
        probs_ = probs;
    }

    /**
     * \brief Default virtual destructor
     */
    virtual ~NoiseBase() = default;

    // getters
    /**
     * \brief Local dimension
     *
     * \return Local dimension
     */
    idx get_d() const { return d_; };

    /**
     * \brief Vector of noise operators
     *
     * \return Vector of noise operators
     */
    std::vector<cmat> get_Ks() const { return Ks_; }

    /**
     * \brief Vector of probabilities corresponding to each noise operator
     *
     * \return Probability vector
     */
    std::vector<double> get_probs() const {
        if (generated_ ||
            std::is_same<NoiseType::StateIndependent, noise_type>::value) {
            return probs_;
        } else
            throw exception::CustomException(
                "qpp::Noise::get_probs()",
                "Noise::operator() was not yet invoked");
    }

    /**
     * \brief Index of the last occurring noise element
     *
     * \return Index of the last occurring noise element
     */
    idx get_last_idx() const {
        if (generated_) {
            return i_;
        } else
            throw exception::CustomException(
                "qpp::Noise::get_last_idx()",
                "Noise::operator() was not yet invoked");
    }

    /**
     * \brief Probability of the last occurring noise element
     *
     * \return Probability of the last occurring noise element
     */
    double get_last_p() const {
        if (generated_) {
            return probs_[i_];
        } else
            throw exception::CustomException(
                "qpp::Noise::get_last_p()",
                "Noise::operator() was not yet invoked");
    }

    /**
     * \brief Last occurring noise element
     *
     * \return Last occurring noise element
     */
    cmat get_last_K() const {
        if (generated_) {
            return Ks_[i_];
        } else
            throw exception::CustomException(
                "qpp::Noise::get_last_K()",
                "Noise::operator() was not yet invoked");
    }
    // end getters

    /**
     * \brief Function invocation operator, applies the underlying noise
     * model on qudit \a target of the multi-partite state vector or density
     * matrix \a state
     *
     * \param state Multi-partite state vector or density matrix
     * \param target Qudit index where the noise is applied
     * \return Resulting state vector or density matrix
     */
    virtual cmat operator()(const cmat& state, idx target) const {
        cmat result;
        try {
            compute_probs_(state, std::vector<idx>{target});
            result = compute_state_(state, std::vector<idx>{target});
        } catch (qpp::exception::Exception&) {
            std::cerr << "In qpp::Noise::operator()\n";
            throw;
        }

        return result;
    }

    /**
     * \brief Function invocation operator, applies the underlying correlated
     * noise model on qudits specified by \a target of the multi-partite state
     * vector or density matrix \a state
     *
     * \param state Multi-partite state vector or density matrix
     * \param target Qudit indexes where the correlated noise is applied
     * \return Resulting state vector or density matrix
     */
    virtual cmat operator()(const cmat& state,
                            const std::vector<idx>& target) const {
        cmat result;
        try {
            compute_probs_(state, target);
            result = compute_state_(state, target);
        } catch (qpp::exception::Exception&) {
            std::cerr << "In qpp::Noise::operator()\n";
            throw;
        }

        return result;
    }
}; /* class Noise */

// qubit noise models

/**
 * \class qpp::QubitDepolarizingNoise
 * \brief Qubit depolarizing noise
 */
class QubitDepolarizingNoise : public NoiseBase<NoiseType::StateIndependent> {
  public:
    /**
     * \brief Qubit depolarizing noise constructor
     *
     * \param p Noise probability
     */
    explicit QubitDepolarizingNoise(double p)
        : NoiseBase({Gates::get_instance().Id2, Gates::get_instance().X,
                     Gates::get_instance().Y, Gates::get_instance().Z},
                    {1 - p, p / 3, p / 3, p / 3}) {
        // EXCEPTION CHECKS

        if (p < 0 || p > 1)
            throw exception::OutOfRange(
                "qpp::QubitDepolarizingNoise::QubitDepolarizingNoise()");
        // END EXCEPTION CHECKS
    }
}; /* class QubitDepolarizingNoise */

/**
 * \class qpp::QubitPhaseFlipNoise
 * \brief Qubit phase flip (dephasing) noise
 */
class QubitPhaseFlipNoise : public NoiseBase<NoiseType::StateIndependent> {
  public:
    /**
     * \brief Qubit phase flip (dephasing) noise constructor
     *
     * \param p Noise probability
     */
    explicit QubitPhaseFlipNoise(double p)
        : NoiseBase({Gates::get_instance().Id2, Gates::get_instance().Z},
                    {1 - p, p}) {
        // EXCEPTION CHECKS

        if (p < 0 || p > 1)
            throw exception::OutOfRange(
                "qpp::QubitPhaseFlipNoise::QubitPhaseFlipNoise()");
        // END EXCEPTION CHECKS
    }
}; /* class QubitPhaseFlipNoise */

/**
 * \class qpp::QubitBitFlipNoise
 * \brief Qubit bit flip noise
 */
class QubitBitFlipNoise : public NoiseBase<NoiseType::StateIndependent> {
  public:
    /**
     * \brief Qubit bit flip noise constructor
     *
     * \param p Noise probability
     */
    explicit QubitBitFlipNoise(double p)
        : NoiseBase({Gates::get_instance().Id2, Gates::get_instance().X},
                    {1 - p, p}) {
        // EXCEPTION CHECKS

        if (p < 0 || p > 1)
            throw exception::OutOfRange(
                "qpp::QubitBitFlipNoise::QubitBitFlipNoise()");
        // END EXCEPTION CHECKS
    }
}; /* class QubitBitFlipNoise */

/**
 * \class qpp::QubitBitPhaseFlipNoise
 * \brief Qubit bit-phase flip (dephasing) noise
 */
class QubitBitPhaseFlipNoise : public NoiseBase<NoiseType::StateIndependent> {
  public:
    /**
     * \brief Qubit bit-phase flip noise constructor
     *
     * \param p Noise probability
     */
    explicit QubitBitPhaseFlipNoise(double p)
        : NoiseBase({Gates::get_instance().Id2, Gates::get_instance().Y},
                    {1 - p, p}) {
        // EXCEPTION CHECKS

        if (p < 0 || p > 1)
            throw exception::OutOfRange(
                "qpp::QubitBitPhaseFlipNoise::QubitBitPhaseFlipNoise()");
        // END EXCEPTION CHECKS
    }
}; /* class QubitBitPhaseFlipNoise */

/**
 * \class qpp::QubitAmplitudeDampingNoise
 * \brief Qubit amplitude damping noise, as described in Nielsen and Chuang
 */
class QubitAmplitudeDampingNoise : public NoiseBase<NoiseType::StateDependent> {
  public:
    /**
     * \brief Qubit amplitude damping noise constructor
     *
     * \param gamma Amplitude damping probability
     */
    explicit QubitAmplitudeDampingNoise(double gamma)
        : NoiseBase(std::vector<cmat>{
              ((cmat(2, 2)) << 1, 0, 0, std::sqrt(gamma)).finished(),
              ((cmat(2, 2)) << 0, std::sqrt(1 - gamma), 0, 0).finished()}) {
        // EXCEPTION CHECKS

        if (gamma < 0 || gamma > 1)
            throw exception::OutOfRange("qpp::QubitAmplitudeDampingNoise::"
                                        "QubitAmplitudeDampingNoise()");
        // END EXCEPTION CHECKS
    }
}; /* class QubitAmplitudeDampingNoise */

/**
 * \class qpp::QubitPhaseDampingNoise
 * \brief Qubit phase damping noise, as described in Nielsen and Chuang
 */
class QubitPhaseDampingNoise : public NoiseBase<NoiseType::StateDependent> {
  public:
    /**
     * \brief Qubit phase damping noise constructor
     *
     * \param gamma Phase damping probability
     */
    explicit QubitPhaseDampingNoise(double lambda)
        : NoiseBase(std::vector<cmat>{
              ((cmat(2, 2)) << 1, 0, 0, std::sqrt(1 - lambda)).finished(),
              ((cmat(2, 2)) << 0, 0, 0, std::sqrt(lambda)).finished()}) {
        // EXCEPTION CHECKS

        if (lambda < 0 || lambda > 1)
            throw exception::OutOfRange("qpp::QubitPhaseDampingNoise::"
                                        "QubitPhaseDampingNoise()");
        // END EXCEPTION CHECKS
    }
}; /* class QubitPhaseDampingNoise */

// qudit noise models

/**
 * \class qpp::QuditDepolarizingNoise
 * \brief Qudit depolarizing noise
 */
class QuditDepolarizingNoise : public NoiseBase<NoiseType::StateIndependent> {
    std::vector<cmat> fill_Ks_(idx d) const {
        std::vector<cmat> Ks(d * d);
        idx cnt = 0;
        for (idx i = 0; i < d; ++i)
            for (idx j = 0; j < d; ++j)
                Ks[cnt++] = powm(Gates::get_instance().Xd(d), i) *
                            powm(Gates::get_instance().Zd(d), j);

        return Ks;
    }
    std::vector<double> fill_probs_(double p, idx d) const {
        std::vector<double> probs(d * d);
        probs[0] = 1 - p;
        for (idx i = 1; i < d * d; ++i)
            probs[i] = p / (d - 1) * (d - 1);

        return probs;
    }
    const idx d_;

  public:
    /**
     * \brief Qudit depolarizing noise constructor
     *
     * \param p Noise probability
     * \param d Subsystem dimension
     */
    explicit QuditDepolarizingNoise(double p, idx d)
        : NoiseBase(fill_Ks_(d), fill_probs_(p, d)), d_{d} {
        // EXCEPTION CHECKS

        if (d < 2 || p < 0 || p > 1)
            throw exception::OutOfRange(
                "qpp::QuditDepolarizingNoise::QuditDepolarizingNoise()");
        // END EXCEPTION CHECKS
    }
}; /* class QuditDepolarizingNoise */

} /* namespace qpp */

#endif /* CLASSES_NOISE_H_ */
