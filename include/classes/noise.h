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
 * \class qpp::INoise
 * \brief Base pure virtual class (interface) for all noise models, derive your
 * particular noise model from this class
 */
class INoise {
  protected:
    const idx d_;               ///< qudit dimension
    std::vector<double> probs_; ///< probabilities
    std::vector<cmat> Ks_;      ///< Kraus operators
  public:
    /**
     * \brief Constructs a noise instance
     *
     * \param d Qudit dimension
     * \param probs Probabilities of each noise (Kraus) operator
     * \param Ks Vector of noise (Kraus) operators that specify the noise
     */
    explicit INoise(idx d, const std::vector<double>& probs,
                    const std::vector<cmat>& Ks)
        : d_{d}, probs_{probs}, Ks_{Ks} {}

    /**
     * \brief Default virtual destructor
     */
    virtual ~INoise() = default;

    /**
     * \brief Conversion operator to a complex matrix
     *
     * \return Complex matrix
     */
    virtual operator cmat() const = 0;

    // getters
    /**
     * \brief Local dimension
     *
     * \return Local dimension
     */
    idx get_d() const { return d_; };

    /**
     * \brief Vector of probabilities corresponding to each noise operator
     *
     * \return Probability vector
     */
    std::vector<double> get_probs() const { return probs_; }

    /**
     * \brief Vector of noise operators
     *
     * \return Vector of noise operators
     */
    std::vector<cmat> get_Ks() const { return Ks_; }
    // end getters
}; /* class INoise */

// template method pattern
inline INoise::operator cmat() const {

    std::discrete_distribution<idx> dd{std::begin(probs_), std::end(probs_)};
    auto gen =
#ifdef NO_THREAD_LOCAL_
        RandomDevices::get_instance().get_prng();
#else
        RandomDevices::get_thread_local_instance().get_prng();
#endif
    return Ks_[dd(gen)];
} /* INoise::operator cmat() const */

/**
 * \class qpp::QubitDepolarizingNoise
 * \brief Qubit depolarizing noise
 */
class QubitDepolarizingNoise : INoise {
  public:
    /**
     * \brief Qubit depolarizing noise constructor
     *
     * \param p Noise probability
     */
    explicit QubitDepolarizingNoise(double p)
        : INoise(2, {1 - p, p / 3, p / 3, p / 3},
                 {Gates::get_instance().Id2, Gates::get_instance().X,
                  Gates::get_instance().Y, Gates::get_instance().Z}) {
        // EXCEPTION CHECKS

        if (p < 0 || p > 1)
            throw exception::OutOfRange(
                "qpp::QubitDepolarizingNoise::QubitDepolarizingNoise()");
        // END EXCEPTION CHECKS
    }

    /**
     * \brief Complex matrix conversion operator override
     *
     * \return Complex matrix
     */
    operator cmat() const override { return INoise::operator cmat(); }
}; /* class QubitDepolarizingNoise */

/**
 * \class qpp::QubitDephasingNoise
 * \brief Qubit dephasing noise
 */
class QubitDephasingNoise : INoise {
  public:
    /**
     * \brief Qubit dephasing noise constructor
     *
     * \param p Noise probability
     */
    explicit QubitDephasingNoise(double p)
        : INoise(2, {1 - p, p},
                 {Gates::get_instance().Id2, Gates::get_instance().Z}) {
        // EXCEPTION CHECKS

        if (p < 0 || p > 1)
            throw exception::OutOfRange(
                "qpp::QubitDephasingNoise::QubitDephasingNoise()");
        // END EXCEPTION CHECKS
    }

    /**
     * \brief Complex matrix conversion operator override
     *
     * \return Complex matrix
     */
    operator cmat() const override { return INoise::operator cmat(); }
}; /* class QubitDephasingNoise */

/**
 * \class qpp::QubitAmplitudeDampingNoise
 * \brief Qubit amplitude damping noise
 */
class QubitAmplitudeDampingNoise : INoise {
  public:
    /**
     * \brief Qubit amplitude damping noise constructor
     *
     * \param gamma Amplitude damping probability
     * \param psi Reference to state vector
     * \param i Qubit index
     */
    explicit QubitAmplitudeDampingNoise(double gamma, const ket& psi, idx i)
        : INoise(2, std::vector<double>(2, 0), std::vector<cmat>(2)) {
        // EXCEPTION CHECKS

        if (gamma < 0 || gamma > 1)
            throw exception::OutOfRange("qpp::QubitAmplitudeDampingNoise::"
                                        "QubitAmplitudeDampingNoise()");

        cmat K0(2, 2);
        cmat K1(2, 2);

        K0 << 1, 0, 0, std::sqrt(gamma);
        K1 << 0, std::sqrt(1 - gamma), 0, 0;

        double p = 0;

        try {
            idx n = internal::get_num_subsys(static_cast<idx>(psi.rows()), d_);
            cmat rho_i = ptrace(psi, {i}, std::vector<idx>(d_, n));
            p = trace(K0 * rho_i * adjoint(K0)).real();
            if (p < qpp::eps)
                p = 0;
            probs_[0] = p;
            probs_[1] = 1 - p;
            Ks_[0] = K0;
            Ks_[1] = K1;
        } catch (qpp::exception::Exception&) {
            std::cerr << "In "
                         "qpp::QubitAmplitudeDampingNoise::"
                         "QubitAmplitudeDampingNoise()\n";
            throw;
        }
        // END EXCEPTION CHECKS
    }

    /**
     * \brief Complex matrix conversion operator override
     *
     * \return Complex matrix
     */
    operator cmat() const override { return INoise::operator cmat(); }
}; /* class QubitAmplitudeDampingNoise */

/**
 * \class qpp::QubitPhaseDampingNoise
 * \brief Qubit phase damping noise
 */
class QubitPhaseDampingNoise : INoise {
  public:
    /**
     * \brief Qubit phase damping noise constructor
     *
     * \param gamma Phase damping probability
     * \param psi Reference to state vector
     * \param i Qubit index
     */
    explicit QubitPhaseDampingNoise(double lambda, const ket& psi, idx i)
        : INoise(2, std::vector<double>(2, 0), std::vector<cmat>(2)) {
        // EXCEPTION CHECKS

        if (lambda < 0 || lambda > 1)
            throw exception::OutOfRange("qpp::QubitPhaseDampingNoise::"
                                        "QubitPhaseDampingNoise()");

        cmat K0(2, 2);
        cmat K1(2, 2);

        K0 << 1, 0, 0, std::sqrt(1 - lambda);
        K1 << 0, 0, 0, std::sqrt(lambda);

        double p = 0;

        try {
            idx n = internal::get_num_subsys(static_cast<idx>(psi.rows()), d_);
            cmat rho_i = ptrace(psi, {i}, std::vector<idx>(d_, n));
            p = trace(K0 * rho_i * adjoint(K0)).real();
            if (p < qpp::eps)
                p = 0;
            probs_[0] = p;
            probs_[1] = 1 - p;
            Ks_[0] = K0;
            Ks_[1] = K1;
        } catch (qpp::exception::Exception&) {
            std::cerr << "In "
                         "qpp::QubitPhaseDampingNoise::"
                         "QubitPhaseDampingNoise()\n";
            throw;
        }
        // END EXCEPTION CHECKS
    }

    /**
     * \brief Complex matrix conversion operator override
     *
     * \return Complex matrix
     */
    operator cmat() const override { return INoise::operator cmat(); }
}; /* class QubitPhaseDampingNoise */

} /* namespace qpp */

#endif /* CLASSES_NOISE_H_ */
