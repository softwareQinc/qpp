/*
 * stat.h
 *
 *  Created on: Dec 17, 2013
 *      Author: vlad
 */

#ifndef STAT_H_
#define STAT_H_

#include <random>

// statistical distributions etc

namespace qpp {
    namespace stat {

        extern std::random_device _rd; // use for seeding
        extern std::mt19937 _rng; // our random number generator

        // light wrappers around C++11 statistical distributions

        class NormalDistribution {
        public:
            std::normal_distribution<> _d;

            NormalDistribution(double mean = 0, double sigma = 1) :
            _d(std::normal_distribution<>(mean, sigma)) {
            }

            double sample() {
                return _d(_rng);
            }
        };

        class UniformRealDistribution {
        public:
            std::uniform_real_distribution<> _d;

            UniformRealDistribution(double a = 0, double b = 1) :
            _d(std::uniform_real_distribution<>(a, b)) {
            }

            double sample() {
                return _d(_rng);
            }
        };

        class DiscreteDistribution {
        public:
            std::discrete_distribution<size_t> _d;

            DiscreteDistribution(std::initializer_list<double> weights) :
            _d(weights) {
            }

            template<typename InputIterator>
            DiscreteDistribution(InputIterator first, InputIterator last) :
            _d(first, last) {
            }

            DiscreteDistribution(std::vector<double> weights) :
            _d(weights.begin(), weights.end()) {
            }

            size_t sample() {
                return _d(_rng);
            }
        };

        //TODO: Add constructor for cmat

        class DiscreteDistributionFromComplex {
        public:
            std::discrete_distribution<size_t> _d;

            DiscreteDistributionFromComplex(
                    std::initializer_list<types::cplx> amplitudes) :
            _d() {
                std::vector<double> weights;
                for (auto i : amplitudes)
                    weights.push_back(std::abs(i));
                std::discrete_distribution<size_t> tmp(weights.begin(), weights.end());
                _d = tmp;
            }

            template<typename InputIterator>
            DiscreteDistributionFromComplex(InputIterator first, InputIterator last) {
            }

            DiscreteDistributionFromComplex(std::vector<double> weights) :
            _d(weights.begin(), weights.end()) {
            }

            size_t sample() {
                return _d(_rng);
            }
        };

    }
}

#endif /* STAT_H_ */
