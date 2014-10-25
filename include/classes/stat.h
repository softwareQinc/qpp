/*
 * stat.h
 *
 *  Created on: Dec 17, 2013
 *      Author: vlad
 */

#ifndef STAT_H_
#define STAT_H_

// statistical distributions etc

namespace qpp
{

// light wrappers around C++11 statistical distributions

template<typename T = double>
class NormalDistribution
{
protected:
	std::normal_distribution<T> _d;

public:
	NormalDistribution(T mean = 0, T sigma = 1) :
			_d(std::normal_distribution<T>(mean, sigma))
	{
	}

	T sample()
	{
		return _d(RandomDevices::get_instance()._rng);
	}
};

template<typename T = double>
class UniformRealDistribution
{
protected:
	std::uniform_real_distribution<T> _d;

public:
	UniformRealDistribution(T a = 0, T b = 1) :
			_d(std::uniform_real_distribution<T>(a, b))
	{
	}

	T sample()
	{
		return _d(RandomDevices::get_instance()._rng);
	}
};

template<typename T = int>
class UniformIntegerDistribution
{
protected:
	std::uniform_int_distribution<T> _d;

public:
	UniformIntegerDistribution(T a = std::numeric_limits<T>::min(), T b =
			std::numeric_limits<T>::max()) :
			_d(std::uniform_int_distribution<T>(a, b))
	{
	}

	T sample()
	{
		return _d(RandomDevices::get_instance()._rng);
	}
};

template<typename T = std::size_t>
class DiscreteDistribution
{
protected:
	std::discrete_distribution<T> _d;

public:
	template<typename InputIterator>
	DiscreteDistribution(InputIterator first, InputIterator last) :
			_d(first, last)
	{
	}

	DiscreteDistribution(std::initializer_list<double> weights) :
			_d(weights)
	{
	}

	DiscreteDistribution(std::vector<double> weights) :
			_d(std::begin(weights), std::end(weights))
	{
	}

	T sample()
	{
		return _d(RandomDevices::get_instance()._rng);
	}

	std::vector<double> probabilities() const
	{
		return _d.probabilities();
	}
};

template<typename T = std::size_t>
class DiscreteDistributionAbsSquare
{
protected:
	std::discrete_distribution<T> _d;

	template<typename InputIterator>
	std::vector<double> cplx2weights(InputIterator first,
			InputIterator last) const
	{
		std::vector<double> weights(last - first);
		std::transform(first, last, std::begin(weights),
				[](const cplx &z)->double
				{	return std::pow(std::abs(z),2);});
		return weights;
	}

public:
	template<typename InputIterator>
	DiscreteDistributionAbsSquare(InputIterator first, InputIterator last) :
			_d { }
	{
		std::vector<double> weights = cplx2weights(first, last);
		std::discrete_distribution<T> tmp(std::begin(weights),
				std::end(weights));
		_d = tmp;
	}

	DiscreteDistributionAbsSquare(std::initializer_list<cplx> amplitudes) :
			_d { }
	{
		std::vector<double> weights = cplx2weights(std::begin(amplitudes),
				std::end(amplitudes));
		std::discrete_distribution<T> tmp(std::begin(weights),
				std::end(weights));
		_d = tmp;
	}

	DiscreteDistributionAbsSquare(std::vector<cplx> amplitudes) :
			_d { }
	{
		std::vector<double> weights = cplx2weights(std::begin(amplitudes),
				std::end(amplitudes));
		std::discrete_distribution<T> tmp(std::begin(weights),
				std::end(weights));
		_d = tmp;
	}

	template<typename Derived>
	DiscreteDistributionAbsSquare(const Eigen::MatrixBase<Derived> &V) :
			_d { }
	{
		const DynMat<typename Derived::Scalar> & rV = V;
		// check zero-size
		if (!internal::_check_nonzero_size(rV))
			throw Exception("DiscreteDistributionAbsSquare::"
					"DiscreteDistributionAbsSquare",
					Exception::Type::ZERO_SIZE);

		// check vector
		if (!internal::_check_vector(rV))
			throw Exception("DiscreteDistributionAbsSquare::"
					"DiscreteDistributionAbsSquare",
					Exception::Type::MATRIX_NOT_VECTOR);

		std::vector<double> weights = cplx2weights(rV.data(),
				rV.data() + rV.size());
		std::discrete_distribution<T> tmp(std::begin(weights),
				std::end(weights));
		_d = tmp;

	}

	T sample()
	{
		return _d(RandomDevices::get_instance()._rng);
	}

	std::vector<double> probabilities() const
	{
		return _d.probabilities();
	}
};

} /* namespace qpp */

#endif /* STAT_H_ */
