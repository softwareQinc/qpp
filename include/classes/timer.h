/*
 * timer.h
 *
 *  Created on: Apr 1, 2014
 *      Author: vlad
 */

#ifndef TIMER_H_
#define TIMER_H_

namespace qpp
{

/**
 * \brief Measures time
 *
 * Uses a std::chrono::steady_clock.
 * It is not affected by wall clock changes during runtime.
 */
class Timer
{
protected:
	std::chrono::steady_clock::time_point _start, _end;

public:
	/**
	 * \brief Constructs the instance with the current time
	 * as the starting point
	 */
	Timer() :
			_start(std::chrono::steady_clock::now()), _end(_start)
	{
	}

	/**
	 * \brief Resets the chronometer
	 *
	 * Resets the starting/ending point to the current time
	 */
	void tic()
	{
		_start = _end = std::chrono::steady_clock::now();
	}

	/**
	 * \brief Stops the chronometer
	 *
	 * Set the current time as the ending point
	 */
	void toc()
	{
		_end = std::chrono::steady_clock::now();
	}

	/**
	 * \brief Time passed in seconds
	 *
	 * @return Number of seconds that passed between the instantiation/reset
	 * and invocation of qpp::Timer::toc()
	 */
	double seconds() const
	{
		return std::chrono::duration_cast<std::chrono::duration<double>>(
				_end - _start).count();
	}

	/**
	 * \brief Overload for std::ostream operators
	 *
	 * @param os Output stream
	 * @param rhs Timer instance
	 * @return Writes to the output stream the number of seconds that passed
	 * between the instantiation/reset and invocation of qpp::Timer::toc().
	 */
	friend std::ostream& operator<<(std::ostream &os, const Timer& rhs)
	{
		return os << rhs.seconds();
	}
}; /* class Timer */

} /* namespace qpp */

#endif /* TIMER_H_ */
