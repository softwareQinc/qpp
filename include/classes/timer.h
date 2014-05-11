/*
 * timer.h
 *
 *  Created on: Apr 1, 2014
 *      Author: vlad
 */

#ifndef TIMER_H_
#define TIMER_H_

#include <chrono>
#include <ostream>

namespace qpp
{

class Timer
{
protected:
	std::chrono::steady_clock::time_point _start, _end;

public:
	Timer() :
			_start(std::chrono::steady_clock::now()), _end(_start)
	{
	}

	void tic()
	{
		_start = _end = std::chrono::steady_clock::now();
	}

	void toc()
	{
		_end = std::chrono::steady_clock::now();
	}

	double seconds() const
	{
		return std::chrono::duration_cast<std::chrono::duration<double>>(
				_end - _start).count();
	}

	friend std::ostream& operator<<(std::ostream &os, const Timer& rhs)
	{
		return os << rhs.seconds();
	}

	virtual ~Timer() = default;
};

} /* namespace qpp */

#endif /* TIMER_H_ */
