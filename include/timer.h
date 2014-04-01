/*
 * timer.h
 *
 *  Created on: Apr 1, 2014
 *      Author: vlad
 */

#ifndef TIMER_H_
#define TIMER_H_

#include <chrono>

namespace qpp
{

class Timer
{
protected:
	std::chrono::high_resolution_clock::time_point _start, _end;
public:
	Timer() :
			_start(std::chrono::high_resolution_clock::now()), _end(_start)
	{
	}

	void toc()
	{
		_end = std::chrono::high_resolution_clock::now();
	}

	void tic()
	{
		_start = _end = std::chrono::high_resolution_clock::now();
	}

	double seconds()
	{
		return std::chrono::duration_cast<std::chrono::duration<double>>(
				_end - _start).count();
	}

	virtual ~Timer() = default;
};

} /* namespace qpp */

#endif /* TIMER_H_ */
