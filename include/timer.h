/*
 * timer.h
 *
 *  Created on: Apr 1, 2014
 *      Author: vlad
 */

#ifndef TIMER_H_
#define TIMER_H_

#include <ctime>

namespace qpp
{

class Timer
{
	clock_t _start, _end;
	bool _toc;
public:
	Timer() :
			_start(clock()), _end(_start), _toc(false)
	{
	}

	void toc()
	{
		_end = clock();
		_toc = true;
	}

	void reset()
	{
		_start = _end = clock();
		_toc = false;
	}

	double ticks()
	{
		return _toc ? _end - _start : -1;
	}

	double secs()
	{
		return _toc ? ticks() / CLOCKS_PER_SEC : -1;
	}
	virtual ~Timer() = default;
};

} /* namespace qpp */

#endif /* TIMER_H_ */
