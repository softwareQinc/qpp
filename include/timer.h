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
	clock_t start;
public:
	Timer() :
			start(clock())
	{
	}

	void reset()
	{
		start = clock();
	}
	double ticks()
	{
		return clock() - start;
	}
	double secs()
	{
		return (double) (clock() - start) / CLOCKS_PER_SEC;
	}
	virtual ~Timer() = default;
};

} /* namespace qpp */

#endif /* TIMER_H_ */
