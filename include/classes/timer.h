/*
 * Quantum++
 *
 * Copyright (c) 2013 - 2015 Vlad Gheorghiu (vgheorgh@gmail.com)
 *
 * This file is part of Quantum++.
 *
 * Quantum++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Quantum++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Quantum++.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
* \file classes/timer.h
* \brief Timing
*/

#ifndef CLASSES_TIMER_H_
#define CLASSES_TIMER_H_

namespace qpp
{

/**
* \class qpp::Timer
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
    * \brief Constructs an instance with the current time
    * as the starting point
    */
    Timer() noexcept :
            _start{std::chrono::steady_clock::now()}, _end{_start}
    {
    }

    /**
    * \brief Resets the chronometer
    *
    * Resets the starting/ending point to the current time
    */
    void tic() noexcept
    {
        _start = _end = std::chrono::steady_clock::now();
    }

    /**
    * \brief Stops the chronometer
    *
    * Set the current time as the ending point
    *
    * \return Current instance
    */
    const Timer& toc() noexcept
    {
        _end = std::chrono::steady_clock::now();
        return *this;
    }

    /**
    * \brief Time passed in seconds
    *
    * \return Number of seconds that passed between the instantiation/reset
    * and invocation of qpp::Timer::toc()
    */
    double seconds() const noexcept
    {
        return std::chrono::duration_cast<std::chrono::duration<double>>(
                _end - _start).count();
    }

    /**
    * \brief Overload for std::ostream operators
    *
    * \param os Output stream
    * \param rhs Timer instance
    * \return Writes to the output stream the number of seconds that passed
    * between the instantiation/reset and invocation of qpp::Timer::toc().
    */
    friend std::ostream& operator<<(std::ostream& os, const Timer& rhs)
    {
        return os << rhs.seconds();
    }

    /**
    * \brief Default virtual destructor
    */
    virtual ~Timer() = default;
}; /* class Timer */

} /* namespace qpp */

#endif /* CLASSES_TIMER_H_ */
