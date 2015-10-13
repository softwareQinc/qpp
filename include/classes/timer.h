/*
 * Quantum++
 *
 * Copyright (c) 2013 - 2016 Vlad Gheorghiu (vgheorgh@gmail.com)
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
* \brief Chronometer
*
* \tparam T Tics duration, default is std::chrono::duration<double, 1>
* i.e. seconds in double precision
* \tparam CLOCK_T Clock's type, default is std::chrono::steady_clock,
* not affected by wall clock changes during runtime
*/
template<typename T = std::chrono::duration<double>,
        typename CLOCK_T = std::chrono::steady_clock>
class Timer : public IDisplay
{
protected:
    typename CLOCK_T::time_point _start, _end;

public:
    /**
    * \brief Constructs an instance with the current time
    * as the starting point
    */
    Timer() noexcept :
            _start{CLOCK_T::now()}, _end{_start}
    {
    }

    /**
    * \brief Resets the chronometer
    *
    * Resets the starting/ending point to the current time
    */
    void tic() noexcept
    {
        _start = _end = CLOCK_T::now();
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
        _end = CLOCK_T::now();
        return *this;
    }

    /**
    * \brief Time passed in the duration specified by T
    *
    * \return Number of tics (specified by T) that passed between the
    * instantiation/reset and invocation of qpp::Timer::toc()
    */
    double tics() const noexcept
    {
        return std::chrono::duration_cast<T>(_end - _start).count();
    }

    /**
    * \brief Duration specified by U
    *
    * \tparam U Duration, default is T, which defaults to
    * std::chrono::duration<double, 1>, i.e. seconds in double precision
    *
    * \return Duration that passed between the
    * instantiation/reset and invocation of qpp::Timer::toc()
    */
    template<typename U = T>
    U get_duration()
    {
        return std::chrono::duration_cast<U>(_end - _start);
    }

    /**
    * \brief Default copy constructor
    */
    Timer(const Timer&) = default;

    /**
    * \brief Default move constructor
    */
    Timer(Timer&&) = default;

    /**
    * \brief Default copy assignment operator
    */
    Timer& operator=(const Timer&) = default;

    /**
    * \brief Default move assignment operator
    */
    Timer& operator=(Timer&&) = default;

    /**
    * \brief Default virtual destructor
    */
    virtual ~Timer() = default;

private:
    /**
    * \brief qpp::IDisplay::display() override
    *
    * \param os Output stream
    * \return Writes to the output stream the number of tics (specified by T)
    * that passed between the instantiation/reset and invocation
    * of qpp::Timer::toc().
    */
    std::ostream& display(std::ostream& os) const override
    {
        return os << tics();
    }
}; /* class Timer */

} /* namespace qpp */

#endif /* CLASSES_TIMER_H_ */
