/*
 * Quantum++
 *
 * Copyright (c) 2013 - 2017 Vlad Gheorghiu (vgheorgh@gmail.com)
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

#include <chrono>
#include <thread> // for std::this_thread::sleep_for()
#include "gtest/gtest.h"
#include "qpp.h"

using namespace qpp;

// Unit testing "classes/timer.h"

// All test below test std::chrono::steady_clock timers,
// i.e Timer<T, CLOCK_T = std::chrono::steady_clock>

// the precision should be at least within 0.02s, i.e. 20 ms or 2000 micros
/******************************************************************************/
/// BEGIN template<typename U = T> U qpp::Timer::get_duration() const noexcept
TEST(qpp_Timer_get_duration, AllTests)
{
    using namespace std::chrono;

    Timer<> t1; // default duration is double
    std::this_thread::sleep_for(seconds(1)); // sleep 1 second
    t1.toc(); // get current time snap

    auto duration_t1_s = t1.get_duration(); // in seconds
    EXPECT_NEAR(duration_t1_s.count(), 1, 0.02); // within 0.02s

    auto duration_t1_ms = t1.get_duration<milliseconds>(); // in milli seconds
    EXPECT_NEAR(duration_t1_ms.count(), 1000, 20); // within 0.02s

    Timer<std::chrono::microseconds> t2; // in micro seconds
    std::this_thread::sleep_for(microseconds(100000)); // sleep 0.1 seconds
    t2.toc(); // get current time snap

    auto duration_t2_micros = t2.get_duration();
    EXPECT_NEAR(duration_t2_micros.count(), 100000, 20000); // within 0.02s
}
/******************************************************************************/
/// BEGIN void qpp::Timer::tic() noexcept
///
///       double qpp::Timer::tics() const noexcept
///
///       double qpp::Timer::toc() const noexcept
TEST(qpp_Timer_tic_tics_toc, AllTests)
{
    using namespace std::chrono;

    Timer<> t;
    std::this_thread::sleep_for(milliseconds(100)); // sleep 0.1 seconds
    t.tic(); // reset the timer
    t.toc(); // get current time snap
    EXPECT_NEAR(t.tics(), 0, 0.02); // within 0.02s

    t.tic(); // reset
    // sleep an additional 0.1s (without reseting the timer)
    std::this_thread::sleep_for(milliseconds(100));
    t.toc(); // get current time snap
    EXPECT_NEAR(t.tics(), 0.1, 0.02); // within 0.02s

    std::this_thread::sleep_for(milliseconds(100));
    t.toc(); // get current time snap
    EXPECT_NEAR(t.tics(), 0.2, 0.02); // within 0.02s

    t.tic(); // reset
    std::this_thread::sleep_for(milliseconds(100));
    t.toc(); // get current time snap
    EXPECT_NEAR(t.tics(), 0.1, 0.02); // within 0.02s
}
/******************************************************************************/
