#include <chrono>
#include <thread> // for std::this_thread::sleep_for()
#include "gtest/gtest.h"
#include "qpp.h"

using namespace qpp;

// Unit testing "classes/timer.hpp"

// All test below test std::chrono::steady_clock timers,
// i.e Timer<T, CLOCK_T = std::chrono::steady_clock>

// the precision should be at least within 0.05s <=> 50 ms <=> 50000 micros
/******************************************************************************/
/// BEGIN template <typename U = T> U qpp::Timer::get_duration() const noexcept
TEST(qpp_Timer_get_duration, AllTests) {
    using namespace std::chrono;

    Timer<> t1;                              // default duration is double
    std::this_thread::sleep_for(seconds(1)); // sleep 1 second
    t1.toc();                                // get current time snap

    auto duration_t1_s = t1.get_duration();      // in seconds
    EXPECT_NEAR(duration_t1_s.count(), 1, 0.05); // within 0.05s

    auto duration_t1_ms = t1.get_duration<milliseconds>(); // in milli seconds
    EXPECT_NEAR(static_cast<double>(duration_t1_ms.count()), 1000,
                50); // within 0.05s

    Timer<std::chrono::microseconds> t2;               // in micro seconds
    std::this_thread::sleep_for(microseconds(100000)); // sleep 0.1 seconds
    t2.toc();                                          // get current time snap

    auto duration_t2_micros = t2.get_duration();
    EXPECT_NEAR(static_cast<double>(duration_t2_micros.count()), 100000,
                50000); // within 0.05s
}
/******************************************************************************/
/// BEGIN Timer& qpp::Timer::tic() noexcept
///
///       double qpp::Timer::tics() const noexcept
///
///       Timer& qpp::Timer::toc() noexcept
TEST(qpp_Timer_tic_tics_toc, AllTests) {
    using namespace std::chrono;

    Timer<> t;
    std::this_thread::sleep_for(milliseconds(100)); // sleep 0.1 seconds
    t.tic();                                        // reset the timer
    t.toc();                                        // get current time snap
    EXPECT_NEAR(t.tics(), 0, 0.05);                 // within 0.05s

    t.tic(); // reset
    // sleep an additional 0.1s (without reseting the timer)
    std::this_thread::sleep_for(milliseconds(100));
    t.toc();                          // get current time snap
    EXPECT_NEAR(t.tics(), 0.1, 0.05); // within 0.05s

    std::this_thread::sleep_for(milliseconds(100));
    t.toc();                          // get current time snap
    EXPECT_NEAR(t.tics(), 0.2, 0.05); // within 0.05s

    t.tic(); // reset
    std::this_thread::sleep_for(milliseconds(100));
    t.toc();                          // get current time snap
    EXPECT_NEAR(t.tics(), 0.1, 0.05); // within 0.05s
}
/******************************************************************************/
