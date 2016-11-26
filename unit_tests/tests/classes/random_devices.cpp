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

#include <sstream>
#include "gtest/gtest.h"
#include "qpp.h"

using namespace qpp;

// Unit testing "classes/random_devices.h"

/******************************************************************************/
/// BEGIN qpp::RandomDevices::load(std::istream& is)
///
///       qpp::RandomDevices::save(std::ostream& os) const
TEST(qpp_RandomDevices_load_save, AllTests)
{
    // save the state of the PRNG
    std::stringstream ss;
    qpp::rdevs.save(ss);

    // generate some random matrix
    cmat A1 = rand<cmat>(4, 4);
    // generate some random index
    idx i1 = randidx(0, 100);
    // generate some random biging
    bigint b1 = rand(static_cast<bigint>(-100), 100);
    // finally generate some random double
    double d1 = rand(-100.0, 100.0);

    // load the state of the PRNG
    ss.seekg(0);
    qpp::rdevs.load(ss);

    // generate again some random matrix
    cmat A2 = rand<cmat>(4, 4);
    // generate again some random index
    idx i2 = randidx(0, 100);
    // generate again some random biging
    bigint b2 = rand(static_cast<bigint>(-100), 100);
    // finally generate again some random double
    double d2 = rand(-100.0, 100.0);

    // make sure we reproduce the randomness
    EXPECT_EQ(0, norm(A1 - A2));
    EXPECT_EQ(i1, i2);
    EXPECT_EQ(b1, b2);
    EXPECT_EQ(d1, d2);
}
/******************************************************************************/
