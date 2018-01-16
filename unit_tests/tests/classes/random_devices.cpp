/*
 * This file is part of Quantum++.
 *
 * MIT License
 *
 * Copyright (c) 2013 - 2018 Vlad Gheorghiu (vgheorgh@gmail.com)
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
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
TEST(qpp_RandomDevices_load_save, AllTests) {
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
