#include <sstream>

#include "gtest/gtest.h"

#include "qpp/qpp.hpp"

using namespace qpp;

// Unit testing "classes/random_devices.hpp"

/// BEGIN std::mt19937& RandomDevices::get_prng()
TEST(qpp_RandomDevices_get_prng, AllTests) {}

/// BEGIN RandomDevices::load(std::istream& is)
///
///       RandomDevices::save(std::ostream& os) const
TEST(qpp_RandomDevices_load_save, AllTests) {
    // save the state of the PRNG
    std::stringstream ss;
    rdevs.save(ss);

    // generate some random matrix
    cmat A1 = rand<cmat>(4, 4);
    // generate some random index
    idx i1 = randidx(0, 100);
    // generate some random biging
    bigint b1 = rand(static_cast<bigint>(-100), static_cast<bigint>(100));
    // finally generate some random realT
    realT d1 = rand(static_cast<realT>(-100.0), static_cast<realT>(100.0));

    // load the state of the PRNG
    ss.seekg(0);
    rdevs.load(ss);

    // generate again some random matrix
    cmat A2 = rand<cmat>(4, 4);
    // generate again some random index
    idx i2 = randidx(0, 100);
    // generate again some random biging
    bigint b2 = rand(static_cast<bigint>(-100), static_cast<bigint>(100));
    // finally generate again some random realT
    realT d2 = rand(static_cast<realT>(-100.0), static_cast<realT>(100.0));

    // make sure we reproduce the randomness
    EXPECT_EQ(0, norm(A1 - A2));
    EXPECT_EQ(i1, i2);
    EXPECT_EQ(b1, b2);
    EXPECT_EQ(d1, d2);
}
