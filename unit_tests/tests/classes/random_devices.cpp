#include <sstream>
#include "gtest/gtest.h"
#include "qpp.h"

using namespace qpp;

// Unit testing "classes/random_devices.h"

/******************************************************************************/
/// BEGIN std::mt19937& qpp::RandomDevices::get_prng()
TEST(qpp_RandomDevices_get_prng, AllTests) {}
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
