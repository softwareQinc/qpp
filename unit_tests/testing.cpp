#include <gtest/gtest.h>

#include <qpp.h>

using namespace qpp;

// TODO: write your unit tests here

// ********** qpp::sum() **********
TEST(qpp_sum_test, PositiveNumbers)
{
    std::vector<int> v{0, 1, 2, 3};
    EXPECT_EQ (6, qpp::sum(v.begin(), v.end()));
}

TEST(qpp_sum_test, NegativeNumbers)
{
    std::vector<int> v{0, -1, -2, -3};
    EXPECT_EQ (-6, qpp::sum(v.begin(), v.end()));
}

TEST(qpp_sum_test, MixedNumbers)
{
    std::vector<int> v{ -3, -2, -1, 0, 1, 2};
    EXPECT_EQ (-3, qpp::sum(v.begin(), v.end()));
}
// ********** END qpp::sum() **********

// ********** qpp::prod() **********
TEST(qpp_prod_test, PositiveNumbers)
{
    std::vector<int> v{1, 2, 3, 4};
    EXPECT_EQ (24, qpp::prod(v.begin(), v.end()));
}
// ********** END qpp::prod() **********

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
