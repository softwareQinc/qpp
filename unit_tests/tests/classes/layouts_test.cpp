#include <vector>

#include "gtest/gtest.h"

#include "qpp/qpp.hpp"

using namespace qpp;

// Unit testing "qpp/classes/layouts.hpp"

/// BEGIN std::vector<idx> Lattice::get_dims() const override
TEST(qpp_Lattice_get_dims, AllTests) {
    const Lattice lattice({2, 3, 4});
    EXPECT_EQ((std::vector<idx>{2, 3, 4}), lattice.get_dims());

    const Lattice line(5);
    EXPECT_EQ((std::vector<idx>{5}), line.get_dims());
    EXPECT_THROW(Lattice(std::vector<idx>{}), exception::ZeroSize);
}

/// BEGIN idx Lattice::operator()(const std::vector<idx>& xs) const override
TEST(qpp_Lattice_functor, Vector) {
    const Lattice lattice({2, 3, 4});
    EXPECT_EQ(0u, lattice(std::vector<idx>{0, 0, 0}));
    EXPECT_EQ(6u, lattice(std::vector<idx>{0, 1, 2}));
    EXPECT_EQ(23u, lattice(std::vector<idx>{1, 2, 3}));

    EXPECT_THROW(lattice(std::vector<idx>{1, 2}), exception::SizeMismatch);
    EXPECT_THROW(lattice(std::vector<idx>{2, 0, 0}), exception::OutOfRange);
}

/// BEGIN template <typename... Ts> idx Lattice::operator()(Ts... xs) const
TEST(qpp_Lattice_functor, Variadic) {
    const Lattice lattice(2, 3, 4);
    EXPECT_EQ(0u, lattice(0, 0, 0));
    EXPECT_EQ(6u, lattice(0, 1, 2));
    EXPECT_EQ(23u, lattice(1, 2, 3));
    EXPECT_EQ(lattice(std::vector<idx>{1, 0, 2}), lattice(1, 0, 2));

    EXPECT_THROW(lattice(1, 2), exception::SizeMismatch);
    EXPECT_THROW(lattice(1, 3, 0), exception::OutOfRange);
}

/// BEGIN std::vector<idx> Lattice::to_coordinates(idx i) const override
TEST(qpp_Lattice_to_coordinates, AllTests) {
    const Lattice lattice({2, 3, 4});
    EXPECT_EQ((std::vector<idx>{0, 0, 0}), lattice.to_coordinates(0));
    EXPECT_EQ((std::vector<idx>{0, 1, 2}), lattice.to_coordinates(6));
    EXPECT_EQ((std::vector<idx>{1, 2, 3}), lattice.to_coordinates(23));

    for (idx i = 0; i < 24; ++i) {
        EXPECT_EQ(i, lattice(lattice.to_coordinates(i)));
    }
    EXPECT_THROW(lattice.to_coordinates(24), exception::OutOfRange);
}

/// BEGIN idx PeriodicBoundaryLattice::operator()(
///       const std::vector<idx>& xs) const override
TEST(qpp_PeriodicBoundaryLattice_functor, AllTests) {
    const PeriodicBoundaryLattice lattice({2, 3, 4});
    EXPECT_EQ(23u, lattice(std::vector<idx>{1, 2, 3}));
    EXPECT_EQ(23u, lattice(std::vector<idx>{3, 8, 7}));
    EXPECT_EQ(0u, lattice(std::vector<idx>{2, 3, 4}));
    EXPECT_EQ(23u, lattice(3, 8, 7));

    EXPECT_THROW(lattice(std::vector<idx>{1, 2}), exception::SizeMismatch);
}
