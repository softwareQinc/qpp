#include <vector>

#include "gtest/gtest.h"

#include "qpp/qpp.hpp"

using namespace qpp;

// Unit testing "classes/layouts.hpp"

/// BEGIN std::vector<idx> Lattice::get_dims() const override
TEST(qpp_Lattice_get_dims, AllTests) {}

/// BEGIN idx Lattice::operator()(const std::vector<idx>& xs) const override
TEST(qpp_Lattice_functor, Vector) {}

/// BEGIN template <typename... Ts> idx Lattice::operator()(Ts... xs) const
TEST(qpp_Lattice_functor, Variadic) {}

/// BEGIN std::vector<idx> Lattice::to_coordinates(idx i) const override
TEST(qpp_Lattice_to_coordinates, AllTests) {}

/// BEGIN idx PeriodicBoundaryLattice::operator()(
///       const std::vector<idx>& xs) const override
TEST(qpp_PeriodicBoundaryLattice_functor, AllTests) {}
