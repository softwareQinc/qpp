#include "gtest/gtest.h"

#include "qpp/qpp.hpp"

using namespace qpp;

// Unit testing "qpp/classes/reversible.hpp"

/// BEGIN BitCircuit& BitCircuit::CNOT(const std::vector<idx>& pos)
TEST(qpp_BitCircuit_CNOT, AllTests) {}

/// BEGIN BitCircuit& BitCircuit::FRED(const std::vector<idx>& pos)
TEST(qpp_BitCircuit_FRED, AllTests) {}

/// BEGIN idx BitCircuit::get_gate_count(std::optional<std::string> name =
///       std::nullopt) const
TEST(qpp_BitCircuit_get_gate_count, TotalGateCount) {}
TEST(qpp_BitCircuit_get_gate_count, SpecificGateCount) {}

/// BEGIN idx BitCircuit::get_gate_depth(std::optional<std::string> name =
///       std::nullopt) const
TEST(qpp_BitCircuit_get_gate_depth, TotalGateDepth) {}
TEST(qpp_BitCircuit_get_gate_depth, SpecificGateDepth) {}

/// BEGIN BitCircuit& BitCircuit::NOT(idx pos)
TEST(qpp_BitCircuit_NOT, AllTests) {}

/// BEGIN BitCircuit& BitCircuit::clear() noexcept override
TEST(qpp_BitCircuit_reset, AllTests) {}

/// BEGIN BitCircuit& BitCircuit::SWAP(const std::vector<idx>& pos)
TEST(qpp_BitCircuit_SWAP, AllTests) {}

/// BEGIN std::string BitCircuit::to_JSON(
///       bool enclosed_in_curly_brackets = true) const override
TEST(qpp_BitCircuit_to_JSON, AllTests) {}

/// BEGIN std::string BitCircuit::to_string(char zero = '0', char one = '1')
///       const override
TEST(qpp_BitCircuit_to_string, AllTests) {}

/// BEGIN BitCircuit& BitCircuit::TOF(const std::vector<idx>& pos)
TEST(qpp_BitCircuit_TOF, AllTests) {}

/// BEGIN BitCircuit& BitCircuit::X(idx pos)
TEST(qpp_BitCircuit_X, AllTests) {}

/// BEGIN bool DynamicBitset::all() const noexcept
TEST(qpp_DynamicBitset_all, AllTests) {}

/// BEGIN bool DynamicBitset::any() const noexcept
TEST(qpp_DynamicBitset_any, AllTests) {}

/// BEGIN idx DynamicBitset::count() const noexcept
TEST(qpp_DynamicBitset_count, AllTests) {}

/// BEGIN const storage_type& DynamicBitset::data() const
TEST(qpp_DynamicBitset_data, AllTests) {}

/// BEGIN DynamicBitset& DynamicBitset::flip() noexcept
TEST(qpp_DynamicBitset_flip, AllBits) {}

/// BEGIN DynamicBitset& DynamicBitset::flip(idx pos)
TEST(qpp_DynamicBitset_flip, SpecificBit) {}

/// BEGIN bool DynamicBitset::get(idx pos) const noexcept
TEST(qpp_DynamicBitset_get, AllTests) {}

/// BEGIN bool DynamicBitset::none() const noexcept
TEST(qpp_DynamicBitset_none, AllTests) {}

/// BEGIN bool DynamicBitset::operator!=(const DynamicBitset& rhs) const
///       noexcept
TEST(qpp_DynamicBitset_operator_noneq, AllTests) {}

/// BEGIN idx DynamicBitset::operator-(const DynamicBitset& rhs) const
///       noexcept
TEST(qpp_DynamicBitset_operator_minus, AllTests) {}

/// BEGIN bool DynamicBitset::operator==(const DynamicBitset& rhs) const
///       noexcept
TEST(qpp_DynamicBitset_operator_eq, AllTests) {}

/// BEGIN DynamicBitset& DynamicBitset::rand(realT p = 0.5)
TEST(qpp_DynamicBitset_rand_default, AllBits) {}

/// BEGIN DynamicBitset& DynamicBitset::rand(idx pos, realT p = 0.5)
TEST(qpp_DynamicBitset_rand, SpecificBit) {}

/// BEGIN virtual DynamicBitset& DynamicBitset::reset() noexcept
TEST(qpp_DynamicBitset_reset, AllBits) {}

/// BEGIN DynamicBitset& DynamicBitset::reset(idx pos)
TEST(qpp_DynamicBitset_reset, SpecificBit) {}

/// BEGIN DynamicBitset& DynamicBitset::set() noexcept
TEST(qpp_DynamicBitset_set, AllBits) {}

/// BEGIN DynamicBitset& DynamicBitset::set(idx pos, bool value = true)
TEST(qpp_DynamicBitset_set, SpecificBit) {}

/// BEGIN idx DynamicBitset::size() const noexcept
TEST(qpp_DynamicBitset_size, AllTests) {}

/// BEGIN idx DynamicBitset::storage_size() const noexcept
TEST(qpp_DynamicBitset_storage_size, AllTests) {}

/// BEGIN virtual std::string DynamicBitset::to_string(char zero = '0',
///       char one = '1') const
TEST(qpp_DynamicBitset_to_string, AllTests) {}
