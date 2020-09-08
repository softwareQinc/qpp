#include "gtest/gtest.h"
#include "qpp.h"

using namespace qpp;

// Unit testing "classes/reversible.hpp"

/******************************************************************************/
/// BEGIN Bit_circuit& qpp::Bit_circuit::CNOT(const std::vector<idx>& pos)
TEST(qpp_Bit_circuit_CNOT, AllTests) {}
/******************************************************************************/
/// BEGIN Bit_circuit& qpp::Bit_circuit::FRED(const std::vector<idx>& pos)
TEST(qpp_Bit_circuit_FRED, AllTests) {}
/******************************************************************************/
/// BEGIN idx qpp::Bit_circuit::get_gate_count() const
TEST(qpp_Bit_circuit_get_gate_count, AllGates) {}
/******************************************************************************/
/// BEGIN idx qpp::Bit_circuit::get_gate_count(const std::string& name) const
TEST(qpp_Bit_circuit_get_gate_count, SpecificGate) {}
/******************************************************************************/
/// BEGIN idx qpp::Bit_circuit::get_gate_depth() const
TEST(qpp_Bit_circuit_get_gate_depth, AllGates) {}
/******************************************************************************/
/// BEGIN idx qpp::Bit_circuit::get_gate_depth(const std::string& name) const
TEST(qpp_Bit_circuit_get_gate_depth, SpecificGate) {}
/******************************************************************************/
/// BEGIN Bit_circuit& qpp::Bit_circuit::NOT(idx pos)
TEST(qpp_Bit_circuit_NOT, AllTests) {}
/******************************************************************************/
/// BEGIN Bit_circuit& qpp::Bit_circuit::reset() noexcept
TEST(qpp_Bit_circuit_reset, AllTests) {}
/******************************************************************************/
/// BEGIN Bit_circuit& qpp::Bit_circuit::SWAP(const std::vector<idx>& pos)
TEST(qpp_Bit_circuit_SWAP, AllTests) {}
/******************************************************************************/
/// BEGIN std::string qpp::Bit_circuit::to_JSON(
///       bool enclosed_in_curly_brackets = true) const override
TEST(qpp_Bit_circuit_to_JSON, AllTests) {}
/******************************************************************************/
/// BEGIN Bit_circuit& qpp::Bit_circuit::TOF(const std::vector<idx>& pos)
TEST(qpp_Bit_circuit_TOF, AllTests) {}
/******************************************************************************/
/// BEGIN Bit_circuit& qpp::Bit_circuit::X(idx pos)
TEST(qpp_Bit_circuit_X, AllTests) {}
/******************************************************************************/

/******************************************************************************/
/// BEGIN bool qpp::Dynamic_bitset::all() const noexcept
TEST(qpp_Dynamic_bitset_all, AllTests) {}
/******************************************************************************/
/// BEGIN bool qpp::Dynamic_bitset::any() const noexcept
TEST(qpp_Dynamic_bitset_any, AllTests) {}
/******************************************************************************/
/// BEGIN idx qpp::Dynamic_bitset::count() const noexcept
TEST(qpp_Dynamic_bitset_count, AllTests) {}
/******************************************************************************/
/// BEGIN const storage_type& qpp::Dynamic_bitset::data() const
TEST(qpp_Dynamic_bitset_data, AllTests) {}
/******************************************************************************/
/// BEGIN Dynamic_bitset& qpp::Dynamic_bitset::flip() noexcept
TEST(qpp_Dynamic_bitset_flip, AllBits) {}
/******************************************************************************/
/// BEGIN Dynamic_bitset& qpp::Dynamic_bitset::flip(idx pos)
TEST(qpp_Dynamic_bitset_flip, SpecificBit) {}
/******************************************************************************/
/// BEGIN bool qpp::Dynamic_bitset::get(idx pos) const noexcept
TEST(qpp_Dynamic_bitset_get, AllTests) {}
/******************************************************************************/
/// BEGIN bool qpp::Dynamic_bitset::none() const noexcept
TEST(qpp_Dynamic_bitset_none, AllTests) {}
/******************************************************************************/
/// BEGIN bool qpp::Dynamic_bitset::operator!=(const Dynamic_bitset& rhs) const
///       noexcept
TEST(qpp_Dynamic_bitset_operator_noneq, AllTests) {}
/******************************************************************************/
/// BEGIN idx qpp::Dynamic_bitset::operator-(const Dynamic_bitset& rhs) const
///       noexcept
TEST(qpp_Dynamic_bitset_operator_minus, AllTests) {}
/******************************************************************************/
/// BEGIN bool qpp::Dynamic_bitset::operator==(const Dynamic_bitset& rhs) const
///       noexcept
TEST(qpp_Dynamic_bitset_operator_eq, AllTests) {}
/******************************************************************************/
/// BEGIN Dynamic_bitset& qpp::Dynamic_bitset::rand(double p = 0.5)
TEST(qpp_Dynamic_bitset_rand_default, AllBits) {}
/******************************************************************************/
/// BEGIN Dynamic_bitset& qpp::Dynamic_bitset::rand(idx pos, double p = 0.5)
TEST(qpp_Dynamic_bitset_rand, SpecificBit) {}
/******************************************************************************/
/// BEGIN Dynamic_bitset& qpp::Dynamic_bitset::reset() noexcept
TEST(qpp_Dynamic_bitset_reset, AllBits) {}
/******************************************************************************/
/// BEGIN Dynamic_bitset& qpp::Dynamic_bitset::reset(idx pos)
TEST(qpp_Dynamic_bitset_reset, SpecificBit) {}
/******************************************************************************/
/// BEGIN Dynamic_bitset& qpp::Dynamic_bitset::set() noexcept
TEST(qpp_Dynamic_bitset_set, AllBits) {}
/******************************************************************************/
/// BEGIN Dynamic_bitset& qpp::Dynamic_bitset::set(idx pos, bool value = true)
TEST(qpp_Dynamic_bitset_set, SpecificBit) {}
/******************************************************************************/
/// BEGIN idx qpp::Dynamic_bitset::size() const noexcept
TEST(qpp_Dynamic_bitset_size, AllTests) {}
/******************************************************************************/
/// BEGIN idx qpp::Dynamic_bitset::storage_size() const noexcept
TEST(qpp_Dynamic_bitset_storage_size, AllTests) {}
/******************************************************************************/
/// BEGIN template <class CharT = char, class Traits = std::char_traits<CharT>,
///       class Allocator = std::allocator<CharT>>
///       std::basic_string<CharT, Traits, Allocator>
///       qpp::Dynamic_bitset::to_string(CharT zero = CharT('0'), CharT one =
///       CharT('1')) const
TEST(qpp_Dynamic_bitset_to_string, AllTests) {}
/******************************************************************************/
