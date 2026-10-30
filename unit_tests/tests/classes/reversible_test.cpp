#include "gtest/gtest.h"

#include "qpp/qpp.hpp"

using namespace qpp;

// Unit testing "qpp/classes/reversible.hpp"

/// BEGIN BitCircuit& BitCircuit::CNOT(idx ctrl, idx target)
TEST(qpp_BitCircuit_CNOT, AllTests) {
    BitCircuit c("100");
    c.CNOT(0, 1);
    EXPECT_EQ("110", c.to_string());
    c.CNOT(0, 1);
    EXPECT_EQ("100", c.to_string());
    EXPECT_EQ(2u, c.get_gate_count("CNOT"));
}

/// BEGIN BitCircuit& BitCircuit::FRED(idx i, idx j, idx k)
TEST(qpp_BitCircuit_FRED, AllTests) {
    BitCircuit c("110");
    c.FRED(0, 1, 2);
    EXPECT_EQ("101", c.to_string());
    BitCircuit inactive("010");
    inactive.FRED(0, 1, 2);
    EXPECT_EQ("010", inactive.to_string());
}

/// BEGIN idx BitCircuit::get_gate_count(std::optional<std::string> name =
///       std::nullopt) const
TEST(qpp_BitCircuit_get_gate_count, TotalGateCount) {
    BitCircuit c(3);
    EXPECT_EQ(0u, c.get_gate_count());
    c.NOT(0).NOT(1).CNOT(0, 2);
    EXPECT_EQ(3u, c.get_gate_count());
}

TEST(qpp_BitCircuit_get_gate_count, SpecificGateCount) {
    BitCircuit c(3);
    c.NOT(0).X(1).CNOT(0, 2).TOF(0, 1, 2);
    EXPECT_EQ(2u, c.get_gate_count("NOT"));
    EXPECT_EQ(2u, c.get_gate_count("X"));
    EXPECT_EQ(1u, c.get_gate_count("CNOT"));
    EXPECT_EQ(1u, c.get_gate_count("TOF"));
    EXPECT_EQ(0u, c.get_gate_count("FRED"));
}

/// BEGIN idx BitCircuit::get_gate_depth(std::optional<std::string> name =
///       std::nullopt) const
TEST(qpp_BitCircuit_get_gate_depth, TotalGateDepth) {
    BitCircuit c(2);
    EXPECT_EQ(0u, c.get_gate_depth());
    c.NOT(0);
    EXPECT_EQ(1u, c.get_gate_depth());
    c.NOT(0);
    EXPECT_EQ(2u, c.get_gate_depth());
}

TEST(qpp_BitCircuit_get_gate_depth, SpecificGateDepth) {
    BitCircuit c(3);
    c.NOT(0).NOT(0).CNOT(0, 1);
    EXPECT_EQ(2u, c.get_gate_depth("NOT"));
    EXPECT_EQ(2u, c.get_gate_depth("X"));
    EXPECT_EQ(1u, c.get_gate_depth("CNOT"));
    EXPECT_EQ(0u, c.get_gate_depth("SWAP"));
}

/// BEGIN BitCircuit& BitCircuit::NOT(idx i)
TEST(qpp_BitCircuit_NOT, AllTests) {
    BitCircuit c(3);
    c.NOT(1);
    EXPECT_EQ("010", c.to_string());
    c.NOT(1);
    EXPECT_EQ("000", c.to_string());
}

/// BEGIN BitCircuit& BitCircuit::reset() noexcept override
TEST(qpp_BitCircuit_reset, AllTests) {
    BitCircuit c("101");
    c.NOT(0).CNOT(0, 1);
    c.reset();
    EXPECT_EQ("000", c.to_string());
    EXPECT_EQ(0u, c.get_gate_count());
    EXPECT_EQ(0u, c.get_gate_depth());
}

/// BEGIN BitCircuit& BitCircuit::SWAP(idx i, idx j)
TEST(qpp_BitCircuit_SWAP, AllTests) {
    BitCircuit c("100");
    c.SWAP(0, 2);
    EXPECT_EQ("001", c.to_string());
    c.SWAP(0, 2);
    EXPECT_EQ("100", c.to_string());
}

/// BEGIN std::string BitCircuit::to_JSON(
///       bool enclosed_in_curly_brackets = true) const override
TEST(qpp_BitCircuit_to_JSON, AllTests) {
    BitCircuit c("101");
    const std::string json = c.to_JSON();
    EXPECT_EQ('{', json.front());
    EXPECT_EQ('}', json.back());
    EXPECT_NE(std::string::npos, json.find("\"bit state\" : \"101\""));
    EXPECT_NE(std::string::npos, json.find("\"Hamming weight\" : 2"));
    EXPECT_EQ(json.substr(1, json.size() - 2), c.to_JSON(false));
}

/// BEGIN std::string BitCircuit::to_string(char zero = '0', char one = '1')
///       const override
TEST(qpp_BitCircuit_to_string, AllTests) {
    BitCircuit c("1010");
    EXPECT_EQ("1010", c.to_string());
    EXPECT_EQ("X.X.", c.to_string('.', 'X'));
}

/// BEGIN BitCircuit& BitCircuit::TOF(idx i, idx j, idx k)
TEST(qpp_BitCircuit_TOF, AllTests) {
    BitCircuit c("110");
    c.TOF(0, 1, 2);
    EXPECT_EQ("111", c.to_string());
    c.TOF(0, 1, 2);
    EXPECT_EQ("110", c.to_string());
    BitCircuit inactive("100");
    inactive.TOF(0, 1, 2);
    EXPECT_EQ("100", inactive.to_string());
}

/// BEGIN BitCircuit& BitCircuit::X(idx i)
TEST(qpp_BitCircuit_X, AllTests) {
    BitCircuit c(2);
    c.X(0);
    EXPECT_EQ("10", c.to_string());
    EXPECT_EQ(1u, c.get_gate_count("NOT"));
    c.X(0);
    EXPECT_EQ("00", c.to_string());
}

/// BEGIN bool DynamicBitset::all() const noexcept
TEST(qpp_DynamicBitset_all, AllTests) {
    DynamicBitset bits(8);
    EXPECT_FALSE(bits.all());
    bits.set();
    EXPECT_TRUE(bits.all());
    bits.reset(3);
    EXPECT_FALSE(bits.all());
}

/// BEGIN bool DynamicBitset::any() const noexcept
TEST(qpp_DynamicBitset_any, AllTests) {
    DynamicBitset bits(5);
    EXPECT_FALSE(bits.any());
    bits.set(2);
    EXPECT_TRUE(bits.any());
    bits.reset(2);
    EXPECT_FALSE(bits.any());
}

/// BEGIN idx DynamicBitset::count() const noexcept
TEST(qpp_DynamicBitset_count, AllTests) {
    DynamicBitset bits("10101");
    EXPECT_EQ(3u, bits.count());
    bits.reset(2);
    EXPECT_EQ(2u, bits.count());
}

/// BEGIN const storage_type& DynamicBitset::data() const
TEST(qpp_DynamicBitset_data, AllTests) {
    DynamicBitset bits(5);
    EXPECT_EQ(bits.storage_size(), bits.data().size());
    EXPECT_EQ(0u, bits.data()[0]);
    bits.set(0).set(3);
    EXPECT_EQ(9u, bits.data()[0]);
}

/// BEGIN DynamicBitset& DynamicBitset::flip() noexcept
TEST(qpp_DynamicBitset_flip, AllBits) {
    DynamicBitset bits(5);
    bits.flip();
    for (idx i = 0; i < bits.size(); ++i) {
        EXPECT_TRUE(bits.get(i));
    }
    bits.flip();
    EXPECT_EQ("00000", bits.to_string());
}

/// BEGIN DynamicBitset& DynamicBitset::flip(idx pos)
TEST(qpp_DynamicBitset_flip, SpecificBit) {
    DynamicBitset bits(5);
    bits.flip(2);
    EXPECT_EQ("00100", bits.to_string());
    bits.flip(2);
    EXPECT_EQ("00000", bits.to_string());
}

/// BEGIN bool DynamicBitset::get(idx pos) const noexcept
TEST(qpp_DynamicBitset_get, AllTests) {
    DynamicBitset bits("10100");
    EXPECT_FALSE(bits.get(0));
    EXPECT_FALSE(bits.get(1));
    EXPECT_TRUE(bits.get(2));
    EXPECT_FALSE(bits.get(3));
    EXPECT_TRUE(bits.get(4));
}

/// BEGIN bool DynamicBitset::none() const noexcept
TEST(qpp_DynamicBitset_none, AllTests) {
    DynamicBitset bits(5);
    EXPECT_TRUE(bits.none());
    bits.set(4);
    EXPECT_FALSE(bits.none());
    bits.reset();
    EXPECT_TRUE(bits.none());
}

/// BEGIN bool DynamicBitset::operator!=(const DynamicBitset& rhs) const
///       noexcept
TEST(qpp_DynamicBitset_operator_noneq, AllTests) {
    DynamicBitset a("101"), b("100");
    EXPECT_TRUE(a != b);
    b.set(0);
    EXPECT_FALSE(a != b);
}

/// BEGIN idx DynamicBitset::operator-(const DynamicBitset& rhs) const
///       noexcept
TEST(qpp_DynamicBitset_operator_minus, AllTests) {
    DynamicBitset a("10110"), b("00101");
    EXPECT_EQ(3u, a - b);
    EXPECT_EQ(3u, b - a);
    EXPECT_EQ(0u, a - a);
}

/// BEGIN bool DynamicBitset::operator==(const DynamicBitset& rhs) const
///       noexcept
TEST(qpp_DynamicBitset_operator_eq, AllTests) {
    DynamicBitset a("101"), b("101"), c("100");
    EXPECT_TRUE(a == b);
    EXPECT_FALSE(a == c);
}

/// BEGIN DynamicBitset& DynamicBitset::rand(realT p = 0.5)
TEST(qpp_DynamicBitset_rand_default, AllBits) {
    DynamicBitset bits(20);
    bits.rand(1.);
    EXPECT_EQ(20u, bits.count());
    bits.rand(0.);
    EXPECT_TRUE(bits.none());
}

/// BEGIN DynamicBitset& DynamicBitset::rand(idx pos, realT p = 0.5)
TEST(qpp_DynamicBitset_rand, SpecificBit) {
    DynamicBitset bits(5);
    bits.rand(2, 1);
    EXPECT_EQ("00100", bits.to_string());
    bits.rand(2, 0);
    EXPECT_EQ("00000", bits.to_string());
}

/// BEGIN virtual DynamicBitset& DynamicBitset::reset() noexcept
TEST(qpp_DynamicBitset_reset, AllBits) {
    DynamicBitset bits("11111");
    bits.reset();
    EXPECT_TRUE(bits.none());
    EXPECT_EQ("00000", bits.to_string());
}

/// BEGIN DynamicBitset& DynamicBitset::reset(idx pos)
TEST(qpp_DynamicBitset_reset, SpecificBit) {
    DynamicBitset bits("111");
    bits.reset(1);
    EXPECT_EQ("101", bits.to_string());
}

/// BEGIN DynamicBitset& DynamicBitset::set() noexcept
TEST(qpp_DynamicBitset_set, AllBits) {
    DynamicBitset bits(5);
    bits.set();
    EXPECT_EQ(5u, bits.count());
    EXPECT_TRUE(bits.all());
}

/// BEGIN DynamicBitset& DynamicBitset::set(idx pos, bool value = true)
TEST(qpp_DynamicBitset_set, SpecificBit) {
    DynamicBitset bits(5);
    bits.set(3);
    EXPECT_EQ("01000", bits.to_string());
    bits.set(3, false);
    EXPECT_EQ("00000", bits.to_string());
}

/// BEGIN idx DynamicBitset::size() const noexcept
TEST(qpp_DynamicBitset_size, AllTests) {
    EXPECT_EQ(0u, DynamicBitset(0).size());
    EXPECT_EQ(5u, DynamicBitset(5).size());
    EXPECT_EQ(5u, DynamicBitset("10101").size());
}

/// BEGIN idx DynamicBitset::storage_size() const noexcept
TEST(qpp_DynamicBitset_storage_size, AllTests) {
    DynamicBitset bits(65);
    EXPECT_EQ(65u / (sizeof(DynamicBitset::value_type) * CHAR_BIT) + 1,
              bits.storage_size());
    EXPECT_EQ(bits.storage_size(), bits.data().size());
}

/// BEGIN virtual std::string DynamicBitset::to_string(char zero = '0',
///       char one = '1') const
TEST(qpp_DynamicBitset_to_string, AllTests) {
    DynamicBitset bits("1010");
    EXPECT_EQ("1010", bits.to_string());
    EXPECT_EQ("X.X.", bits.to_string('.', 'X'));
}
