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

/**
* \file classes/reversible.h
* \brief Support for classical reversible circuits
*/

#ifndef CLASSES_REVERSIBLE_H
#define CLASSES_REVERSIBLE_H

namespace qpp {
/**
* \class qpp::Dynamic_bitset
* \brief Dynamic bitset class, allows the specification of the number of bits
* at runtime (unlike std::bitset<N>)
*/
class Dynamic_bitset : public IDisplay {
  public:
    using value_type = unsigned int; ///< Type of the storage elements
    using storage_type = std::vector<value_type>; ///< Type of the storage
  protected:
    idx storage_size_;          ///< Storage size
    idx N_;                     ///< Number of bits
    std::vector<value_type> v_; ///< Storage space

    /**
    * \brief Index of the \a pos bit in the storage space
    *
    * \param pos Bit location
    * \return Index of the \a pos bit in the storage space
    */
    idx index_(idx pos) const { return pos / (sizeof(value_type) * CHAR_BIT); }

    /**
    * \brief Offset of the \a pos bit in the storage space relative to its index
    *
    * \param pos Bit location
    * \return Offset of the \a pos bit in the storage space relative to its
    * index
    */
    idx offset_(idx pos) const { return pos % (sizeof(value_type) * CHAR_BIT); }

  public:
    /**
    * \brief Constructor, initializes all bits to false (zero)
    *
    * \param N Number of bits in the bitset
    */
    Dynamic_bitset(idx N)
        : storage_size_{N / (sizeof(value_type) * CHAR_BIT) + 1}, N_{N},
          v_(storage_size_) {}

    /* getters */

    /**
    * \brief Raw storage space of the bitset
    *
    * \return Const reference to the underlying storage space
    */
    const storage_type& data() const { return v_; }

    /**
    * \brief Number of bits stored in the bitset
    *
    * \return Number of bits stored in the bitset
    */
    idx size() const noexcept { return N_; }

    /**
    * \brief Size of the underlying storage space (in units of value_type,
    * unsigned int by default)
    *
    * \return Size of the underlying storage space
    */
    idx storage_size() const noexcept { return storage_size_; }

    /**
    * \brief Number of bits set to one in the bitset (Hamming weight)
    *
    * \return Hamming weight
    */
    idx count() const noexcept {
        std::size_t result = 0;
        for (idx i = 0; i < size(); ++i) {
            if (this->get(i))
                ++result;
        }

        return result;
    }

    /**
    * \brief The value of the bit at position \a pos
    *
    * \param pos Position in the bitset
    * \return The value of the bit at position \a pos
    */
    bool get(idx pos) const noexcept {
        return 1 & (v_[index_(pos)] >> offset_(pos));
    }

    /**
    * \brief Checks whether none of the bits are set
    *
    * \return True if none of the bits are set
    */
    bool none() const noexcept {
        bool result = true;
        for (idx i = 0; i < storage_size(); ++i) {
            if (v_[i]) {
                return false;
            }
        }

        return result;
    }

    /**
    * \brief Checks whether all bits are set
    *
    * \return True if all of the bits are set
    */
    bool all() const noexcept {
        bool result = true;
        for (idx i = 0; i < storage_size(); ++i) {
            if (~v_[i]) {
                return false;
            }
        }

        return result;
    }

    /**
    * \brief Checks whether any bit is set
    *
    * \return True if any of the bits is set
    */
    bool any() const noexcept { return !(this->none()); }

    /* setters */
    /**
    * \brief Sets the bit at position \a pos
    *
    * \param pos Position in the bitset
    * \param value Bit value
    * \return Reference to the current instance
    */
    Dynamic_bitset& set(idx pos, bool value = true) {
        value ? v_[index_(pos)] |= (1 << offset_(pos))
              : v_[index_(pos)] &= ~(1 << offset_(pos));

        //        v_[index_(pos)] &= ~(!value << offset_(pos));

        return *this;
    }

    /**
    * \brief Set all bits to true
    *
    * \return Reference to the current instance
    */
    Dynamic_bitset& set() noexcept {
        for (idx i = 0; i < storage_size(); ++i) {
            v_[i] = ~0;
        }

        return *this;
    }

    /**
    * \brief Sets the bit at position \a pos according to a Bernoulli(p)
    * distribution
    *
    * \param pos Position in the bitset
    * \param p Probability
    * \return Reference to the current instance
    */
    Dynamic_bitset& rand(idx pos, double p = 0.5) {
        std::random_device rd;
        std::mt19937 gen{rd()};
        std::bernoulli_distribution d{p};

        this->set(pos, d(gen));

        return *this;
    }

    /**
    * \brief Sets all bits according to a Bernoulli(p) distribution
    *
    * \param p Probability
    * \return Reference to the current instance
    */
    Dynamic_bitset& rand(double p = 0.5) {
        for (idx i = 0; i < size(); ++i) {
            this->rand(i, p);
        }

        return *this;
    }

    /**
    * \brief Sets the bit at position \a pos to false
    *
    * \param pos Position in the bitset
    * \return Reference to the current instance
    */
    Dynamic_bitset& reset(idx pos) {
        v_[index_(pos)] &= ~(1 << offset_(pos));

        return *this;
    }

    /**
    * \brief Sets all bits to false
    *
    * \return Reference to the current instance
    */
    Dynamic_bitset& reset() noexcept {
        for (idx i = 0; i < storage_size(); ++i) {
            v_[i] = 0;
        }

        return *this;
    }

    /**
    * \brief Flips the bit at position \a pos
    *
    * \param pos Position in the bitset
    * \return Reference to the current instance
    */
    Dynamic_bitset& flip(idx pos) {
        v_[index_(pos)] ^= 1 << (offset_(pos));

        return *this;
    }

    /**
    * \brief Flips all bits
    *
    * \return Reference to the current instance
    */
    Dynamic_bitset& flip() noexcept {
        for (idx i = 0; i < storage_size(); ++i) {
            v_[i] = ~v_[i];
        }

        return *this;
    }

    /* operators */
    /**
    * \brief Equality operator
    *
    * \param rhs Dynamic_bitset against which the equality is being tested
    * \return True if the bitsets are equal (bit by bit), false otherwise
    */
    bool operator==(const Dynamic_bitset& rhs) const noexcept {
        assert(this->size() == rhs.size());
        bool result = true;
        idx n = std::min(this->storage_size(), rhs.storage_size());
        for (idx i = 0; i < n; ++i) {
            if (v_[i] != rhs.v_[i]) {
                return false;
            }
        }

        return result;
    }

    /**
    * \brief Inequality operator
    *
    * \param rhs Dynamic_bitset against which the inequality is being tested
    * \return True if the bitsets are not equal (bit by bit), false otherwise
    */
    bool operator!=(const Dynamic_bitset& rhs) const noexcept {
        return !(*this == rhs);
    }

    /**
    * \brief Number of places the two bitsets differ (Hamming distance)
    *
    * \param rhs Dynamic_bitset against which the the Hamming distance is
    * computed
    * \return Hamming distance
    */
    idx operator-(const Dynamic_bitset& rhs) const noexcept {
        idx result = 0;
        for (idx i = 0; i < size(); ++i) {
            if (this->get(i) != rhs.get(i))
                ++result;
        }

        return result;
    }

    /* input/output */
    /**
    * \brief String representation
    *
    * \tparam CharT String character type
    * \tparam Traits String traits
    * \tparam Allocator String Allocator
    * \param zero Character representing the zero
    * \param one Character representing the one
    * \return The bitset as a string
    */
    template <class CharT = char, class Traits = std::char_traits<CharT>,
              class Allocator = std::allocator<CharT>>
    std::basic_string<CharT, Traits, Allocator>
    to_string(CharT zero = CharT('0'), CharT one = CharT('1')) const {
        std::basic_string<CharT, Traits, Allocator> result;
        idx bitset_size = this->size();
        result.resize(bitset_size);

        for (idx i = bitset_size; i-- > 0;) {
            if (!this->get(i)) {
                result[bitset_size - i - 1] = zero;
            } else {
                result[bitset_size - i - 1] = one;
            }
        }

        return result;
    }

  private:
    /**
    * \brief qpp::IDisplay::display() override, displays the bitset bit by bit
    *
    * \param os Output stream
    * \return Reference to the output stream
    */
    std::ostream& display(std::ostream& os) const override {
        idx size = this->size();
        for (idx i = size; i-- > 0;) {
            os << this->get(i);
        }

        return os;
    }
}; /* class Dynamic_bitset */

/**
* \class qpp::Bit_circuit
* \brief Classical reversible circuit simulator
*/
class Bit_circuit : public Dynamic_bitset {
  public:
    struct Gate_count {
        // 1 bit gates
        idx NOT = 0;
        idx& X = NOT;

        // 2 bit gates
        idx CNOT = 0;
        idx SWAP = 0;

        // 3 bit gates
        idx FRED = 0;
        idx TOF = 0;
    } gate_count{}; ///< Gate counters

    /**
    * \brief Inherited constructor
    */
    using Dynamic_bitset::Dynamic_bitset;

    /**
    * \brief Conversion constructor, used to initialize a qpp::Bit_circuit with
    * a qpp::Dynamic_bitset
    *
    * \param dynamic_bitset Dynamic bitset
    */
    Bit_circuit(const Dynamic_bitset& dynamic_bitset)
        : Dynamic_bitset{dynamic_bitset} {};

    /**
    * \brief Bit flip
    * \see qpp::Bit_circuit::NOT()
    *
    * \param pos Bit position in the circuit
    * \return Reference to the current instance
    */
    Bit_circuit& X(idx pos) {
        this->flip(pos);
        ++gate_count.X;

        return *this;
    }

    /**
    * \brief Bit flip
    * \see qpp::Bit_circuit::X()
    *
    * \param pos Bit position in the circuit
    * \return Reference to the current instance
    */
    Bit_circuit& NOT(idx pos) {
        this->flip(pos);
        ++gate_count.NOT;

        return *this;
    }

    /**
    * \brief Controlled-NOT
    *
    * \param pos Bit position in the circuit
    * \return Reference to the current instance
    */
    Bit_circuit& CNOT(const std::vector<idx>& pos) {
        v_[index_(pos[1])] ^= (1 & (v_[index_(pos[0])] >> offset_(pos[0])))
                              << offset_(pos[1]);
        ++gate_count.CNOT;

        return *this;
    }

    /**
    * \brief Toffoli gate
    *
    * \param pos Bit positions in the circuit, in the order
    * control-control-target
    * \return Reference to the current instance
    */
    Bit_circuit& TOF(const std::vector<idx>& pos) {
        v_[index_(pos[2])] ^= ((1 & (v_[index_(pos[1])] >> offset_(pos[1]))) &
                               (1 & (v_[index_(pos[0])] >> offset_(pos[0]))))
                              << offset_(pos[2]);
        ++gate_count.TOF;

        return *this;
    }

    /**
    * \brief Swap bits
    *
    * \param pos Bit positions in the circuit
    * \return Reference to the current instance
    */
    Bit_circuit& SWAP(const std::vector<idx>& pos) {
        if (this->get(pos[0]) != this->get(pos[1])) {
            this->X(pos[0]);
            this->X(pos[1]);
        }
        ++gate_count.SWAP;

        return *this;
    }

    /**
    * \brief Fredkin gate (Controlled-SWAP)
    *
    * \param pos Bit positions in the circuit, in the order
    * control-target-target
    * \return Reference to the current instance
    */
    Bit_circuit& FRED(const std::vector<idx>& pos) {
        if (this->get(pos[0])) {
            this->SWAP({pos[1], pos[2]});
        }
        ++gate_count.FRED;

        return *this;
    }

    /**
    * \brief Reset the circuit all zero, clear all gates
    *
    * \return Reference to the current instance
    */
    Bit_circuit& reset() noexcept {
        gate_count.NOT = gate_count.X = 0;
        gate_count.CNOT = gate_count.SWAP = 0;
        gate_count.FRED = gate_count.TOF = 0;
        Dynamic_bitset::reset();

        return *this;
    }
}; /* class Bit_circuit */

} /* namespace qpp */

#endif /* CLASSES_REVERSIBLE_H */
