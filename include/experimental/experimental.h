/*
 * Quantum++
 *
 * Copyright (c) 2013 - 2017 Vlad Gheorghiu (vgheorgh@gmail.com)
 *
 * This file is part of Quantum++.
 *
 * Quantum++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Quantum++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Quantum++.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
* \file experimental/experimental.h
* \brief Experimental/test functions/classes
*/

#ifndef EXPERIMENTAL_EXPERIMENTAL_H_
#define EXPERIMENTAL_EXPERIMENTAL_H_

namespace qpp
{
/**
* \namespace qpp::experimental
* \brief Experimental/test functions/classes, do not use or modify
*/
namespace experimental
{

class Dynamic_bitset
{
public:
    using value_type = unsigned int;
protected:
    idx storage_size_; // storage size
    idx num_bits_; // number of bits
    std::vector<value_type> v_; // storage space

    idx get_index_(idx pos) const
    {
        return pos / (sizeof(value_type) * CHAR_BIT);
    }

    value_type get_offset_(idx pos) const
    {
        return pos % (sizeof(value_type) * CHAR_BIT);
    }
public:
    // by default all bits are set to zero
    Dynamic_bitset(idx num_bits):
            storage_size_{num_bits / (sizeof(value_type) * CHAR_BIT) + 1},
            num_bits_{num_bits},
            v_(storage_size_, 0)
    {}

    ///////////////////// getters /////////////////////

    // get the storage
    const std::vector<value_type>& data() const
    {
        return v_;
    }

    // get the size in bits
    std::size_t get_bit_size() const
    {
        return num_bits_;
    }

    // get the storage size (in units of value_type)
    std::size_t get_storage_size() const
    {
        return storage_size_;
    }

    // count the bits set to true
    // (i.e. computes the Hamming distance)
    std::size_t count() const noexcept
    {
        std::size_t result = 0;
        for (idx i = 0 ; i < get_bit_size(); ++i)
        {
            if (this->get(i))
                ++result;
        }

        return result;
    }

    // returns the bit at position pos
    bool get(idx pos) const
    {
        return 1 & (v_[get_index_(pos)] >> get_offset_(pos));
    }

    // returns true if none of the bits are set
    bool none() const noexcept
    {
        bool result = true;
        for (idx i = 0; i < get_storage_size(); ++i)
        {
            if (v_[i])
            {
                return false;
            }
        }

        return result;
    }

    // returns true if all bits are set to true
    bool all() const noexcept
    {
        bool result = true;
        for (idx i = 0; i < get_storage_size(); ++i)
        {
            if (~v_[i])
            {
                return false;
            }
        }

        return result;
    }

    // returns true if any bit is set to true
    bool any() const noexcept
    {
        return !(this->none());
    }

    ///////////////////// setters /////////////////////

    // set the bit to a specific value
    Dynamic_bitset& set(idx pos, bool value = true)
    {
        value ?  v_[get_index_(pos)] |=
                         (1 << get_offset_(pos)) :
                v_[get_index_(pos)] &= ~(1 << get_offset_(pos));

//        v_[get_index_(pos)] &= ~(!value << get_offset_(pos));

        return *this;
    }

    // sets all bits to true
    Dynamic_bitset& set() noexcept
    {
        for (idx i = 0 ; i < get_storage_size(); ++i)
        {
            v_[i] = ~0;
        }

        return *this;
    }

    // set the bit at pos according to a random Bernoulli(p) distribution
    Dynamic_bitset& rand(idx pos, double p = 0.5)
    {
        std::random_device rd;
        std::mt19937 gen{rd()};
        std::bernoulli_distribution d{p};

        this->set(pos, d(gen));

        return *this;
    }

    // set all bits according to a random Bernoulli(p) distribution
    Dynamic_bitset& rand(double p = 0.5)
    {
        for (idx i = 0 ; i < get_bit_size(); ++i)
        {
            this->rand(i, p);
        }

        return *this;
    }

    // set the bit false
    Dynamic_bitset& reset(idx pos)
    {
        v_[get_index_(pos)] &= ~(1 << get_offset_(pos));

        return *this;
    }

    // set all bits to 0
    Dynamic_bitset& reset() noexcept
    {
        for (idx i = 0 ; i < get_storage_size(); ++i)
        {
            v_[i] = 0;
        }

        return *this;
    }

    // flips the bit
    Dynamic_bitset& flip(idx pos)
    {
        v_[get_index_(pos)] ^= 1 << (get_offset_(pos));

        return *this;
    }

    // flip all bits
    Dynamic_bitset& flip() noexcept
    {
        for (idx i = 0 ; i < get_storage_size(); ++i)
        {
            v_[i] = ~v_[i];
        }

        return *this;
    }

    ///////////////////// operators /////////////////////

    bool operator==(const Dynamic_bitset& rhs) const noexcept
    {
        assert(this->get_bit_size() == rhs.get_bit_size());
        bool result = true;
        idx n = std::min(this->get_storage_size(), rhs.get_storage_size());
        for (idx i = 0 ; i < n; ++i)
        {
            if (v_[i] != rhs.v_[i])
            {
                return false;
            }
        }

        return result;
    }

    bool operator!=(const Dynamic_bitset& rhs) const noexcept
    {
        return !(*this == rhs);
    }

    ///////////////////// input/output /////////////////////

    friend
    std::ostream& operator<<(std::ostream& os, const Dynamic_bitset& rhs)
    {
        for (idx i = rhs.get_bit_size(); i-- > 0;)
        {
            os << rhs.get(i);
        }

        return os;
    }
};

// reversible circuit tool
class Bit_circuit: public Dynamic_bitset
{
public:
    // gate counters
    struct Gate_count
    {
        // 1 bit gates
        idx NOT = 0;
        idx& X = NOT;

        // 2 bit gates
        idx CNOT = 0;
        idx SWAP = 0;

        // 3 bit gates
        idx FRED = 0;
        idx TOF = 0;
    } gate_count{};

    // inherit the constructor
    using Dynamic_bitset::Dynamic_bitset;

    // flips the bit
    Bit_circuit& X(idx pos)
    {
        this->flip(pos);
        ++gate_count.X;

        return *this;
    }

    // flips the bit
    Bit_circuit& NOT(idx pos)
    {
        this->flip(pos);
        ++gate_count.NOT;

        return *this;
    }

    // controlled-NOT control-target
    Bit_circuit& CNOT(const std::vector<idx>& pos)
    {
        v_[get_index_(pos[1])] ^=
                (1 & (v_[get_index_(pos[0])] >> get_offset_(pos[0])))
                        << get_offset_(pos[1]);
        ++gate_count.CNOT;

        return *this;
    }

    // Toffoli control-control-target
    Bit_circuit& TOF(const std::vector<idx>& pos)
    {
        v_[get_index_(pos[2])] ^=
                ( (1 & (v_[get_index_(pos[1])] >> get_offset_(pos[1])))
                  &
                  (1 & (v_[get_index_(pos[0])] >> get_offset_(pos[0])))
                ) << get_offset_(pos[2]);

        ++gate_count.TOF;
        return *this;
    }

    // SWAP 2 bits
    Bit_circuit& SWAP(const std::vector<idx>& pos)
    {
        if ( this->get(pos[0]) != this->get(pos[1]) )
        {
            this->X(pos[0]);
            this->X(pos[1]);
        }

        ++gate_count.SWAP;
        return *this;
    }

    // Fredkin gate
    Bit_circuit& FRED(const std::vector<idx>& pos)
    {
        if ( this->get(pos[0]) )
        {
            this->SWAP({pos[1], pos[2]});
        }

        ++gate_count.FRED;
        return *this;
    }

    // reset the circuit all zero, clear all gates
    Bit_circuit& reset() noexcept
    {
        gate_count.NOT = gate_count.X = 0;
        gate_count.CNOT = gate_count.SWAP = 0;
        gate_count.FRED = gate_count.TOF = 0;

        return *this;
    }
};

} /* namespace experimental */
} /* namespace qpp */

#endif /* EXPERIMENTAL_EXPERIMENTAL_H_ */
