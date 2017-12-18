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
* \file experimental/experimental.h
* \brief Experimental/test functions/classes
*/

#ifndef EXPERIMENTAL_EXPERIMENTAL_H_
#define EXPERIMENTAL_EXPERIMENTAL_H_

#include <algorithm>
#include <cassert>
#include <climits>
#include <cstddef>
#include <random>
#include <utility>
#include <vector>

using idx = std::size_t;

namespace qpp
{
/**
* \namespace qpp::experimental
* \brief Experimental/test functions/classes, do not use or modify
*/
namespace experimental
{
/**
* \class qpp::Dynamic_bitset
* \brief Dynamic bitset class, allows the specification of the number of bits
* at runtime (unlike std::bitset<N>)
*/
class Dynamic_bitset
{
public:
    using value_type = unsigned int; ///< Type of the storage elements
    using storage_type = std::vector<value_type>; ///< Type of the storage
protected:
    idx storage_size_;          ///< Storage size
    idx N_;                     ///< Number of bits
    std::vector<value_type> v_; ///< Storage space

    /**
    * \brief Index of the \a pos bit in the storage space
    * \param pos Bit location
    * \return Index of the \a pos bit in the storage space
    */
    idx index_(idx pos) const
    {
        return pos / (sizeof(value_type) * CHAR_BIT);
    }

    /**
    * \brief Offset of the \a pos bit in the storage space relative to its index
    * \param pos Bit location
    * \return Offset of the \a pos bit in the storage space relative to its
    * index
    */
    idx offset_(idx pos) const
    {
        return pos % (sizeof(value_type) * CHAR_BIT);
    }

public:
    /**
    * \brief Constructor, initializes all bits to false (zero)
    * \param N Number of bits in the bitset
    */
    Dynamic_bitset(idx N) :
            storage_size_{N / (sizeof(value_type) * CHAR_BIT) + 1},
            N_{N},
            v_(storage_size_)
    {}

    /* getters */

    /**
    * \brief Raw storage space of the bitset
    * \return Const reference to the underlying storage space
    */
    const storage_type& data() const
    {
        return v_;
    }

    /**
    * \brief Number of bits stored in the bitset
    * \return Number of bits
    */
    idx size() const
    {
        return N_;
    }

    /**
    * \brief Size of the underlying storage space (in units of value_type,
    * unsigned int by default)
    * \return Size of the underlying storage space
    */
    idx storage_size() const
    {
        return storage_size_;
    }

    // count the bits set to true
    // (i.e. computes the Hamming distance)
    /**
    * \brief
    * \return
    */
    idx count() const noexcept
    {
        std::size_t result = 0;
        for (idx i = 0; i < size(); ++i)
        {
            if (this->get(i))
                ++result;
        }

        return result;
    }

    // returns the bit at position pos
    /**
    * \brief
    * \param pos
    * \return
    */
    bool get(idx pos) const
    {
        return 1 & (v_[index_(pos)] >> offset_(pos));
    }

    // returns true if none of the bits are set
    /**
    * \brief
    * \return
    */
    bool none() const noexcept
    {
        bool result = true;
        for (idx i = 0; i < storage_size(); ++i)
        {
            if (v_[i])
            {
                return false;
            }
        }

        return result;
    }

    // returns true if all bits are set to true
    /**
    * \brief
    * \return
    */
    bool all() const noexcept
    {
        bool result = true;
        for (idx i = 0; i < storage_size(); ++i)
        {
            if (~v_[i])
            {
                return false;
            }
        }

        return result;
    }

    // returns true if any bit is set to true
    /**
    * \brief
    * \return
    */
    bool any() const noexcept
    {
        return !(this->none());
    }

    /* setters */
    // set the bit to a specific value
    /**
    * \brief
    * \param pos
    * \param value
    * \return
    */
    Dynamic_bitset& set(idx pos, bool value = true)
    {
        value ? v_[index_(pos)] |=
                        (1 << offset_(pos)) :
                v_[index_(pos)] &= ~(1 << offset_(pos));

//        v_[index_(pos)] &= ~(!value << offset_(pos));

        return *this;
    }

    // sets all bits to true
    /**
    * \brief
    * \return
    */
    Dynamic_bitset& set() noexcept
    {
        for (idx i = 0; i < storage_size(); ++i)
        {
            v_[i] = ~0;
        }

        return *this;
    }

    // set the bit according to a random Bernoulli(p) distribution
    /**
    * \brief
    * \param pos
    * \param p
    * \return
    */
    Dynamic_bitset& rand(idx pos, double p = 0.5)
    {
        std::random_device rd;
        std::mt19937 gen{rd()};
        std::bernoulli_distribution d{p};

        this->set(pos, d(gen));

        return *this;
    }

    // set all bits according to a random Bernoulli(p) distribution
    /**
    * \brief
    * \param p
    * \return
    */
    Dynamic_bitset& rand(double p = 0.5)
    {
        for (idx i = 0; i < size(); ++i)
        {
            this->rand(i, p);
        }

        return *this;
    }

    // set the bit false
    /**
    * \brief
    * \param pos
    * \return
    */
    Dynamic_bitset& reset(idx pos)
    {
        v_[index_(pos)] &= ~(1 << offset_(pos));

        return *this;
    }

    // set all bits to 0
    /**
    * \brief
    * \return
    */
    Dynamic_bitset& reset() noexcept
    {
        for (idx i = 0; i < storage_size(); ++i)
        {
            v_[i] = 0;
        }

        return *this;
    }

    // flips the bit
    /**
    * \brief
    * \param pos
    * \return
    */
    Dynamic_bitset& flip(idx pos)
    {
        v_[index_(pos)] ^= 1 << (offset_(pos));

        return *this;
    }

    // flip all bits
    /**
    * \brief
    * \return
    */
    Dynamic_bitset& flip() noexcept
    {
        for (idx i = 0; i < storage_size(); ++i)
        {
            v_[i] = ~v_[i];
        }

        return *this;
    }

    /* operators */
    /**
    * \brief
    * \param rhs
    * \return
    */
    bool operator==(const Dynamic_bitset& rhs) const noexcept
    {
        assert(this->size() == rhs.size());
        bool result = true;
        idx n = std::min(this->storage_size(), rhs.storage_size());
        for (idx i = 0; i < n; ++i)
        {
            if (v_[i] != rhs.v_[i])
            {
                return false;
            }
        }

        return result;
    }

    /**
    * \brief
    * \param rhs
    * \return
    */
    bool operator!=(const Dynamic_bitset& rhs) const noexcept
    {
        return !(*this == rhs);
    }

    /* input/output */
    /**
    * \brief
    * \param os
    * \param rhs
    * \return
    */
    friend
    std::ostream& operator<<(std::ostream& os, const Dynamic_bitset& rhs)
    {
        for (idx i = rhs.size(); i-- > 0;)
        {
            os << rhs.get(i);
        }

        return os;
    }

    /**
    * \brief
    * \tparam CharT
    * \tparam Traits
    * \tparam Allocator
    * \param zero
    * \param one
    * \return
    */
    template<
            class CharT = char,
            class Traits = std::char_traits<CharT>,
            class Allocator = std::allocator<CharT>
    >
    std::basic_string<CharT, Traits, Allocator>
    to_string(CharT zero = CharT('0'), CharT one = CharT('1')) const
    {
        std::basic_string<CharT, Traits, Allocator> result;
        idx bitset_size = this->size();
        result.resize(bitset_size);

        for (idx i = bitset_size; i-- > 0;)
        {
            if (!this->get(i))
            {
                result[bitset_size - i - 1] = zero;
            } else
            {
                result[bitset_size - i - 1] = one;
            }
        }

        return result;
    }
};

/**
* \class qpp::Bit_circuit
* \brief Classical reversible circuit simulator
*/
class Bit_circuit : public Dynamic_bitset
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
        v_[index_(pos[1])] ^=
                (1 & (v_[index_(pos[0])] >> offset_(pos[0])))
                        << offset_(pos[1]);
        ++gate_count.CNOT;

        return *this;
    }

    // Toffoli control-control-target
    Bit_circuit& TOF(const std::vector<idx>& pos)
    {
        v_[index_(pos[2])] ^=
                ((1 & (v_[index_(pos[1])] >> offset_(pos[1])))
                 &
                 (1 & (v_[index_(pos[0])] >> offset_(pos[0])))
                ) << offset_(pos[2]);

        ++gate_count.TOF;
        return *this;
    }

    // SWAP 2 bits
    Bit_circuit& SWAP(const std::vector<idx>& pos)
    {
        if (this->get(pos[0]) != this->get(pos[1]))
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
        if (this->get(pos[0]))
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
