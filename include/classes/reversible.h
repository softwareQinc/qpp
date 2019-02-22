/*
 * This file is part of Quantum++.
 *
 * MIT License
 *
 * Copyright (c) 2013 - 2019 Vlad Gheorghiu (vgheorgh@gmail.com)
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
    using value_type = unsigned int; ///< type of the storage elements
    using storage_type = std::vector<value_type>; ///< type of the storage
  protected:
    idx storage_size_;          ///< storage size
    idx n_;                     ///< number of bits
    std::vector<value_type> v_; ///< storage space

    /**
     * \brief Index of the \a pos bit in the storage space
     *
     * \param pos Bit location
     * \return Index of the \a pos bit in the storage space
     */
    idx index_(idx pos) const { return pos / (sizeof(value_type) * CHAR_BIT); }

    /**
     * \brief Offset of the \a pos bit in the storage space relative to its
     * index
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
     * \param n Number of bits in the bitset
     */
    explicit Dynamic_bitset(idx n)
        : storage_size_{n / (sizeof(value_type) * CHAR_BIT) + 1}, n_{n},
          v_(storage_size_) {}

    /**
     * \brief Default virtual destructor
     */
    virtual ~Dynamic_bitset() = default;

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
    idx size() const noexcept { return n_; }

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
        idx result = 0;
        idx bitset_size = this->size();
        for (idx i = 0; i < bitset_size; ++i) {
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
        idx bitset_storage_size = this->storage_size();
        for (idx i = 0; i < bitset_storage_size; ++i) {
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
        idx bitset_storage_size = this->storage_size();
        for (idx i = 0; i < bitset_storage_size; ++i) {
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
        idx bitset_storage_size = this->storage_size();
        for (idx i = 0; i < bitset_storage_size; ++i) {
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
        idx bitset_size = this->size();
        for (idx i = 0; i < bitset_size; ++i) {
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
        idx bitset_storage_size = this->storage_size();
        for (idx i = 0; i < bitset_storage_size; ++i) {
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
        idx bitset_storage_size = this->storage_size();
        for (idx i = 0; i < bitset_storage_size; ++i) {
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
        idx bitset_size = this->size();
        for (idx i = 0; i < bitset_size; ++i) {
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
     * \param os Output stream passed by reference
     * \return Reference to the output stream
     */
    std::ostream& display(std::ostream& os) const override {
        idx bitset_size = this->size();
        for (idx i = bitset_size; i-- > 0;) {
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
    /**
     * \class qpp::Bit_circuit::Gate_count
     * \brief Gate counters
     * \see qpp::Bit_circuit::Gate_depth
     */
    class Gate_count {
      public:
        /**
         * \brief Resets the gate count
         */
        void reset() { NOT = CNOT = SWAP = FRED = TOF = total = 0; }
        // 1 bit gates
        idx NOT = 0;
        idx& X = NOT;

        // 2 bit gates
        idx CNOT = 0;
        idx SWAP = 0;

        // 3 bit gates
        idx FRED = 0;
        idx TOF = 0;

        // total count
        idx total = 0;
    } gate_count{}; ///< gate counters

    /**
     * \class qpp::Bit_circuit::Gate_depth
     * \brief Gate depth
     * \see qpp::Bit_circuit::Gate_counter
     */
    class Gate_depth {
        friend Bit_circuit;
        const Bit_circuit* bit_circuit_{};
        Dynamic_bitset bNOT, bCNOT, bSWAP, bFRED, bTOF, btotal;

      public:
        /**
         * \brief Constructs the gate depth class out of a reversible bit
         * circuit
         *
         * \note The reversible bit circuit must be an lvalue
         * \see qpp::Gate_depth(Gate_depth&&)
         *
         * \param bit_circuit Reversible bit circuit
         */
        explicit Gate_depth(const Bit_circuit& bit_circuit)
            : bit_circuit_{std::addressof(bit_circuit)},
              bNOT{bit_circuit.size()}, bCNOT{bit_circuit.size()},
              bSWAP{bit_circuit.size()}, bFRED{bit_circuit.size()},
              bTOF{bit_circuit.size()}, btotal{bit_circuit.size()} {}

        // silence -Weffc++ class has pointer data members
        /**
         * \brief Default copy constructor
         */
        Gate_depth(const Gate_depth&) = default;

        // silence -Weffc++ class has pointer data members
        /**
         * \brief Default copy assignment operator
         *
         * \return Reference to the current instance
         */
        Gate_depth& operator=(const Gate_depth&) = default;

        /**
         * \brief Disables rvalue QCircuit
         */
        Gate_depth(Gate_depth&&) = delete;

        /**
         * \brief Resets the gate depth
         */
        void reset() {
            NOT = CNOT = SWAP = FRED = TOF = total = 0;
            idx n = bit_circuit_->size();
            bNOT = bCNOT = bSWAP = bFRED = bTOF = btotal = Dynamic_bitset(n);
        }

        // 1 bit gates
        idx NOT = 0;
        idx& X = NOT;

        // 2 bit gates
        idx CNOT = 0;
        idx SWAP = 0;

        // 3 bit gates
        idx FRED = 0;
        idx TOF = 0;

        // total depth
        idx total = 0;
    } gate_depth{*this}; ///< gate depths

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
    explicit Bit_circuit(const Dynamic_bitset& dynamic_bitset)
        : Dynamic_bitset{dynamic_bitset} {};

    /**
     * \brief Bit flip
     * \see qpp::Bit_circuit::NOT()
     *
     * \param i Bit position in the circuit
     * \return Reference to the current instance
     */
    Bit_circuit& X(idx i) {
        this->NOT(i);

        return *this;
    }

    /**
     * \brief Bit flip
     * \see qpp::Bit_circuit::X()
     *
     * \param i Bit position in the circuit
     * \return Reference to the current instance
     */
    Bit_circuit& NOT(idx i) {
        this->flip(i);
        ++gate_count.NOT;
        ++gate_count.total;

        // compute the depth
        if (gate_count.NOT == 1)
            gate_depth.NOT = 1;
        // apply the gate
        gate_depth.bNOT.flip(i);
        // check whether gates overlap
        if (gate_depth.bNOT.get(i) == false) {
            // reset the b vector
            gate_depth.bNOT = Dynamic_bitset{n_};
            // set to true the locations of the last gate
            gate_depth.bNOT.set(i);
            ++gate_depth.NOT;
        }

        // compute the total depth
        if (gate_count.total == 1)
            gate_depth.total = 1;
        // apply the gate
        gate_depth.btotal.flip(i);
        // check whether gates overlap
        if (gate_depth.btotal.get(i) == false) {
            // reset the b vector
            gate_depth.btotal = Dynamic_bitset{n_};
            // set to true the locations of the last gate
            gate_depth.btotal.set(i);
            ++gate_depth.total;
        }

        return *this;
    }

    /**
     * \brief Controlled-NOT
     *
     * \param ctrl Control bit index
     * \param target Target bit index
     * \return Reference to the current instance
     */
    Bit_circuit& CNOT(idx ctrl, idx target) {
        v_[index_(target)] ^= (1 & (v_[index_(ctrl)] >> offset_(ctrl)))
                              << offset_(target);
        ++gate_count.CNOT;
        ++gate_count.total;

        // compute the depth
        if (gate_count.CNOT == 1)
            gate_depth.CNOT = 1;
        // apply the gate
        gate_depth.bCNOT.flip(ctrl).flip(target);
        // check whether gates overlap
        if (gate_depth.bCNOT.get(ctrl) == false ||
            gate_depth.bCNOT.get(target) == false) {
            // reset the b vector
            gate_depth.bCNOT = Dynamic_bitset{n_};
            // set to true the locations of the last gate
            gate_depth.bCNOT.set(ctrl).set(target);
            ++gate_depth.CNOT;
        }

        // compute the total depth
        if (gate_count.total == 1)
            gate_depth.total = 1;
        // apply the gate
        gate_depth.btotal.flip(ctrl).flip(target);
        // check whether gates overlap
        if (gate_depth.btotal.get(ctrl) == false ||
            gate_depth.btotal.get(target) == false) {
            // reset the b vector
            gate_depth.btotal = Dynamic_bitset{n_};
            // set to true the locations of the last gate
            gate_depth.btotal.set(ctrl).set(target);
            ++gate_depth.total;
        }

        return *this;
    }

    /**
     * \brief Toffoli gate
     *
     * \param i Control first bit index
     * \param j Control second bit index
     * \param k Target bit index
     * \return Reference to the current instance
     */
    Bit_circuit& TOF(idx i, idx j, idx k) {
        v_[index_(k)] ^= ((1 & (v_[index_(j)] >> offset_(j))) &
                          (1 & (v_[index_(i)] >> offset_(i))))
                         << offset_(k);
        ++gate_count.TOF;
        ++gate_count.total;

        // compute the depth
        if (gate_count.TOF == 1)
            gate_depth.TOF = 1;
        // apply the gate
        gate_depth.bTOF.flip(i).flip(j).flip(k);
        // check whether gates overlap
        if (gate_depth.bTOF.get(i) == false ||
            gate_depth.bTOF.get(j) == false ||
            gate_depth.bTOF.get(k) == false) {
            // reset the b vector
            gate_depth.bTOF = Dynamic_bitset{n_};
            // set to true the locations of the last gate
            gate_depth.bTOF.set(i).set(j).set(k);
            ++gate_depth.TOF;
        }

        // compute the total depth
        if (gate_count.total == 1)
            gate_depth.total = 1;
        // apply the gate
        gate_depth.btotal.flip(i).flip(j).flip(k);
        // check whether gates overlap
        if (gate_depth.btotal.get(i) == false ||
            gate_depth.btotal.get(j) == false ||
            gate_depth.btotal.get(k) == false) {
            // reset the b vector
            gate_depth.btotal = Dynamic_bitset{n_};
            // set to true the locations of the last gate
            gate_depth.btotal.set(i).set(j).set(k);
            ++gate_depth.total;
        }

        return *this;
    }

    /**
     * \brief Swap bits
     *
     * \param i Bit index
     * \param j Bit index
     * \return Reference to the current instance
     */
    Bit_circuit& SWAP(idx i, idx j) {
        if (this->get(i) != this->get(j)) {
            this->X(i);
            this->X(j);
        }
        ++gate_count.SWAP;
        ++gate_count.total;

        // compute the depth
        if (gate_count.SWAP == 1)
            gate_depth.SWAP = 1;
        // apply the gate
        gate_depth.bSWAP.flip(i).flip(j);
        // check whether gates overlap
        if (gate_depth.bSWAP.get(i) == false ||
            gate_depth.bSWAP.get(j) == false) {
            // reset the b vector
            gate_depth.bSWAP = Dynamic_bitset{n_};
            // set to true the locations of the last gate
            gate_depth.bSWAP.set(i).set(j);
            ++gate_depth.SWAP;
        }

        // compute the total depth
        if (gate_count.total == 1)
            gate_depth.total = 1;
        // apply the gate
        gate_depth.btotal.flip(i).flip(j);
        // check whether gates overlap
        if (gate_depth.btotal.get(i) == false ||
            gate_depth.btotal.get(j) == false) {
            // reset the b vector
            gate_depth.btotal = Dynamic_bitset{n_};
            // set to true the locations of the last gate
            gate_depth.btotal.set(i).set(j);
            ++gate_depth.total;
        }

        return *this;
    }

    /**
     * \brief Fredkin gate (Controlled-SWAP)
     *
     * \param i Control bit index
     * \param j Target first bit index
     * \param k Target second bit index
     * \return Reference to the current instance
     */
    Bit_circuit& FRED(idx i, idx j, idx k) {
        if (this->get(i)) {
            this->SWAP(j, k);
        }
        ++gate_count.FRED;
        ++gate_count.total;

        // compute the depth
        if (gate_count.FRED == 1)
            gate_depth.FRED = 1;
        // apply the gate
        gate_depth.bFRED.flip(i).flip(j).flip(k);
        // check whether gates overlap
        if (gate_depth.bFRED.get(i) == false ||
            gate_depth.bFRED.get(j) == false ||
            gate_depth.bFRED.get(k) == false) {
            // reset the b vector
            gate_depth.bFRED = Dynamic_bitset{n_};
            // set to true the locations of the last gate
            gate_depth.bFRED.set(i).set(j).set(k);
            ++gate_depth.FRED;
        }

        // compute the total depth
        if (gate_count.total == 1)
            gate_depth.total = 1;
        // apply the gate
        gate_depth.btotal.flip(i).flip(j).flip(k);
        // check whether gates overlap
        if (gate_depth.btotal.get(i) == false ||
            gate_depth.btotal.get(j) == false ||
            gate_depth.btotal.get(k) == false) {
            // reset the b vector
            gate_depth.btotal = Dynamic_bitset{n_};
            // set to true the locations of the last gate
            gate_depth.btotal.set(i).set(j).set(k);
            ++gate_depth.total;
        }

        return *this;
    }

    /**
     * \brief Reset the circuit all zero, clear all gates
     *
     * \return Reference to the current instance
     */
    Bit_circuit& reset() noexcept {
        gate_count.reset();
        gate_depth.reset();
        Dynamic_bitset::reset();

        return *this;
    }
}; /* class Bit_circuit */

} /* namespace qpp */

#endif /* CLASSES_REVERSIBLE_H */
