/*
 * This file is part of Quantum++.
 *
 * MIT License
 *
 * Copyright (c) 2013 - 2020 Vlad Gheorghiu.
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

#ifndef CLASSES_REVERSIBLE_HPP_
#define CLASSES_REVERSIBLE_HPP_

namespace qpp {
/**
 * \class qpp::Dynamic_bitset
 * \brief Dynamic bitset class, allows the specification of the number of bits
 * at runtime
 * \note The interface mimics std::bitset<>
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
    inline static idx index_(idx pos) {
        return pos / (sizeof(value_type) * CHAR_BIT);
    }

    /**
     * \brief Offset of the \a pos bit in the storage space relative to its
     * index
     *
     * \param pos Bit location
     * \return Offset of the \a pos bit in the storage space relative to its
     * index
     */
    inline static idx offset_(idx pos) {
        return pos % (sizeof(value_type) * CHAR_BIT);
    }

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
    ~Dynamic_bitset() override = default;

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
        idx bitset_size = size();
        for (idx i = 0; i < bitset_size; ++i) {
            if (get(i))
                ++result;
        }

        return result;
    }

    /**
     * \brief The value of the bit at position \a pos
     *
     * \param pos Position in the bitset
     * \return Value of the bit at position \a pos
     */
    bool get(idx pos) const noexcept {
        return static_cast<value_type>(1) & (v_[index_(pos)] >> offset_(pos));
    }

    /**
     * \brief Checks whether none of the bits are set
     *
     * \return True if none of the bits are set
     */
    bool none() const noexcept {
        bool result = true;
        idx bitset_storage_size = storage_size();
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
        idx bitset_storage_size = storage_size();
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
    bool any() const noexcept { return !(none()); }

    /* setters */
    /**
     * \brief Sets the bit at position \a pos
     *
     * \param pos Position in the bitset
     * \param value Bit value
     * \return Reference to the current instance
     */
    Dynamic_bitset& set(idx pos, bool value = true) {
        value
            ? v_[index_(pos)] |= (static_cast<value_type>(1) << offset_(pos))
            : v_[index_(pos)] &= ~(static_cast<value_type>(1) << offset_(pos));

        //        v_[index_(pos)] &= ~(!value << offset_(pos));

        return *this;
    }

    /**
     * \brief Set all bits to true
     *
     * \return Reference to the current instance
     */
    Dynamic_bitset& set() noexcept {
        idx bitset_storage_size = storage_size();
        for (idx i = 0; i < bitset_storage_size; ++i) {
            v_[i] = ~static_cast<value_type>(0);
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

        set(pos, d(gen));

        return *this;
    }

    /**
     * \brief Sets all bits according to a Bernoulli(p) distribution
     *
     * \param p Probability
     * \return Reference to the current instance
     */
    Dynamic_bitset& rand(double p = 0.5) {
        idx bitset_size = size();
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
        v_[index_(pos)] &= ~(static_cast<value_type>(1) << offset_(pos));

        return *this;
    }

    /**
     * \brief Sets all bits to false
     *
     * \return Reference to the current instance
     */
    Dynamic_bitset& reset() noexcept {
        idx bitset_storage_size = storage_size();
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
        v_[index_(pos)] ^= static_cast<value_type>(1) << (offset_(pos));

        return *this;
    }

    /**
     * \brief Flips all bits
     *
     * \return Reference to the current instance
     */
    Dynamic_bitset& flip() noexcept {
        idx bitset_storage_size = storage_size();
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
        assert(size() == rhs.size());
        bool result = true;
        idx n = std::min(storage_size(), rhs.storage_size());
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
        idx bitset_size = size();
        for (idx i = 0; i < bitset_size; ++i) {
            if (get(i) != rhs.get(i))
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
     * \return Bitset as a string
     */
    template <class CharT = char, class Traits = std::char_traits<CharT>,
              class Allocator = std::allocator<CharT>>
    std::basic_string<CharT, Traits, Allocator>
    to_string(CharT zero = CharT('0'), CharT one = CharT('1')) const {
        std::basic_string<CharT, Traits, Allocator> result;
        idx bitset_size = size();
        result.resize(bitset_size);

        for (idx i = bitset_size; i-- > 0;) {
            if (!get(i)) {
                result[bitset_size - i - 1] = zero;
            } else {
                result[bitset_size - i - 1] = one;
            }
        }

        return result;
    }

  protected:
    /**
     * \brief qpp::IDisplay::display() override
     *
     * Writes to the output stream a textual representation of the bitset
     * (displays the bitset bit by bit_
     *
     * \param os Output stream passed by reference
     * \return Reference to the output stream
     */
    std::ostream& display(std::ostream& os) const override {
        idx bitset_size = size();
        for (idx i = bitset_size; i-- > 0;) {
            os << get(i);
        }

        return os;
    }
}; /* class Dynamic_bitset */

/**
 * \class qpp::Bit_circuit
 * \brief Classical reversible circuit simulator
 */
class Bit_circuit : public Dynamic_bitset, public IJSON {
    std::unordered_map<std::string, idx> count_{}; ///< gate counts
    std::unordered_map<std::string, idx> depth_{}; ///< gate depths
    Dynamic_bitset bNOT_, bCNOT_, bSWAP_, bTOF_, bFRED_,
        btotal_; ///< used for depth calculations

  public:
    /**
     * \brief Constructs a bit circuit instance
     *
     * \param n Number of classical bits
     */
    explicit Bit_circuit(idx n)
        : Dynamic_bitset{n}, bNOT_{n}, bCNOT_{n}, bSWAP_{n}, bTOF_{n},
          bFRED_{n}, btotal_{n} {}
    /**
     * \brief Conversion constructor, used to initialize a qpp::Bit_circuit with
     * a qpp::Dynamic_bitset
     *
     * \param dynamic_bitset Dynamic bitset
     */
    explicit Bit_circuit(const Dynamic_bitset& dynamic_bitset)
        : Dynamic_bitset{dynamic_bitset}, bNOT_{size()}, bCNOT_{size()},
          bSWAP_{size()}, bTOF_{size()}, bFRED_{size()}, btotal_{size()} {}

    /**
     * \brief Bit flip
     * \see qpp::Bit_circuit::NOT()
     *
     * \param i Bit position in the circuit
     * \return Reference to the current instance
     */
    Bit_circuit& X(idx i) {
        NOT(i);

        return *this;
    }

    /**
     * \brief Default virtual destructor
     */
    ~Bit_circuit() override = default;

    /**
     * \brief Bit flip
     * \see qpp::Bit_circuit::X()
     *
     * \param i Bit position in the circuit
     * \return Reference to the current instance
     */
    Bit_circuit& NOT(idx i) {
        flip(i);
        ++count_["NOT"];
        ++count_[__FILE__ "__total__"];

        // compute the depth
        if (count_["NOT"] == 1)
            depth_["NOT"] = 1;
        // apply the gate
        bNOT_.flip(i);
        // check whether gates overlap
        if (!bNOT_.get(i)) {
            // reset the b bitset
            bNOT_ = Dynamic_bitset{n_};
            // set to true the locations of the last gate
            bNOT_.set(i);
            ++depth_["NOT"];
        }

        // compute the total depth
        if (count_[__FILE__ "__total__"] == 1)
            depth_[__FILE__ "__total__"] = 1;
        // apply the gate
        btotal_.flip(i);
        // check whether gates overlap
        if (!btotal_.get(i)) {
            // reset the b bitset
            btotal_ = Dynamic_bitset{n_};
            // set to true the locations of the last gate
            btotal_.set(i);
            ++depth_[__FILE__ "__total__"];
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
        v_[index_(target)] ^=
            (static_cast<value_type>(1) & (v_[index_(ctrl)] >> offset_(ctrl)))
            << offset_(target);
        ++count_["CNOT"];
        ++count_[__FILE__ "__total__"];

        // compute the depth
        if (count_["CNOT"] == 1)
            depth_["CNOT"] = 1;
        // apply the gate
        bCNOT_.flip(ctrl).flip(target);
        // check whether gates overlap
        if (!bCNOT_.get(ctrl) || !bCNOT_.get(target)) {
            // reset the b bitset
            bCNOT_ = Dynamic_bitset{n_};
            // set to true the locations of the last gate
            bCNOT_.set(ctrl).set(target);
            ++depth_["CNOT"];
        }

        // compute the total depth
        if (count_[__FILE__ "__total__"] == 1)
            depth_[__FILE__ "__total__"] = 1;
        // apply the gate
        btotal_.flip(ctrl).flip(target);
        // check whether gates overlap
        if (!btotal_.get(ctrl) || !btotal_.get(target)) {
            // reset the b bitset
            btotal_ = Dynamic_bitset{n_};
            // set to true the locations of the last gate
            btotal_.set(ctrl).set(target);
            ++depth_[__FILE__ "__total__"];
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
        v_[index_(k)] ^=
            ((static_cast<value_type>(1) & (v_[index_(j)] >> offset_(j))) &
             (static_cast<value_type>(1) & (v_[index_(i)] >> offset_(i))))
            << offset_(k);
        ++count_["TOF"];
        ++count_[__FILE__ "__total__"];

        // compute the depth
        if (count_["TOF"] == 1)
            depth_["TOF"] = 1;
        // apply the gate
        bTOF_.flip(i).flip(j).flip(k);
        // check whether gates overlap
        if (!bTOF_.get(i) || !bTOF_.get(j) || !bTOF_.get(k)) {
            // reset the b bitset
            bTOF_ = Dynamic_bitset{n_};
            // set to true the locations of the last gate
            bTOF_.set(i).set(j).set(k);
            ++depth_["TOF"];
        }

        // compute the total depth
        if (count_[__FILE__ "__total__"] == 1)
            depth_[__FILE__ "__total__"] = 1;
        // apply the gate
        btotal_.flip(i).flip(j).flip(k);
        // check whether gates overlap
        if (!btotal_.get(i) || !btotal_.get(j) || !btotal_.get(k)) {
            // reset the b bitset
            btotal_ = Dynamic_bitset{n_};
            // set to true the locations of the last gate
            btotal_.set(i).set(j).set(k);
            ++depth_[__FILE__ "__total__"];
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
        if (get(i) != get(j)) {
            X(i);
            X(j);
        }
        ++count_["SWAP"];
        ++count_[__FILE__ "__total__"];

        // compute the depth
        if (count_["SWAP"] == 1)
            depth_["SWAP"] = 1;
        // apply the gate
        bSWAP_.flip(i).flip(j);
        // check whether gates overlap
        if (!bSWAP_.get(i) || !bSWAP_.get(j)) {
            // reset the b bitset
            bSWAP_ = Dynamic_bitset{n_};
            // set to true the locations of the last gate
            bSWAP_.set(i).set(j);
            ++depth_["SWAP"];
        }

        // compute the total depth
        if (count_[__FILE__ "__total__"] == 1)
            depth_[__FILE__ "__total__"] = 1;
        // apply the gate
        btotal_.flip(i).flip(j);
        // check whether gates overlap
        if (!btotal_.get(i) || !btotal_.get(j)) {
            // reset the b bitset
            btotal_ = Dynamic_bitset{n_};
            // set to true the locations of the last gate
            btotal_.set(i).set(j);
            ++depth_[__FILE__ "__total__"];
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
        if (get(i)) {
            SWAP(j, k);
        }
        ++count_["FRED"];
        ++count_[__FILE__ "__total__"];

        // compute the depth
        if (count_["FRED"] == 1)
            depth_["FRED"] = 1;
        // apply the gate
        bFRED_.flip(i).flip(j).flip(k);
        // check whether gates overlap
        if (!bFRED_.get(i) || !bFRED_.get(j) || !bFRED_.get(k)) {
            // reset the b bitset
            bFRED_ = Dynamic_bitset{n_};
            // set to true the locations of the last gate
            bFRED_.set(i).set(j).set(k);
            ++depth_["FRED"];
        }

        // compute the total depth
        if (count_[__FILE__ "__total__"] == 1)
            depth_[__FILE__ "__total__"] = 1;
        // apply the gate
        btotal_.flip(i).flip(j).flip(k);
        // check whether gates overlap
        if (!btotal_.get(i) || !btotal_.get(j) || !btotal_.get(k)) {
            // reset the b bitset
            btotal_ = Dynamic_bitset{n_};
            // set to true the locations of the last gate
            btotal_.set(i).set(j).set(k);
            ++depth_[__FILE__ "__total__"];
        }

        return *this;
    }

    /**
     * \brief Reset the circuit all zero, clear all gates
     *
     * \return Reference to the current instance
     */
    Bit_circuit& reset() noexcept {
        count_ = {};
        depth_ = {};
        Dynamic_bitset::reset();

        return *this;
    }

    // getters
    /**
     * \brief Bit circuit gate count
     *
     * \param name Gate name. Possible names are NOT (X), CNOT, SWAP, TOF, FRED.
     * \return Gate count
     */
    idx get_gate_count(const std::string& name) const {
        idx result;

        // EXCEPTION CHECKS

        try {
            if (name == "X")
                result = count_.at("NOT");
            else
                result = count_.at(name);
        } catch (...) {
            return 0;
        }
        // END EXCEPTION CHECKS

        return result;
    }

    /**
     * \brief Bit circuit total gate count
     *
     * \return Total gate count
     */
    idx get_gate_count() const {
        idx result;

        // EXCEPTION CHECKS

        try {
            result = count_.at(__FILE__ "__total__");
        } catch (...) {
            return 0;
        }
        // END EXCEPTION CHECKS

        return result;
    }

    /**
     * \brief Bit circuit gate depth
     *
     * \param name Gate name. Possible names are NOT (X), CNOT, SWAP, TOF, FRED.
     * \return Gate depth
     */
    idx get_gate_depth(const std::string& name) const {
        idx result;

        // EXCEPTION CHECKS

        try {
            if (name == "X")
                result = depth_.at("NOT");
            else
                result = depth_.at(name);
        } catch (...) {
            return 0;
        }
        // END EXCEPTION CHECKS

        return result;
    }

    /**
     * \brief Bit circuit total gate depth
     *
     * \return Total gate depth
     */
    idx get_gate_depth() const {
        idx result;

        // EXCEPTION CHECKS

        try {
            result = depth_.at(__FILE__ "__total__");
        } catch (...) {
            return 0;
        }
        // END EXCEPTION CHECKS

        return result;
    }
    // end getters

    /**
     * \brief qpp::IJSON::to_JSON() override
     *
     * Displays the bit circuit in JSON format
     *
     * \param enclosed_in_curly_brackets If true, encloses the result in curly
     * brackets
     * \return String containing the JSON representation of the state of the
     * engine
     */
    std::string to_JSON(bool enclosed_in_curly_brackets = true) const override {
        std::string result;

        if (enclosed_in_curly_brackets)
            result += "{";

        result += "\"n\" : " + std::to_string(n_);
        result +=
            ", \"total gate count\" : " + std::to_string(get_gate_count());
        result +=
            ", \"total gate depth\" : " + std::to_string(get_gate_depth());
        result += R"(, "bit state" : ")" + this->to_string() + '\"';
        result += ", \"Hamming weight\" : " + std::to_string(count());

        if (enclosed_in_curly_brackets)
            result += "}";

        return result;
    }

  private:
    /**
     * \brief qpp::IDisplay::display() override
     *
     * Writes to the output stream a textual representation of the bit circuit
     *
     * \param os Output stream passed by reference
     * \return Reference to the output stream
     */
    std::ostream& display(std::ostream& os) const override {
        os << "n = " << n_ << '\n';
        os << "total gate count: " << get_gate_count() << '\n';
        os << "total gate depth: " << get_gate_depth() << '\n';
        os << "bit state: " << this->to_string() << '\n';
        os << "Hamming weight: " << count();

        return os;
    }

}; /* class Bit_circuit */

} /* namespace qpp */

#endif /* CLASSES_REVERSIBLE_HPP_ */
