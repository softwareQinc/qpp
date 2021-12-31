/*
 * This file is part of qasmtools.
 *
 * Copyright (c) 2019 - 2022 softwareQ Inc. All rights reserved.
 *
 * MIT License
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
 * \file qasmtools/ast/decl.hpp
 * \brief OpenQASM declarations
 */

#pragma once

#include "stmt.hpp"

#include <list>

namespace qasmtools {
namespace ast {

#if USE_OPENQASM2_SPECS
static const std::set<std::string_view> qelib_defs{
    "u3", "u2",   "u1",  "cx",  "id",  "u0",  "x",  "y",  "z",
    "h",  "s",    "sdg", "t",   "tdg", "rx",  "ry", "rz", "cz",
    "cy", "swap", "ch",  "ccx", "crz", "cu1", "cu3"};
#else
/**
 * \brief Qiskit definitions include r and cswap gates, see
 * qasmtools/parser/preprocessor.hpp
 */
static const std::set<std::string_view> qelib_defs{
    "u3", "u2", "u1",   "cx", "id",  "u0",    "x",   "y",   "z",
    "h",  "s",  "sdg",  "t",  "tdg", "r",     "rx",  "ry",  "rz",
    "cz", "cy", "swap", "ch", "ccx", "cswap", "crz", "cu1", "cu3"};
#endif

/**
 * \brief Tests whether identifier is part of the standard OpenQASM qelib or not
 *
 * \param id Identifier
 * \return True if \a id is part of the standard OpenQASM qelib, false otherwise
 */
inline bool is_std_qelib(const std::string& id) {
    return qelib_defs.find(id) != qelib_defs.end();
}

/**
 * \class qasmtools::ast::Decl
 * \brief Base class for OpenQASM declarations
 *
 * Declarations are attribute classes as they can occur in different
 * statement contexts. To avoid diamond inheritance, any derived declaration
 * should also inherit from a statement class
 */
class Decl {
  protected:
    symbol id_; ///< the name of the declaration

  public:
    Decl(symbol id) : id_(id) {}
    virtual ~Decl() = default;

    /**
     * \brief Return the name being declared
     *
     * \return Constant reference to the identifier
     */
    const symbol& id() { return id_; }
};

/**
 * \class qasmtools::ast::GateDecl
 * \brief Class for gate declarations
 * \see qasmtools::ast::Stmt
 * \see qasmtools::ast::Decl
 */
class GateDecl final : public Stmt, public Decl {
    bool opaque_;                  ///< whether the declaration is opaque
    std::vector<symbol> c_params_; ///< classical parameters
    std::vector<symbol> q_params_; ///< quantum parameters
    std::list<ptr<Gate>> body_;    ///< gate body

  public:
    /**
     * \brief Constructs a gate declaration
     *
     * \param pos The source position
     * \param id The gate identifier
     * \param c_params List of classical parameters
     * \param q_params List of quantum parameters
     * \param body List of gate statements
     */
    GateDecl(parser::Position pos, symbol id, bool opaque,
             std::vector<symbol> c_params, std::vector<symbol> q_params,
             std::list<ptr<Gate>>&& body)
        : Stmt(pos), Decl(id), opaque_(opaque), c_params_(c_params),
          q_params_(q_params), body_(std::move(body)) {}

    /**
     * \brief Protected heap-allocated construction
     */
    static ptr<GateDecl> create(parser::Position pos, symbol id, bool opaque,
                                std::vector<symbol> c_params,
                                std::vector<symbol> q_params,
                                std::list<ptr<Gate>>&& body) {
        return std::make_unique<GateDecl>(pos, id, opaque, c_params, q_params,
                                          std::move(body));
    }

    /**
     * \brief Whether the declaration is opaque
     *
     * \return true is the declaration is opaque
     */
    bool is_opaque() { return opaque_; }

    /**
     * \brief Get the classical parameter list
     *
     * \return Reference to the list of classical parameter names
     */
    std::vector<symbol>& c_params() { return c_params_; }

    /**
     * \brief Get the quantum parameter list
     *
     * \return Reference to the list of quantum parameter names
     */
    std::vector<symbol>& q_params() { return q_params_; }

    /**
     * \brief Get the gate body
     *
     * \return Reference to the body of the gate as a list of gate statements
     */
    std::list<ptr<Gate>>& body() { return body_; }

    /**
     * \brief Apply a function to each statement of the gate
     *
     * \param f A void function taking a reference to a Gate
     */
    void foreach_stmt(std::function<void(Gate&)> f) {
        for (auto it = body_.begin(); it != body_.end(); it++)
            f(**it);
    }

    /**
     * \brief Get an iterator to the beginning of the body
     *
     * \return std::list iterator
     */
    std::list<ptr<Gate>>::iterator begin() { return body_.begin(); }

    /**
     * \brief Get an iterator to the end of the body
     *
     * \return std::list iterator
     */
    std::list<ptr<Gate>>::iterator end() { return body_.end(); }

    void accept(Visitor& visitor) override { visitor.visit(*this); }
    std::ostream& pretty_print(std::ostream& os,
                               bool suppress_std) const override {
        if (suppress_std && is_std_qelib(id_))
            return os;

        os << (opaque_ ? "opaque " : "gate ") << id_;
        if (c_params_.size() > 0) {
            os << "(";
            for (auto it = c_params_.begin(); it != c_params_.end(); it++) {
                os << (it == c_params_.begin() ? "" : ",") << *it;
            }
            os << ")";
        }
        os << " ";
        for (auto it = q_params_.begin(); it != q_params_.end(); it++) {
            os << (it == q_params_.begin() ? "" : ",") << *it;
        }
        if (opaque_) {
            os << ";\n";
        } else {
            os << " {\n";
            for (auto it = body_.begin(); it != body_.end(); it++) {
                os << "\t" << **it;
            }
            os << "}\n";
        }
        return os;
    }

  protected:
    GateDecl* clone() const override {
        std::list<ptr<Gate>> tmp;
        for (auto it = body_.begin(); it != body_.end(); it++) {
            tmp.emplace_back(object::clone(**it));
        }
        return new GateDecl(pos_, id_, opaque_, c_params_, q_params_,
                            std::move(tmp));
    }
};

/**
 * \class qasmtools::ast::OracleDecl
 * \brief Class for oracle declarations
 * \see qasmtools::ast::Decl
 */
class OracleDecl final : public Stmt, public Decl {
    std::vector<symbol> params_; ///< quantum parameters
    symbol fname_;               ///< filename of external declaration

  public:
    /**
     * \brief Constructs an oracle declaration
     *
     * \param pos The source position
     * \param id The gate identifier
     * \param params List of quantum parameters
     * \param fname Filename defining the classical logic
     */
    OracleDecl(parser::Position pos, symbol id, std::vector<symbol> params,
               symbol fname)
        : Stmt(pos), Decl(id), params_(params), fname_(fname) {}

    /**
     * \brief Protected heap-allocated construction
     */
    static ptr<OracleDecl> create(parser::Position pos, symbol id,
                                  std::vector<symbol> params, symbol fname) {
        return std::make_unique<OracleDecl>(pos, id, params, fname);
    }

    /**
     * \brief Get the parameter list
     *
     * \return Reference to the list of quantum parameter names
     */
    std::vector<symbol>& params() { return params_; }

    /**
     * \brief Get the filename
     *
     * \return Constant reference to the filename
     */
    const symbol& fname() { return fname_; }

    void accept(Visitor& visitor) override { visitor.visit(*this); }
    std::ostream& pretty_print(std::ostream& os, bool) const override {
        os << "oracle " << id_ << " ";
        for (auto it = params_.begin(); it != params_.end(); it++) {
            os << (it == params_.begin() ? "" : ",") << *it;
        }
        os << " { \"" << fname_ << "\" }\n";
        return os;
    }

  protected:
    OracleDecl* clone() const override {
        return new OracleDecl(pos_, id_, params_, fname_);
    }
};

/**
 * \class qasmtools::ast::RegisterDecl
 * \brief Class for register declarations
 * \see qasmtools::ast::Decl
 */
class RegisterDecl final : public Stmt, public Decl {
    bool quantum_; ///< whether the register is quantum
    int size_;     ///< the size of the register

  public:
    /**
     * \brief Constructs a register declaration
     *
     * \param pos The source position
     * \param id The register identifier
     * \param quantum whether the register is a quantum register
     * \param size the size of the register
     */
    RegisterDecl(parser::Position pos, symbol id, bool quantum, int size)
        : Stmt(pos), Decl(id), quantum_(quantum), size_(size) {}

    /**
     * \brief Protected heap-allocated construction
     */
    static ptr<RegisterDecl> create(parser::Position pos, symbol id,
                                    bool quantum, int size) {
        return std::make_unique<RegisterDecl>(pos, id, quantum, size);
    }

    /**
     * \brief Whether the register is quantum or classical
     *
     * \return true if the register is quantum
     */
    bool is_quantum() { return quantum_; }

    /**
     * \brief Get the size of the register
     *
     * \return The size of the register
     */
    int size() { return size_; }

    void accept(Visitor& visitor) override { visitor.visit(*this); }
    std::ostream& pretty_print(std::ostream& os, bool) const override {
        os << (quantum_ ? "qreg " : "creg ") << id_ << "[" << size_ << "];\n";
        return os;
    }

  protected:
    RegisterDecl* clone() const override {
        return new RegisterDecl(pos_, id_, quantum_, size_);
    }
};

/**
 * \class qasmtools::ast::AncillaDecl
 * \brief Class for local register declarations
 * \see qasmtools::ast::Decl
 */
class AncillaDecl final : public Gate, public Decl {
    bool dirty_; ///< whether the register can be dirty
    int size_;   ///< the size of the register

  public:
    /**
     * \brief Constructs a register declaration
     *
     * \param pos The source position
     * \param id The register identifier
     * \param dirty Whether the register is dirty
     * \param size The size of the register
     */
    AncillaDecl(parser::Position pos, symbol id, bool dirty, int size)
        : Gate(pos), Decl(id), dirty_(dirty), size_(size) {}

    /**
     * \brief Protected heap-allocated construction
     */
    static ptr<AncillaDecl> create(parser::Position pos, symbol id, bool dirty,
                                   int size) {
        return std::make_unique<AncillaDecl>(pos, id, dirty, size);
    }

    /**
     * \brief Whether the register is dirty
     *
     * \return true if the register is dirty
     */
    bool is_dirty() { return dirty_; }

    /**
     * \brief Get the size of the register
     *
     * \return The size of the register
     */
    int size() { return size_; }

    void accept(Visitor& visitor) override { visitor.visit(*this); }
    std::ostream& pretty_print(std::ostream& os, bool) const override {
        if (dirty_)
            os << "dirty ";
        os << "ancilla " << id_ << "[" << size_ << "];\n";
        return os;
    }

  protected:
    AncillaDecl* clone() const override {
        return new AncillaDecl(pos_, id_, dirty_, size_);
    }
};

} // namespace ast
} // namespace qasmtools
