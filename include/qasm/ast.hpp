/*
 * This file is part of Quantum++.
 *
 * Copyright (c) 2013 - 2021 softwareQ Inc. All rights reserved.
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
 * \file qasm/ast.hpp
 * \brief openQASM 2.0 AST
 */

#ifndef QASM_AST_HPP_
#define QASM_AST_HPP_

namespace qpp {
namespace qasm {

using ident = std::string;

/**
 * \brief Enum of binary operators
 */
enum class BinaryOp {
    Plus,
    Minus,
    Times,
    Divide,
    Pow,
};

/**
 * \brief Enum of unary operators
 */
enum class UnaryOp {
    Neg,
    Sin,
    Cos,
    Tan,
    Ln,
    Sqrt,
    Exp,
};

static std::unordered_map<ident,
                          std::function<cmat(const std::vector<double>&)>>
    known_matrices{
        ///< generators for various gate constants
        {"cx",
         [](const std::vector<double>&) { return Gates::get_instance().CNOT; }},
        {"id",
         [](const std::vector<double>&) { return Gates::get_instance().Id2; }},
        {"x",
         [](const std::vector<double>&) { return Gates::get_instance().X; }},
        {"y",
         [](const std::vector<double>&) { return Gates::get_instance().Y; }},
        {"z",
         [](const std::vector<double>&) { return Gates::get_instance().Z; }},
        {"h",
         [](const std::vector<double>&) { return Gates::get_instance().H; }},
        {"s",
         [](const std::vector<double>&) { return Gates::get_instance().S; }},
        {"sdg",
         [](const std::vector<double>&) {
             return Gates::get_instance().S.adjoint();
         }},
        {"t",
         [](const std::vector<double>&) { return Gates::get_instance().T; }},
        {"tdg",
         [](const std::vector<double>&) {
             return Gates::get_instance().T.adjoint();
         }},
        {"rx",
         [](const std::vector<double>& args) {
             assert(!args.empty());
             return Gates::get_instance().RX(args[0]);
         }},
        {"rz",
         [](const std::vector<double>& args) {
             assert(!args.empty());
             // note the discrepancy; Qiskit defines it as diag(1, e^{i\phi})
             // we comply to the Qiskit definition (and not the OPENQASM
             // specs); see https://github.com/softwareQinc/qpp/issues/70
             return (std::exp(1_i * args[0] / 2.0) *
                     Gates::get_instance().RZ(args[0]))
                 .eval();
         }},
        {"ry",
         [](const std::vector<double>& args) {
             assert(!args.empty());
             return Gates::get_instance().RY(args[0]);
         }},
        {"cz",
         [](const std::vector<double>&) { return Gates::get_instance().CZ; }},
        {"cy",
         [](const std::vector<double>&) {
             cmat mat{cmat::Identity(4, 4)};
             mat.block(2, 2, 2, 2) = Gates::get_instance().Y;
             return mat;
         }},
        {"swap",
         [](const std::vector<double>&) { return Gates::get_instance().SWAP; }},
        {"ch",
         [](const std::vector<double>&) {
             cmat mat{cmat::Identity(4, 4)};
             mat.block(2, 2, 2, 2) = Gates::get_instance().H;
             return mat;
         }},
        {"ccx",
         [](const std::vector<double>&) { return Gates::get_instance().TOF; }},
        {"crz", [](const std::vector<double>& args) {
             assert(!args.empty());
             cmat mat{cmat::Identity(4, 4)};
             // note the discrepancy; Qiskit defines it as
             // CTRL-diag(1, e^{i\phi}) we comply to the Qiskit definition (and
             // not the OPENQASM specs); see
             // https://github.com/softwareQinc/qpp/issues/99 and
             // https://github.com/softwareQinc/qpp/issues/70
             mat.block(2, 2, 2, 2) = std::exp(1_i * args[0] / 2.0) *
                                     Gates::get_instance().RZ(args[0]);
             return mat;
         }}};

/**
 * \class qpp::qasm::Value
 * \brief Interface class for openQASM values during evaluation
 *
 * Allows environments to contain different types of mappings depending on
 * the type of variable. See deriving classes
 */
class Value {
  public:
    virtual ~Value() = default;
};

/**
 * \class qpp::qasm::Context
 * \brief QCircuit translation context
 *
 * Stores all information relevant to execution (i.e., generation of a
 * QCircuit). Each classical and quantum bit in the AST is mapped to a unique
 * classical or quantum index corresponding to a QCircuit qubit or bit.
 * Measurements are always non-destructive.
 */
class Context {
    QCircuit* circuit_; ///< pointer to the accumulating circuit

    // Hack for MSVC
    using hash_ident_uptr = std::unordered_map<ident, std::unique_ptr<Value>>;
    struct Environment {
        Environment() noexcept : val_(){};
        Environment(Environment&& rhs) noexcept : val_(std::move(rhs.val_)) {}
        hash_ident_uptr val_;
    };

    std::vector<Environment> env_{};     ///< environment stack
    std::list<idx> qubit_pool_{};        ///< pool of unassigned physical qubits
    idx max_bit_ = static_cast<idx>(-1); ///< largest classical bit index
    idx max_qubit_ = static_cast<idx>(-1); ///< largest (virtual) qubit index

    // For controlled contexts
    std::vector<idx> cctrls_{}; ///< classical controls in the current context
    std::vector<idx> cctrl_shift_{}; ///< control shift

  public:
    /**
     * \brief Constructs a translation context with a given target QCircuit
     *
     * \param qc Pointer to a QCircuit object
     */
    explicit Context(QCircuit* qc) : circuit_(qc) {}

    /**
     * \brief Disables copy constructor
     */
    Context(const Context&) = delete;

    /**
     * \brief Disables assignment
     */
    Context& operator=(const Context&) = delete;

    /**
     * \brief The accumulating QCircuit
     *
     * \return Pointer to the QCircuit object
     */
    QCircuit* get_circuit() { return circuit_; }

    /*-------------- (Qu)bit management & allocation ---------------*/

    /**
     * \brief Allocate a classical bit
     *
     * \return Index of a fresh classical bit
     */
    idx alloc_bit() {
        idx ret = ++max_bit_;

        // check if circuit has enough classical bits
        if (max_bit_ >= circuit_->get_nc()) {
            circuit_->add_dit(1, ret);
        }

        return ret;
    }

    /**
     * \brief Allocate a virtual qubit
     *
     * \return Index of a fresh qubit
     */
    idx alloc_qubit() {
        idx ret = ++max_qubit_;

        // check if circuit has enough classical bits
        if (max_qubit_ >= circuit_->get_nq()) {
            circuit_->add_qudit(1, ret);
        }

        return ret;
    }

    /*------------------ Environment management --------------*/
    /**
     * \brief Enter a new scope
     */
    void enter_scope() { env_.emplace_back(); }

    /**
     * \brief Exit current scope
     */
    void exit_scope() { env_.pop_back(); }

    /**
     * \brief Lookup an identifier in the environment
     *
     * Traverses all scopes from inner-most to outer-most and returns
     * the first mapping
     *
     * \param id Const reference to an identifier
     * \param loc Const reference to the location in the source file
     * \return Pointer to a value
     */
    Value* lookup(const ident& id, const Location& loc) {
        for (auto table = env_.rbegin(); table != env_.rend(); table++) {
            auto it = table->val_.find(id);
            if (it != table->val_.end())
                return it->second.get();
        }

        std::stringstream context;
        context << loc << ": Undeclared identifier " << id;
        throw exception::Undeclared("qpp::qasm::Context::lookup()",
                                    context.str());
    }

    /**
     * \brief Set the value of an identifier in the current scope
     * \note Assumes ownership of the value
     *
     * \param id Const reference to an identifier
     * \param val Unique pointer to a value
     */
    void set(const ident& id, std::unique_ptr<Value> val) {
        if (env_.empty())
            env_.emplace_back();
        env_.back().val_[id] = std::move(val);
    }

    /*------------------- Classical controls ----------------*/
    /**
     * \brief Enter a classically controlled context
     * \note QASM syntax does not allow nested classically controlled contexts
     *
     * \param cctrls Const reference to a list of classical controls
     * \param shift Const reference to a list of classical control values
     */
    void set_ccontrols(const std::vector<idx>& cctrls,
                       const std::vector<idx>& shift) {
        cctrls_ = cctrls;
        cctrl_shift_ = shift;
    }

    /**
     * \brief Exit a classically controlled context
     */
    void clear_ccontrols() {
        cctrls_.clear();
        cctrl_shift_.clear();
    }

    /**
     * \brief Whether the current context is classically controlled
     *
     * \return True if and only if the current context is classically controlled
     */
    bool ccontrolled() { return !cctrls_.empty(); }

    /**
     * \brief The current classical controls
     *
     * \return Const reference to the classical controls
     */
    const std::vector<idx>& get_cctrls() { return cctrls_; }

    /**
     * \brief The current classical control values
     *
     * \return Const reference to the classical control values
     */
    const std::vector<idx>& get_shift() { return cctrl_shift_; }
};

/*------------------- Virtual base classes -------------------*/

/**
 * \class qpp::qasm::Statement
 * \brief Base class for openQASM statements
 */
class Statement : public IDisplay {
  protected:
    Location loc_; ///< source location of the statement

  public:
    /**
     * \brief Constructs a statement
     *
     * \param loc Source location
     */
    explicit Statement(const Location& loc) : loc_(loc) {}

    /**
     * \brief Disable copying
     */
    Statement(const Statement&) = delete;

    /**
     * \brief Must be overridden by all derived classes
     *
     * Implements the QCircuit translation by appending the statement
     * to the QCircuit in the given context
     *
     * \param ctx Translation context
     */
    virtual void evaluate(Context& ctx) const = 0;

    /**
     * \brief Formats and prints the statement
     *
     * \param os Output stream passed by reference
     * \param prefix Prefix for the statement
     */
    virtual void pretty_print(std::ostream& os,
                              const std::string& prefix) const {
        os << "(" << loc_ << "):" << prefix << *this;
    }
};
using StatementPtr = std::unique_ptr<Statement>;

/**
 * \class qpp::qasm::Gate
 * \brief Base class for openQASM gates
 */
class Gate : public Statement {
  public:
    /**
     * \brief Constructs a gate statement
     *
     * \param loc Source location
     */
    explicit Gate(const Location& loc) : Statement(loc) {}

    /**
     * \brief Disable copying
     */
    Gate(const Gate&) = delete;
};
using GatePtr = std::unique_ptr<Gate>;

/**
 * \class qpp::qasm::Decl
 * \brief Base class for openQASM declarations
 * \see qpp::qasm::Statement
 */
class Decl : public Statement {
  protected:
    ident id_; ///< declared identifier

  public:
    /**
     * \brief Constructs a declaration statement
     *
     * \param loc Source location
     * \param id Declared identifier
     */
    Decl(const Location& loc, ident id) : Statement(loc), id_(std::move(id)) {}

    /**
     * \brief Disable copying
     */
    Decl(const Decl&) = delete;
};

/**
 * \class qpp::qasm::Expr
 * \brief Base class for openQASM expressions
 */
class Expr : public IDisplay {
  protected:
    Location loc_; ///<  source location of the expression

  public:
    /**
     * \brief Constructs an expression
     *
     * \param loc Source location
     */
    explicit Expr(const Location& loc) : loc_(loc) {}

    /**
     * \brief Disable copying
     */
    Expr(const Expr&) = delete;

    /**
     * \brief Must be overridden by all derived classes
     *
     * Evaluates the expression to a floating point value
     * in the given context
     *
     * \param ctx Translation context
     */
    virtual double evaluate(Context& ctx) const = 0;
};
using ExprPtr = std::unique_ptr<Expr>;

/*-------------------- Non-virtual leaf classes ------------------*/
/**
 * \class qpp::qasm::QASM
 * \brief QASM program class
 * \brief QASM program class
 */
class QASM : public IDisplay {
    int bits_ = -1;                  ///< number of bits
    int qubits_ = -1;                ///< number of qubits
    std::vector<StatementPtr> body_; ///< program body

  public:
    /**
     * \brief Constructs a QASM program
     *
     * \param bits Number of bits
     * \param qubits Number of qubits
     * \param body Rvalue reference to a list of statements
     */
    QASM(int bits, int qubits, std::vector<StatementPtr>&& body)
        : bits_(bits), qubits_(qubits), body_(std::move(body)) {}

    /**
     * \brief Constructs a QCircuit object from the QASM AST
     *
     * \return Unique pointer to a QCircuit object
     */
    std::unique_ptr<QCircuit> to_QCircuit() {
        auto ret = std::unique_ptr<QCircuit>(new QCircuit(qubits_, bits_));
        Context ctx(ret.get());
        for (auto& it : body_) {
            it->evaluate(ctx);
        }

        return ret;
    }

    /**
     * \brief qpp::IDisplay::display() override
     */
    std::ostream& display(std::ostream& os) const override {
        os << "OPENQASM 2.0;\ninclude \"qelib1.inc\";\n";
        for (const auto& it : body_) {
            os << *it;
        }

        return os;
    }

    /**
     * \brief Formats and prints the AST with location annotations
     *
     * \param os Output stream passed by reference
     */
    void pretty_print(std::ostream& os) const {
        os << "(Included header): OPENQASM 2.0;\n(Included header): include "
              "\"qelib1.inc\";\n";
        for (const auto& it : body_) {
            it->pretty_print(os, "\t");
        }
    }
};

/**
 * \class qpp::qasm::Circuit
 * \brief QASM circuit values
 * \see qpp::qasm::Value
 */
class Circuit : public Value {
  public:
    const std::vector<ident>& c_params_; ///< the classical parameters
    const std::vector<ident>& q_params_; ///< the quantum parameters
    const std::vector<GatePtr>& body_;   ///< the circuit body

    Circuit(const std::vector<ident>& c_params,
            const std::vector<ident>& q_params,
            const std::vector<GatePtr>& body)
        : c_params_(c_params), q_params_(q_params), body_(body) {}
};

/**
 * \class qpp::qasm::Register
 * \brief QASM register values
 * \see qpp::qasm::Value
 */
class Register : public Value {
  public:
    bool quantum_;             ///< whether the register is a qreg or creg
    std::vector<idx> indices_; ///< the (virtual) register indices

    Register(bool quantum, std::vector<idx> indices)
        : quantum_(quantum), indices_(std::move(indices)) {}
};

/**
 * \class qpp::qasm::Qubit
 * \brief QASM qubit values
 * \see qpp::qasm::Value
 */
class Qubit : public Value {
  public:
    idx index_; ///< index of the qubit

    explicit Qubit(idx index) : index_(index) {}
};

/**
 * \class qpp::qasm::Number
 * \brief QASM number values
 * \see qpp::qasm::Value
 */
class Number : public Value {
  public:
    double value_; ///< the floating point number

    explicit Number(double value) : value_(value) {}
};

/**
 * \class qpp::qasm::Varinfo
 * \brief Class for variable accesses
 */
class Varinfo : public IDisplay {
    Location loc_; ///< location in source
    ident id_;     ///< the accessed variable
    int offset_; ///< the offset of the register access, or -1 if not derefenced

  public:
    /**
     * \brief Constructs a variable access
     *
     * \param loc The source location
     * \param id The identifier accessed
     * \param offset The register offset (optional)
     */
    Varinfo(const Location& loc, ident id, int offset = -1)
        : loc_(loc), id_(std::move(id)), offset_(offset) {}

    /**
     * \brief Try to interpret the access as a classical register in a given
     * context
     * \note Throws qpp:exception::SemanticError if ill-typed
     *
     * \param ctx The translation context
     * \return List of bit indices corresponding to the variable access
     */
    std::vector<idx> as_creg(Context& ctx) const {
        auto reg = dynamic_cast<const Register*>(ctx.lookup(id_, loc_));

        if (reg == nullptr || reg->quantum_) {
            std::stringstream context;
            context << loc_ << ": Identifier " << id_
                    << " does not refer to a bit register";
            throw exception::SemanticError("qpp::qasm::Varinfo::as_creg()",
                                           context.str());
        }

        if (offset_ == -1) {
            // creg
            return std::vector<idx>{reg->indices_};
        } else {
            // creg deref

            // check register size
            if ((idx) offset_ >= reg->indices_.size()) {
                std::stringstream context;
                context << loc_ << ": Index out of bounds";
                throw exception::SemanticError("qpp::qasm::Varinfo::as_qreg()",
                                               context.str());
            }

            return std::vector<idx>{reg->indices_[offset_]};
        }
    }

    /**
     * \brief Try to interpret the access as a quantum register in a given
     * context
     * \note Throws qpp:exception::SemanticError if ill-typed
     *
     * \param ctx The translation context
     * \return List of virtual qubit indices corresponding to the variable
     * access
     */
    std::vector<idx> as_qreg(Context& ctx) const {
        if (offset_ == -1) {
            // qubit or qreg
            auto tmp = ctx.lookup(id_, loc_);

            if (auto qubit = dynamic_cast<const Qubit*>(tmp)) {
                return std::vector<idx>{qubit->index_};
            } else {
                auto reg = dynamic_cast<const Register*>(tmp);
                // check register type
                if (reg == nullptr || !reg->quantum_) {
                    std::stringstream context;
                    context << loc_ << ": Identifier " << id_
                            << " does not refer to a qubit register";
                    throw exception::SemanticError(
                        "qpp::qasm::Varinfo::as_qreg()", context.str());
                }

                return std::vector<idx>{reg->indices_};
            }
        } else {
            // qreg deref
            auto reg = dynamic_cast<const Register*>(ctx.lookup(id_, loc_));

            // check register type
            if (reg == nullptr || !reg->quantum_) {
                std::stringstream context;
                context << loc_ << ": Identifier " << id_
                        << " does not refer to a qubit register";
                throw exception::SemanticError("qpp::qasm::Varinfo::as_qreg()",
                                               context.str());
            }

            // check register size
            if ((idx) offset_ >= reg->indices_.size()) {
                std::stringstream context;
                context << loc_ << ": Index out of bounds";
                throw exception::SemanticError("qpp::qasm::Varinfo::as_qreg()",
                                               context.str());
            }

            return std::vector<idx>{reg->indices_[offset_]};
        }
    }

    /**
     * \brief qpp::IDisplay::display() override
     */
    std::ostream& display(std::ostream& os) const override {
        os << id_;
        if (offset_ != -1)
            os << "[" << offset_ << "]";

        return os;
    }
};

/**
 * \class qpp::qasm::GateDecl
 * \brief Class for gate declarations
 * \see qpp::qasm::Decl
 */
class GateDecl final : public Decl {
    bool opaque_;                 ///< whether the declaration is opaque
    std::vector<ident> c_params_; ///< classical parameters
    std::vector<ident> q_params_; ///< quantum parameters
    std::vector<GatePtr> body_;   ///< gate body

  public:
    /**
     * \brief Constructs a gate declaration
     *
     * \param loc The source location
     * \param id The gate identifier
     * \param opaque Opaque
     * \param c_params List of classical parameters
     * \param q_params List of quantum parameters
     * \param body List of gate statements
     */
    GateDecl(const Location& loc, ident id, bool opaque,
             std::vector<ident> c_params, std::vector<ident> q_params,
             std::vector<GatePtr> body)
        : Decl(loc, std::move(id)), opaque_(opaque),
          c_params_(std::move(c_params)), q_params_(std::move(q_params)),
          body_(std::move(body)) {}

    /**
     * \brief qpp::qasm::Statement::evaluate() override
     */
    void evaluate(Context& ctx) const override {
        ctx.set(id_, std::unique_ptr<Value>(
                         new Circuit(c_params_, q_params_, body_)));
    }

    /**
     * \brief qpp::IDisplay::display() override
     */
    std::ostream& display(std::ostream& os) const override {
        os << (opaque_ ? "opaque " : "gate ") << id_;
        if (!c_params_.empty()) {
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
            os << "{\n";
            for (const auto& it : body_) {
                os << "\t" << *it;
            }
            os << "}\n";
        }
        return os;
    }

    /**
     * \brief qpp::qasm::Statement::pretty_print() override
     * \note Necessary since textual representation is multi-line
     */
    void pretty_print(std::ostream& os,
                      const std::string& prefix) const override {
        os << "(" << loc_ << "):" << prefix;
        os << (opaque_ ? "opaque " : "gate ") << id_;
        if (!c_params_.empty()) {
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
            auto nprefix = prefix + "\t";
            os << "{\n";
            for (const auto& it : body_) {
                it->pretty_print(os, nprefix);
            }
            os << "(" << loc_ << "):" << prefix << "}\n";
        }
    }
};

/**
 * \class qpp::qasm::RegisterDecl
 * \brief Class for register declarations
 * \see qpp::qasm::Decl
 */
class RegisterDecl final : public Decl {
    bool quantum_; ///< whether the register is quantum
    idx length_;   ///< the length of the register

  public:
    /**
     * \brief Constructs a register declaration
     *
     * \param loc The source location
     * \param id The register identifier
     * \param quantum whether the register is a quantum register
     * \param length the length of the register
     */
    RegisterDecl(const Location& loc, ident id, bool quantum, idx length)
        : Decl(loc, std::move(id)), quantum_(quantum), length_(length) {}

    /**
     * \brief qpp::qasm::Statement::evaluate() override
     */
    void evaluate(Context& ctx) const override {
        std::vector<idx> indices(length_);
        for (idx i = 0; i < length_; i++) {
            indices[i] = quantum_ ? ctx.alloc_qubit() : ctx.alloc_bit();
        }

        ctx.set(id_, std::unique_ptr<Value>(new Register(quantum_, indices)));
    }

    /**
     * \brief qpp::IDisplay::display() override
     */
    std::ostream& display(std::ostream& os) const override {
        os << (quantum_ ? "qreg " : "creg ") << id_ << "[" << length_ << "];\n";
        return os;
    }
};

/**
 * \class qpp::qasm::MeasurementStatement
 * \brief Class for measurement statements
 * \see qpp::qasm::Statement
 */
class MeasureStatement final : public Statement {
    Varinfo q_arg_; ///< the quantum bit|register
    Varinfo c_arg_; ///< the classical bit|register

  public:
    /**
     * \brief Constructs a measurement statement
     *
     * \param loc The source location
     * \param q_arg The quantum argument
     * \param c_arg The classical argument
     */
    MeasureStatement(const Location& loc, Varinfo q_arg, Varinfo c_arg)
        : Statement(loc), q_arg_(std::move(q_arg)), c_arg_(std::move(c_arg)) {}

    /**
     * \brief qpp::qasm::Statement::evaluate() override
     */
    void evaluate(Context& ctx) const override {
        auto q_args = q_arg_.as_qreg(ctx);
        auto c_args = c_arg_.as_creg(ctx);
        auto circuit = ctx.get_circuit();

        // check register lengths
        if (q_args.size() != c_args.size()) {
            std::stringstream context;
            context << loc_ << ": Registers have different lengths";
            throw exception::ParseError(
                "qpp::qasm::MeasureStatement::evaluate()", context.str());
        }

        // apply measurements non-desctructively
        for (idx i = 0; i < q_args.size(); i++) {
            circuit->measureZ(q_args[i], c_args[i], false);
        }
    }

    /**
     * \brief qpp::IDisplay::display() override
     */
    std::ostream& display(std::ostream& os) const override {
        os << "measure " << q_arg_ << " -> " << c_arg_ << ";\n";
        return os;
    }
};

/**
 * \class qpp::qasm::ResetStatement
 * \brief Class for reset statements
 * \see qpp::qasm::Statement
 */
class ResetStatement final : public Statement {
    Varinfo arg_; ///< the qbit|qreg

  public:
    /**
     * \brief Constructs a reset statement
     *
     * \param loc The source location
     * \param arg The argument
     */
    ResetStatement(const Location& loc, Varinfo arg)
        : Statement(loc), arg_(std::move(arg)) {}

    /**
     * \brief qpp::qasm::Statement::evaluate() override
     */
    void evaluate(Context& ctx) const override {
        auto q_args = arg_.as_qreg(ctx);
        auto circuit = ctx.get_circuit();

        circuit->reset(q_args);
    }

    /**
     * \brief qpp::IDisplay::display() override
     */
    std::ostream& display(std::ostream& os) const override {
        os << "reset " << arg_ << ";\n";
        return os;
    }
};

/**
 * \class qpp::qasm::IfStatement
 * \brief Class for if statements
 * \see qpp::qasm::Statement
 */
class IfStatement final : public Statement {
    ident id_;          ///< classical register name
    int value_;         ///< value to check against
    StatementPtr then_; ///< statement to be executed if true

  public:
    /**
     * \brief Constructs an if statement
     *
     * \param loc The source location
     * \param id The classical register name
     * \param value The value to check against
     * \param then Unique pointer to the statement to be executed
     */
    IfStatement(const Location& loc, ident id, int value, StatementPtr then)
        : Statement(loc), id_(std::move(id)), value_(value),
          then_(std::move(then)) {}

    /**
     * \brief qpp::qasm::Statement::evaluate() override
     */
    void evaluate(Context& ctx) const override {
        auto creg = dynamic_cast<const Register*>(ctx.lookup(id_, loc_));

        // check register type
        if (creg == nullptr || creg->quantum_) {
            std::stringstream context;
            context << loc_ << ": Identifier " << id_
                    << " does not refer to a classic register";
            throw exception::SemanticError("qpp::qasm::IfStatement::evaluate()",
                                           context.str());
        }

        // create the shift
        int tmp = value_;
        std::vector<idx> shift(creg->indices_.size(), 0);
        for (idx i = 0; i < creg->indices_.size(); i++) {
            if (tmp % 2 == 0)
                shift[i] = 1;
            tmp >>= 1;
        }

        // apply controls to then branch
        ctx.set_ccontrols(creg->indices_, shift);
        then_->evaluate(ctx);
        ctx.clear_ccontrols();
    }

    /**
     * \brief qpp::IDisplay::display() override
     */
    std::ostream& display(std::ostream& os) const override {
        os << "if (" << id_ << "==" << value_ << ") " << *then_;
        return os;
    }
};

/**
 * \class qpp::qasm::UGate
 * \brief Class for U gates
 * \see qpp::qasm::Gate
 */
class UGate final : public Gate {
    ExprPtr theta_;  ///< theta angle
    ExprPtr phi_;    ///< phi angle
    ExprPtr lambda_; ///< lambda angle

    Varinfo arg_; ///< quantum bit|register
  public:
    /**
     * \brief Constructs a single qubit U gate
     *
     * \param loc The source location
     * \param theta Unique pointer to an angle expression
     * \param phi Unique pointer to an angle expression
     * \param lambda Unique pointer to an angle expression
     * \param arg The quantum bit or register
     */
    UGate(const Location& loc, ExprPtr theta, ExprPtr phi, ExprPtr lambda,
          Varinfo arg)
        : Gate(loc), theta_(std::move(theta)), phi_(std::move(phi)),
          lambda_(std::move(lambda)), arg_(std::move(arg)) {}

    /**
     * \brief qpp::qasm::Statement::evaluate() override
     */
    void evaluate(Context& ctx) const override {
        auto theta = theta_->evaluate(ctx);
        auto phi = phi_->evaluate(ctx);
        auto lambda = lambda_->evaluate(ctx);
        auto args = arg_.as_qreg(ctx);
        auto circuit = ctx.get_circuit();

        // generate the matrix
        cmat u{cmat::Zero(2, 2)};

        // standard QASM spec, as defined in
        // https://arxiv.org/pdf/1707.03429.pdf
        // u << std::cos(theta / 2) * std::exp(-1_i * (phi + lambda) / 2.0),
        //     -(std::sin(theta / 2)) * std::exp(-1_i * (phi - lambda) / 2.0),
        //     std::sin(theta / 2) * std::exp(1_i * (phi - lambda) / 2.0),
        //     std::cos(theta / 2) * std::exp(1_i * (phi + lambda) / 2.0);

        // Qiskit spec, as defined in
        // https://github.com/Qiskit/qiskit-terra/tree/master/qiskit/circuit/library/standard_gates
        // We use these definitions, see
        // https://github.com/vsoftco/qpp/issues/65 for the reasons why.
        u << std::cos(theta / 2),
            -(std::sin(theta / 2)) * std::exp(1_i * lambda),
            std::sin(theta / 2) * std::exp(1_i * phi),
            std::cos(theta / 2) * std::exp(1_i * (phi + lambda));

        // apply the gate
        for (auto i : args) {
            if (ctx.ccontrolled()) {
                circuit->cCTRL(u, ctx.get_cctrls(), i, ctx.get_shift(),
                               "Controlled-U");
            } else {
                circuit->gate(u, i, "U");
            }
        }
    }

    /**
     * \brief qpp::IDisplay::display() override
     */
    std::ostream& display(std::ostream& os) const override {
        os << "U(" << *theta_ << "," << *phi_ << "," << *lambda_ << ") " << arg_
           << ";\n";
        return os;
    }
};

/**
 * \class qpp::qasm::CNOTGate
 * \brief Class for CX gates
 * \see qpp::qasm::Gate
 */
class CNOTGate final : public Gate {
    Varinfo ctrl_; ///< control qubit|qreg
    Varinfo tgt_;  ///< target qubit|qreg

  public:
    /**
     * \brief Constructs a CNOT gate
     *
     * \param loc The source location
     * \param ctrl The control bit or register
     * \param tgt The target bit or register
     */
    CNOTGate(const Location& loc, Varinfo ctrl, Varinfo tgt)
        : Gate(loc), ctrl_(std::move(ctrl)), tgt_(std::move(tgt)) {}

    /**
     * \brief qpp::qasm::Statement::evaluate() override
     */
    void evaluate(Context& ctx) const override {
        auto ctrls = ctrl_.as_qreg(ctx);
        auto tgts = tgt_.as_qreg(ctx);
        auto circuit = ctx.get_circuit();

        // different combinations of controls and targets
        if (ctrls.size() == 1 && tgts.size() == 1) {
            if (ctx.ccontrolled()) {
                std::vector<idx> tmp{ctrls[0], tgts[0]};
                circuit->cCTRL_joint(Gates::get_instance().CNOT,
                                     ctx.get_cctrls(), tmp, ctx.get_shift(),
                                     "CX");
            } else {
                circuit->gate(Gates::get_instance().CNOT, ctrls[0], tgts[0],
                              "CX");
            }
        } else if (ctrls.size() > 1 && tgts.size() == 1) {
            for (idx ctrl : ctrls) {
                if (ctx.ccontrolled()) {
                    std::vector<idx> tmp{ctrl, tgts[0]};
                    circuit->cCTRL_joint(Gates::get_instance().CNOT,
                                         ctx.get_cctrls(), tmp, ctx.get_shift(),
                                         "CX");
                } else {
                    circuit->gate(Gates::get_instance().CNOT, ctrl, tgts[0],
                                  "CX");
                }
            }
        } else if (ctrls.size() == 1 && tgts.size() > 1) {
            for (idx tgt : tgts) {
                if (ctx.ccontrolled()) {
                    std::vector<idx> tmp{ctrls[0], tgt};
                    circuit->cCTRL_joint(Gates::get_instance().CNOT,
                                         ctx.get_cctrls(), tmp, ctx.get_shift(),
                                         "CX");
                } else {
                    circuit->gate(Gates::get_instance().CNOT, ctrls[0], tgt,
                                  "CX");
                }
            }
        } else if (ctrls.size() == tgts.size()) {
            for (idx i = 0; i < ctrls.size(); i++) {
                if (ctx.ccontrolled()) {
                    std::vector<idx> tmp{ctrls[i], tgts[i]};
                    circuit->cCTRL_joint(Gates::get_instance().CNOT,
                                         ctx.get_cctrls(), tmp, ctx.get_shift(),
                                         "CX");
                } else {
                    circuit->gate(Gates::get_instance().CNOT, ctrls[i], tgts[i],
                                  "CX");
                }
            }
        } else {
            std::stringstream context;
            context << loc_ << ": Registers have different lengths";
            throw exception::SemanticError("qpp::qasm::CNOTGate::evaluate()",
                                           context.str());
        }
    }

    /**
     * \brief qpp::IDisplay::display() override
     */
    std::ostream& display(std::ostream& os) const override {
        os << "CX " << ctrl_ << tgt_ << ";\n";
        return os;
    }
};

/**
 * \class qpp::qasm::BarrierGate
 * \brief Class for barrier gates
 * \see qpp::qasm::Gate
 */
class BarrierGate final : public Gate {
    std::vector<Varinfo> args_; ///< list of quantum bits|registers

  public:
    /**
     * \brief Constructs a barrier gate
     *
     * \param loc The source location
     * \param args The quantum bits or registers
     */
    BarrierGate(const Location& loc, std::vector<Varinfo> args)
        : Gate(loc), args_(std::move(args)) {}

    /**
     * \brief qpp::qasm::Statement::evaluate() override
     */
    void evaluate(Context& ctx) const override {
        // just check validity and do nothing
        idx mapping_size = 1;

        for (auto& arg : args_) {
            auto tmp = arg.as_qreg(ctx);
            if (tmp.size() > 1) {
                if (mapping_size == 1) {
                    mapping_size = tmp.size();
                } else if (mapping_size != tmp.size()) {
                    std::stringstream context;
                    context << loc_ << ": Registers have different lengths";
                    throw exception::SemanticError(
                        "qpp::qasm::BarrierGate::evaluate()", context.str());
                }
            }
        }
    }

    /**
     * \brief qpp::IDisplay::display() override
     */
    std::ostream& display(std::ostream& os) const override {
        os << "barrier ";
        for (auto it = args_.begin(); it != args_.end(); it++) {
            os << (it == args_.begin() ? "" : ",") << *it;
        }
        os << ";\n";
        return os;
    }
};

/**
 * \class qpp::qasm::DeclaredGate
 * \brief Class for declared gate applications
 * \see qpp::qasm::Gate
 */
class DeclaredGate final : public Gate {
    ident id_;                    ///< gate identifier
    std::vector<ExprPtr> c_args_; ///< list of classical arguments
    std::vector<Varinfo> q_args_; ///< list of quantum arguments

  public:
    /**
     * \brief Constructs a gate application
     *
     * \param loc The source location
     * \param id The gate name
     * \param c_args Rvalue reference to a list of the classical arguments
     * \param q_args List of the quantum arguments
     */
    DeclaredGate(const Location& loc, ident id, std::vector<ExprPtr>&& c_args,
                 std::vector<Varinfo> q_args)
        : Gate(loc), id_(std::move(id)), c_args_(std::move(c_args)),
          q_args_(std::move(q_args)) {}

    /**
     * \brief qpp::qasm::Statement::evaluate() override
     */
    void evaluate(Context& ctx) const override {
        auto gate = dynamic_cast<const Circuit*>(ctx.lookup(id_, loc_));

        // check gate type
        if (gate == nullptr) {
            std::stringstream context;
            context << loc_ << ": Identifier " << id_
                    << " does not refer to a gate declaration";
            throw exception::SemanticError(
                "qpp::qasm::DeclaredGate::evaluate()", context.str());
        }

        // check argument lengths
        std::stringstream context;
        if (c_args_.size() != gate->c_params_.size()) {
            context << loc_ << ": " << id_ << " expects "
                    << gate->c_params_.size();
            context << " classic arguments, got " << c_args_.size();
            throw exception::SemanticError(
                "qpp::qasm::DeclaredGate::evaluate()", context.str());
        } else if (q_args_.size() != gate->q_params_.size()) {
            context << loc_ << ": " << id_ << " expects "
                    << gate->q_params_.size();
            context << " quantum arguments, got " << q_args_.size();
            throw exception::SemanticError(
                "qpp::qasm::DeclaredGate::evaluate()", context.str());
        }

        // evaluate arguments
        std::vector<double> c_args(c_args_.size());
        std::vector<std::vector<idx>> q_args(q_args_.size());
        for (idx i = 0; i < c_args_.size(); i++)
            c_args[i] = c_args_[i]->evaluate(ctx);
        for (idx i = 0; i < q_args_.size(); i++)
            q_args[i] = q_args_[i].as_qreg(ctx);

        // map gate across registers
        idx mapping_size = 1;
        std::vector<bool> mapped(q_args.size(), false);
        for (idx i = 0; i < q_args.size(); i++) {
            if (q_args[i].size() > 1) {
                mapped[i] = true;
                if (mapping_size == 1) {
                    mapping_size = q_args[i].size();
                } else if (mapping_size != q_args[i].size()) {
                    std::stringstream inner_context;
                    inner_context << loc_
                                  << ": Registers have different lengths";
                    throw exception::SemanticError(
                        "qpp::qasm::DeclaredGate::evaluate()",
                        inner_context.str());
                }
            }
        }

        auto it = known_matrices.find(id_);
        if (it != known_matrices.end()) {
            // apply the known matrix directly
            auto circuit = ctx.get_circuit();
            auto mat = it->second(c_args);

            // map the gate accross registers
            for (idx j = 0; j < mapping_size; j++) {
                // map virtual qubits to physical qubits
                std::vector<idx> mapped_args(q_args.size());
                for (idx i = 0; i < q_args.size(); i++) {
                    mapped_args[i] = mapped[i] ? q_args[i][j] : q_args[i][0];
                }

                // apply (possibly classical controlled) gate
                if (ctx.ccontrolled()) {
                    circuit->cCTRL_joint(mat, ctx.get_cctrls(), mapped_args,
                                         ctx.get_shift(), id_);
                } else {
                    circuit->gate_joint(mat, mapped_args, id_);
                }
            }
        } else {
            // push classical arguments onto a new scope
            ctx.enter_scope();
            for (idx i = 0; i < c_args.size(); i++) {
                ctx.set(gate->c_params_[i],
                        std::unique_ptr<Value>(new Number(c_args[i])));
            }

            // map the gate
            for (idx j = 0; j < mapping_size; j++) {
                ctx.enter_scope();
                for (idx i = 0; i < q_args.size(); i++) {
                    ctx.set(gate->q_params_[i],
                            std::unique_ptr<Value>(new Qubit(
                                mapped[i] ? q_args[i][j] : q_args[i][0])));
                }

                // evaluate the gate
                for (auto&& elem : gate->body_) {
                    elem->evaluate(ctx);
                }

                ctx.exit_scope();
            }
            ctx.exit_scope();
        }
    }

    /**
     * \brief qpp::IDisplay::display() override
     */
    std::ostream& display(std::ostream& os) const override {
        os << id_;
        if (!c_args_.empty()) {
            os << "(";
            for (auto it = c_args_.begin(); it != c_args_.end(); it++) {
                os << (it == c_args_.begin() ? "" : ",") << **it;
            }
            os << ")";
        }
        os << " ";
        for (auto it = q_args_.begin(); it != q_args_.end(); it++) {
            os << (it == q_args_.begin() ? "" : ",") << *it;
        }
        os << ";\n";
        return os;
    }
};

/**
 * \class qpp::qasm::RealExpr
 * \brief Class for floating point literal expressions
 * \see qpp::qasm::Expr
 */
class RealExpr final : public Expr {
    double value_; ///< the floating point value

  public:
    /**
     * \brief Constructs a real literal expression
     *
     * \param loc The source location
     * \param value The floating point value
     */
    RealExpr(const Location& loc, double value) : Expr(loc), value_(value) {}

    /**
     * \brief qpp::qasm::Expr::evaluate() override
     */
    double evaluate(Context& ctx) const override {
        (void) ctx;
        return value_;
    }

    /**
     * \brief qpp::IDisplay::display() override
     */
    std::ostream& display(std::ostream& os) const override {
        os << value_;
        return os;
    }
};

/**
 * \class qpp::qasm::IntExpr
 * \brief Class for integer literal expressions
 * \see qpp::qasm::Expr
 */
class IntExpr final : public Expr {
    int value_; ///< the integer value

  public:
    /**
     * \brief Constructs an integer literal expression
     *
     * \param loc The source location
     * \param value The integer value
     */
    IntExpr(const Location& loc, int value) : Expr(loc), value_(value) {}

    /**
     * \brief qpp::qasm::Expr::evaluate() override
     */
    double evaluate(Context& ctx) const override {
        (void) ctx;
        return (double) value_;
    }

    /**
     * \brief qpp::IDisplay::display() override
     */
    std::ostream& display(std::ostream& os) const override {
        os << value_;
        return os;
    }
};

/**
 * \class qpp::qasm::PiExpr
 * \brief Class for pi constants
 * \see qpp::qasm::Expr
 */
class PiExpr final : public Expr {
  public:
    /**
     * \brief Constructs a pi constant expression
     *
     * \param loc The source location
     */
    explicit PiExpr(const Location& loc) : Expr(loc) {}

    /**
     * \brief qpp::qasm::Expr::evaluate() override
     */
    double evaluate(Context& ctx) const override {
        (void) ctx;
        return qpp::pi;
    }

    /**
     * \brief qpp::IDisplay::display() override
     */
    std::ostream& display(std::ostream& os) const override {
        os << "pi";
        return os;
    }
};

/**
 * \class qpp::qasm::VarExpr
 * \brief Class for variable expressions
 * \see qpp::qasm::Expr
 */
class VarExpr final : public Expr {
    ident value_; ///< the variable identifier

  public:
    /**
     * \brief Constructs a variable expression
     *
     * \param loc The source location
     * \param value The variable identifier
     */
    VarExpr(const Location& loc, ident value)
        : Expr(loc), value_(std::move(value)) {}

    /**
     * \brief qpp::qasm::Expr::evaluate() override
     */
    double evaluate(Context& ctx) const override {
        auto val = dynamic_cast<const Number*>(ctx.lookup(value_, loc_));

        if (val == nullptr) {
            std::stringstream context;
            context << loc_ << ": Identifier " << value_
                    << " does not refer to a classical parameter";
            throw exception::ParseError("qpp::qasm::VarExpr::evaluate()",
                                        context.str());
        }

        return val->value_;
    }

    /**
     * \brief qpp::IDisplay::display() override
     */
    std::ostream& display(std::ostream& os) const override {
        os << value_;
        return os;
    }
};

/**
 * \class qpp::qasm::BExpr
 * \brief Class for binary operator expressions
 * \see qpp::qasm::Expr
 */
class BExpr final : public Expr {
    ExprPtr lexp_; ///< the left sub-expression
    BinaryOp bop_; ///< the binary operator
    ExprPtr rexp_; ///< the right sub-expression

  public:
    /**
     * \brief Constructs a binary expression
     *
     * \param loc The source location
     * \param lexp Unique pointer to the left sub-expression
     * \param bop The binary operator
     * \param rexp Unique pointer to the right sub-expression
     */
    BExpr(const Location& loc, ExprPtr lexp, BinaryOp bop, ExprPtr rexp)
        : Expr(loc), lexp_(std::move(lexp)), bop_(bop), rexp_(std::move(rexp)) {
    }

    /**
     * \brief qpp::qasm::Expr::evaluate() override
     */
    double evaluate(Context& ctx) const override {
        double lval = lexp_->evaluate(ctx);
        double rval = rexp_->evaluate(ctx);

        switch (bop_) {
            case BinaryOp::Plus:
                return lval + rval;
            case BinaryOp::Minus:
                return lval - rval;
            case BinaryOp::Times:
                return lval * rval;
            case BinaryOp::Divide:
                return lval / rval;
            case BinaryOp::Pow:
                return pow(lval, rval);
            default:
                return 0;
        }
    }

    /**
     * \brief qpp::IDisplay::display() override
     */
    std::ostream& display(std::ostream& os) const override {
        os << "(" << *lexp_;
        switch (bop_) {
            case BinaryOp::Plus:
                os << "+";
                break;
            case BinaryOp::Minus:
                os << "-";
                break;
            case BinaryOp::Times:
                os << "*";
                break;
            case BinaryOp::Divide:
                os << "/";
                break;
            case BinaryOp::Pow:
                os << "^";
                break;
        }
        os << *rexp_ << ")";
        return os;
    }
};

/**
 * \class qpp::qasm::UExpr
 * \brief Class for unary operator expressions
 * \see qpp::qasm::Expr
 */
class UExpr final : public Expr {
    UnaryOp uop_; ///< the unary operator
    ExprPtr exp_; ///< the sub-expression

  public:
    /**
     * \brief Constructs a unary expression
     *
     * \param loc The source location
     * \param uop The unary operator
     * \param exp Unique pointer to the sub-expression
     */
    UExpr(const Location& loc, UnaryOp uop, ExprPtr exp)
        : Expr(loc), uop_(uop), exp_(std::move(exp)) {}

    /**
     * \brief qpp::qasm::Expr::evaluate() override
     */
    double evaluate(Context& ctx) const override {
        double val = exp_->evaluate(ctx);

        switch (uop_) {
            case UnaryOp::Neg:
                return -val;
            case UnaryOp::Sin:
                return std::sin(val);
            case UnaryOp::Cos:
                return std::cos(val);
            case UnaryOp::Tan:
                return std::tan(val);
            case UnaryOp::Ln:
                return std::log(val);
            case UnaryOp::Sqrt:
                return std::sqrt(val);
            case UnaryOp::Exp:
                return std::exp(val);
            default:
                return 0;
        }
    }

    /**
     * \brief qpp::IDisplay::display() override
     */
    std::ostream& display(std::ostream& os) const override {
        switch (uop_) {
            case UnaryOp::Neg:
                os << "-" << *exp_;
                break;
            case UnaryOp::Sin:
                os << "sin(" << *exp_ << ")";
                break;
            case UnaryOp::Cos:
                os << "cos(" << *exp_ << ")";
                break;
            case UnaryOp::Tan:
                os << "tan(" << *exp_ << ")";
                break;
            case UnaryOp::Ln:
                os << "ln(" << *exp_ << ")";
                break;
            case UnaryOp::Sqrt:
                os << "sqrt(" << *exp_ << ")";
                break;
            case UnaryOp::Exp:
                os << "exp(" << *exp_ << ")";
                break;
        }
        return os;
    }
};

} /* namespace qasm */
} /* namespace qpp */

#endif /* QASM_AST_HPP_ */
