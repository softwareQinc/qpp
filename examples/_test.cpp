// Used for testing, do not use it as an example
#include <iostream>

#include "qpp.h"
#include "experimental/experimental.h"

using namespace qpp;

template <class T>
class QCirc {
    idx d_;                      // dimension, qubit circuit (d = 2) by default
    idx nq_, nc_;                // quantum bits, classical bits/dits
    ket psi_;                    // current state of the circuit
    std::vector<bool> measured_; // set elements to one if the corresponding
                                 // qubit had been measured
    std::vector<idx> results_;   // get_results of the measurements
    std::vector<idx> bits_;      // store the values of the classical bits/dits

  protected:
    std::vector<idx> update_subsys_(const std::vector<idx>& subsys) {
        std::vector<idx> result = subsys;
        idx subsys_size = subsys.size();
        for (idx i = 0; i < subsys_size; ++i) {
            for (idx m = 0; m < subsys[i]; ++m) {
                if (measured_[m]) { // if the qubit m was measured
                    --result[i];
                }
            }
        }
        std::sort(std::begin(result), std::end(result), std::less<idx>{});
        return result;
    }

  public:
    QCirc(idx nq, idx nc, idx d = 2)
        : d_{d}, nq_{nq}, nc_{nc}, psi_{st.zero(nq_, d_)},
          measured_(nq_, false), results_(nq_, -1), bits_(nc_, 0) {}

    // destructive measurement
    void measure(std::vector<idx> subsys) {
        // sort subsys in decreasing order
        idx subsys_size = subsys.size();
        for (idx i = 0; i < subsys_size; ++i) {
            if (measured_[subsys[i]])
                throw exception::CustomException(
                    "qpp::QCirc::measure()", "Subsystem was measured before");
        }
        // update subsystem labels
        std::vector<idx> subsys_updated = update_subsys_(subsys);
        for (idx i = 0; i < subsys_size; ++i) {
            measured_[subsys[i]] = true;
        }
        auto m = measure_seq(psi_, subsys_updated, d_);
        auto result = std::get<0>(m); // measurement result
        for (idx i = 0; i < subsys_size; ++i)
            results_[subsys[i]] = result[i];
        // update psi_
        psi_ = std::get<2>(m);
        // std::cout << disp(psi_) << '\n';
    }

    // destructive measurement of all remaining qubits
    void measure_all() {
        std::vector<idx> subsys;
        for (idx i = 0; i < nq_; ++i) {
            if (!measured_[i])
                subsys.push_back(i);
        }
        // update subsystem labels
        std::vector<idx> subsys_updated = update_subsys_(subsys);
        if (subsys.size() != 0)
            this->measure(subsys);
        else
            throw exception::CustomException("qpp::QCirc::measure_all()",
                                             "All qubits were measured before");
    }

    void apply(const cmat& gate, const std::vector<idx>& subsys) {
        psi_ = qpp::apply(psi_, gate, update_subsys_(subsys), d_);
    }

    void applyCTRL(const cmat& gate, const std::vector<idx>& ctrl,
                   const std::vector<idx>& target) {
        psi_ = qpp::applyCTRL(psi_, gate, ctrl, target, d_);
    }

    void apply_all(const cmat& gate) {
        for (idx i = 0; i < nq_; ++i)
            if (!measured_[i]) {
                psi_ = qpp::apply(psi_, gate, update_subsys_({i}), d_);
            }
    }

    // resets the circuit
    void reset() {
        psi_ = st.zero(nq_);
        measured_ = std::vector<bool>(nq_, false);
        results_ = std::vector<idx>(nq_, -1);
        bits_ = std::vector<idx>(nc_, 0);
    }

    /* getters */
    // local dimension
    idx d() const noexcept { return d_; }

    // total number of qubits, regardless of being measured or not
    idx get_nq() const noexcept { return nq_; }

    // total number of classical bits/dits
    idx get_nc() const noexcept { return nc_; }

    // total number of qubits AND classical bits/dits
    idx get_size() const noexcept { return nq_ + nc_; }

    // total number of measured qubits
    idx get_num_measured_qubits() const noexcept {
        return std::count(std::begin(measured_), std::end(measured_), true);
    }

    // total number of non-measured qubits
    idx get_num_active_qubits() const noexcept {
        return this->get_nq() - this->get_num_measured_qubits();
    }

    // returns the up-to-date quantum state
    ket get_psi() const { return psi_; }

    // returns the results of the measurements
    std::vector<idx> get_results() const { return results_; }

    // returns the current result of the measurement as an integer,
    // ignores non-measured qubits
    idx get_results_as_N() const {
        std::vector<idx> tmp;
        for (idx i = 0; i < nq_; ++i) {
            if (measured_[i])
                tmp.push_back(results_[i]);
        }

        return multiidx2n(tmp, std::vector<idx>(tmp.size(), d_));
    }

    // returns the classical bits/dits array
    std::vector<idx>& bits() noexcept { return bits_; }
    /* end getters */
};

int main() {
    QCirc<int> qc(10, 10);
    qc.apply_all(gt.H);
    qc.measure({3, 1, 7});
    std::cout << qc.get_num_measured_qubits() << '\n';
    std::cout << qc.get_num_active_qubits() << '\n';
    qc.measure({2, 4, 5, 6, 0});
    std::cout << disp(qc.get_results(), " ") << '\n';
    qc.apply_all(gt.H);
    qc.apply_all(gt.X);
    qc.measure_all();
    // qc.measure_all();
    std::cout << disp(qc.get_results(), "") << '\n';
    std::cout << qc.get_results_as_N() << '\n';
    std::cout << qc.get_nq() << '\n';
    std::cout << qc.get_num_measured_qubits() << '\n';
    std::cout << qc.get_num_active_qubits() << "\n\n";
    std::cout << "\n\n";

    QCirc<int> qc1(2, 0);
    qc1.apply(gt.H, {0});
    qc1.applyCTRL(gt.X, {0}, {1});
    qc1.apply(gt.CNOT, {0, 1});
    qc1.measure({0});

    std::cout << qc1.get_num_measured_qubits() << '\n';
    std::cout << qc1.get_num_active_qubits() << '\n';
    std::cout << disp(qc1.get_results(), " ") << '\n';
    std::cout << disp(qc1.get_psi()) << "\n\n";

    qc1.reset();
    std::cout << qc1.get_num_measured_qubits() << '\n';
    std::cout << qc1.get_num_active_qubits() << '\n';
    std::cout << disp(qc1.get_results(), " ") << '\n';
    std::cout << disp(qc1.get_psi()) << "\n\n";

    QCirc<int> qc2(2, 10, 3);
    qc2.apply(gt.Xd(3), {0});
    qc2.measure({1});
    std::cout << disp(qc2.get_psi()) << "\n\n";
    std::cout << qc2.get_size() << " " << qc2.d() << "\n\n";
}
