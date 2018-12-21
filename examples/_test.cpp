// Used for testing, do not use it as an example
#include <iostream>

#include "qpp.h"
#include "experimental/experimental.h"

using namespace qpp;

template <class T>
class QCirc {
    idx nq_, nc_;              // quantum bits, classical bits
    ket psi_;                  // current state of the circuit
    Dynamic_bitset measured_;  // set elements to one if the corresponding
                               // qubit had been measured
    std::vector<idx> results_; // results of the measurements
    std::vector<idx> update_subsys_(const std::vector<idx>& subsys) {
        std::vector<idx> result = subsys;
        idx subsys_size = subsys.size();
        for (idx i = 0; i < subsys_size; ++i) {
            for (idx m = 0; m < subsys[i]; ++m) {
                if (measured_.get(m)) { // if the qubit m was measured
                    --result[i];
                }
            }
        }
        std::sort(std::begin(result), std::end(result), std::less<idx>{});
        return result;
    }

  public:
    QCirc(idx nq, idx nc)
        : nq_{nq}, nc_{nc}, psi_{st.zero(nq_)}, measured_{nq_},
          results_(nq_, -1) {}

    // destructive measurement
    void measure(std::vector<idx> subsys) {
        // sort subsys in decreasing order
        idx subsys_size = subsys.size();
        for (idx i = 0; i < subsys_size; ++i) {
            if (measured_.get(subsys[i]) == true)
                throw exception::CustomException(
                    "qpp::QCirc::measure()", "Subsystem was measured before");
        }
        // update subsystem labels
        std::vector<idx> subsys_updated = update_subsys_(subsys);
        for (idx i = 0; i < subsys_size; ++i) {
            measured_.set(subsys[i]);
        }
        //        std::cout << disp(subsys, " ") << "\n";
        //        std::cout << disp(subsys_updated, " ") << "\n";
        //        std::cout << measured_ << "\n";
        auto m = measure_seq(psi_, subsys_updated);
        auto result = std::get<0>(m); // measurement result
        for (idx i = 0; i < subsys_size; ++i)
            results_[subsys[i]] = result[i];
        // update psi_
        psi_ = std::get<2>(m);
        // std::cout << disp(psi_) << "\n";
    }

    // destructive measurement of all remaining qubits
    void measure_all() {
        std::vector<idx> subsys;
        for (idx i = 0; i < nq_; ++i) {
            if (!measured_.get(i))
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

    std::vector<idx> results() const { return results_; }

    // returns the current result of the measurement as an integer,
    // ignores non-measured qubits
    idx results_as_int() const {
        std::vector<idx> tmp;
        for (idx i = 0; i < nq_; ++i) {
            if (measured_.get(i))
                tmp.push_back(results_[i]);
        }

        return multiidx2n(tmp, std::vector<idx>(tmp.size(), 2));
    }

    void apply(const cmat& gate, const std::vector<idx>& subsys) {
        psi_ = qpp::apply(psi_, gate, update_subsys_(subsys));
    }

    void applyCTRL(const cmat& gate, const std::vector<idx>& ctrl,
                   const std::vector<idx>& target) {
        psi_ = qpp::applyCTRL(psi_, gate, ctrl, target);
    }

    void apply_all(const cmat& gate) {
        for (idx i = 0; i < nq_; ++i)
            if (!measured_.get(i)) {
                // std::cout << "HERE" << disp(transpose(psi_)) << "\n";
                // std::cout << disp(update_subsys_({i}), " ") << "\n";
                psi_ = qpp::apply(psi_, gate, update_subsys_({i}));
            }
    }

    // total number of qubits, regardless of being measured or not
    idx size() const noexcept
    {
        return nq_;
    }

    // total number of measured qubits
    idx num_measured_qubits() const noexcept{
        return measured_.count();
    }

    // total number of non-measured qubits
    idx num_active_qubits() const noexcept{
        return this->size() - this->num_measured_qubits();
    }

};

int main() {
    QCirc<int> qc(10, 10);
    qc.apply_all(gt.H);
    qc.measure({3, 1, 7});
    std::cout << qc.num_measured_qubits() << std::endl;
    std::cout << qc.num_active_qubits() << std::endl;
    qc.measure({2, 4, 5, 6, 0});
    std::cout << disp(qc.results(), " ") << "\n";
    qc.apply_all(gt.H);
    qc.apply_all(gt.X);
    qc.measure_all();
    //qc.measure_all();
    std::cout << disp(qc.results(), "") << "\n";
    std::cout << qc.results_as_int() << "\n";
    std::cout << qc.size() << "\n";
    std::cout << qc.num_measured_qubits() << std::endl;
    std::cout << qc.num_active_qubits() << std::endl;
}
