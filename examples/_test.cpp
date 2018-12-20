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
                if (measured_.get(m)) // if the qubit was measured
                {
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
        std::vector<idx> subsys_updated = subsys;
        for (idx i = 0; i < subsys_size; ++i) {
            if (measured_.get(subsys[i]) == true)
                throw exception::CustomException(
                    "qpp::QCirc::measure()", "Subsystem was measured before");
        }
        // update subsystem labels
        subsys_updated = update_subsys_(subsys);
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

    std::vector<idx> results() const { return results_; }

    // returns the current result of the measurement as an integer
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

    void apply_all(const cmat& gate) {
        for (idx i = 0; i < nq_; ++i)
            if (!measured_.get(i)) {
                // std::cout << "HERE" << disp(transpose(psi_)) << "\n";
                // std::cout << disp(update_subsys_({i}), " ") << "\n";
                psi_ = qpp::apply(psi_, gate, update_subsys_({i}));
            }
    }
};

int main() {
    QCirc<int> qc(10, 10);
    qc.apply_all(gt.H);
    qc.measure({3, 1, 7});
    qc.measure({2, 4, 5, 6, 0, 9});
    std::cout << disp(qc.results(), " ") << "\n";
    qc.apply_all(gt.H);
    qc.apply_all(gt.X);
    qc.measure({8});
    std::cout << disp(qc.results(), "") << "\n";
    std::cout << qc.results_as_int() << "\n";
}
