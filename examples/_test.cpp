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
  public:
    QCirc(idx nq, idx nc)
        : nq_{nq}, nc_{nc}, psi_{st.zero(nq_)}, measured_{nq_} {}

    // destructive measurement
    void measure(std::vector<idx> subsys) {
        // sort subsys in decreasing order
        idx subsys_size = subsys.size();
        std::sort(std::begin(subsys), std::end(subsys), std::less<idx>{});
        std::vector<idx> subsys_copy = subsys;
        for (idx i = 0; i < subsys_size; ++i) {
            if (measured_.get(subsys[i]) == true)
                throw exception::CustomException(
                    "qpp::QCirc::measure()", "Subsystem was measured before");
            measured_.set(subsys[i], true);
            // re-label subsystems according to what was measured before
            for (idx m = 0; m < subsys[i]; ++m) {
                if (measured_.get(m) == true)
                    --subsys_copy[i];
            }
        }
        auto m = measure_seq(psi_, subsys_copy);
        auto result = std::get<0>(m); // measurement result
        for(idx i = 0; i < subsys_size; ++i)
            results_[subsys[i]] = result[i];
        // update psi_
        psi_ = std::get<2>(m);
    }
};

int main() {
    std::cout << "ok\n";
    QCirc<int> qc(10, 10);
}
