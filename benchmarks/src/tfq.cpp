#define CATCH_CONFIG_RUNNER

#include <iostream>
#include <string>
#include <thread>

#include <catch2/catch_all.hpp>

#ifdef QPP_OPENMP
#include <omp.h>
#endif

#include <qpp/qpp.hpp>

namespace {
qpp::idx nq = 10; // default number of qubits if none provided at runtime
#ifdef QPP_OPENMP
int num_cores = -1; // -1 -> auto-detect CPU core count
#endif
} // namespace

int main(int argc, char* argv[]) {
    Catch::Session session;
    auto cli = session.cli() |
               Catch::Clara::Opt(nq, "nq")["--nq"]("Number of qubits")
#ifdef QPP_OPENMP
               | Catch::Clara::Opt(num_cores, "num_cores")["--num_cores"](
                     "Number of OpenMP threads/cores")
#endif
        ; // keep the ; on this line!
    session.cli(cli);
    session.applyCommandLine(argc, argv);

#ifdef QPP_OPENMP
    if (num_cores <= 0) {
        // Fallback: Use standard C++ to find hardware threads
        // If hardware_concurrency returns 0 (fail), fallback to 1 or 2 safely
        unsigned int hw_threads = std::thread::hardware_concurrency();
        num_cores = (hw_threads > 0) ? static_cast<int>(hw_threads) : 1;

        std::cout << "OpenMP enabled, num_cores=" << num_cores
                  << " auto-detected\n";
    } else {
        std::cout << "OpenMP enabled, num_cores=" << num_cores
                  << " explicitly set\n";
    }
    omp_set_num_threads(num_cores);
#endif

    int result = session.run();
    return result;
}

TEST_CASE("qpp::QCircuit::TFQ() benchmark", "[benchmark][tfq]") {
    REQUIRE(nq > 0);
    BENCHMARK_ADVANCED("TFQ nq=" + std::to_string(nq))
    (Catch::Benchmark::Chronometer meter) {
        // Setup (NOT measured)
        qpp::QCircuit qc{nq};
        qc.TFQ();
        qpp::QEngine qe{qc};

        // Benchmarked portion (executed repeatedly)
        meter.measure([&] {
            qe.reset(); // important: engine must be reset after each iteration
            qe.execute(); // <- measured
        });
    };
}
