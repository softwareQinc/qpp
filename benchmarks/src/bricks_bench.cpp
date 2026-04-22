#define CATCH_CONFIG_RUNNER

#include <string>

#include <catch2/catch_all.hpp>
#include <qpp/qpp.hpp>

#include "openmp_utils.hpp"

namespace {
qpp::idx nq = 10; // default number of qubits if none provided at runtime
qpp::idx nl = 10; // default number of layers if none provided at runtime
} // namespace

int main(int argc, char* argv[]) {
    Catch::Session session;

    auto cli = session.cli() |
               Catch::Clara::Opt(nq, "qubits")["--nq"]("Number of qubits") |
               Catch::Clara::Opt(nl, "layers")["--nl"]("Number of layers");
#ifdef QPP_OPENMP
    cli |= Catch::Clara::Opt(cli_threads, "threads")["--threads"](
        "Number of OpenMP threads/cores (ignored if the environment variable "
        "OMP_NUM_THREADS is set)");
#endif // QPP_OPENMP
    session.cli(cli);

    int returnCode = session.applyCommandLine(argc, argv);
    if (returnCode != 0) {
        return returnCode;
    }

    if (session.config().showHelp()) {
        return 0;
    }

#ifdef QPP_OPENMP
    set_openmp_threads();
#endif // QPP_OPENMP

    return session.run();
}

TEST_CASE("Bricks circuit benchmark", "[benchmark][bricks]") {
    // Setup (NOT measured)
    REQUIRE(nq > 1);
    qpp::QCircuit qc{nq};
    for (qpp::idx i = 0; i < nl; ++i) {
        // layer of 1-qubit random unitaries
        for (qpp::idx q = 0; q < nq; ++q) {
            qc.gate(qpp::randU(), q);
        }
        // followed by a layer of 2-qubit random unitaries
        for (qpp::idx q = 0; q < nq - 1; q += 2) {
            qc.gate(qpp::randU(4), q, q + 1);
        }
    }
    qpp::QEngine qe{qc};

    // Benchmarked portion (executed repeatedly)
    BENCHMARK("Bricks circuit/baseline/psi/nq=" + std::to_string(nq) +
              "/nl=" + std::to_string(nl)) {
        // IMPORTANT: engine must be reset after each iteration
        qe.reset();
        // CRITICAL: Return the result so the compiler doesn't optimize the
        // calculation away.
        return qe.execute();
    };
}
