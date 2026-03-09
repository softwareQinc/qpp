#define CATCH_CONFIG_RUNNER

#include <string>

#include <catch2/catch_all.hpp>
#include <qpp/qpp.hpp>

#include "openmp_utils.hpp"

namespace {
qpp::idx nq = 10; // default number of qubits if none provided at runtime
} // namespace

int main(int argc, char* argv[]) {
    Catch::Session session;

    auto cli = session.cli() |
               Catch::Clara::Opt(nq, "qubits")["--nq"]("Number of qubits");
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

TEST_CASE("qpp::measure_seq() state vector benchmark",
          "[benchmark][measure_seq-state-vector]") {
    // Setup (NOT measured)
    REQUIRE(nq > 0);
    qpp::idx D = qpp::internal::safe_pow<qpp::idx>(2, nq);
    qpp::ket psi = qpp::randket(D); // random D x 1 state vector
    // Test worst-case scenario
    std::vector<qpp::idx> subsys(nq);
    std::iota(subsys.begin(), subsys.end(), 0);

    // Benchmarked portion (executed repeatedly)
    BENCHMARK("measure_seq/baseline/psi/nq=" + std::to_string(nq)) {
        // CRITICAL: Return the result so the compiler doesn't optimize the
        // calculation away.
        return qpp::measure_seq(psi, subsys);
    };
}
