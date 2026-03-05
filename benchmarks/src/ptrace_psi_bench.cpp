#define CATCH_CONFIG_RUNNER

#include <string>

#include <catch2/catch_all.hpp>
#include <qpp/qpp.hpp>

#include "openmp_utils.hpp"

#include "qpp/internal/kernels/qubit/ptrace.hpp"

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

TEST_CASE("qpp::ptrace() state vector benchmark",
          "[benchmark][ptrace-state-vector]") {
    // Setup (NOT measured)
    REQUIRE(nq > 0);
    qpp::idx D = qpp::internal::safe_pow<qpp::idx>(2, nq);
    qpp::ket psi = qpp::randket(D); // random D x 1 state vector
    // Test worst-case scenario
    std::vector<qpp::idx> subsys = {nq - 1};

    std::vector<qpp::idx> dims(nq, 2);

    // Benchmarked portion (executed repeatedly)
    BENCHMARK("Partial trace (psi) nq=" + std::to_string(nq)) {
        // CRITICAL: Return the result so the compiler doesn't optimize the
        // calculation away.
        return qpp::ptrace(psi, subsys);
    };

    // Benchmarked portion (executed repeatedly)
    BENCHMARK("Partial trace qubit optimizations (psi) nq=" +
              std::to_string(nq)) {
        // CRITICAL: Return the result so the compiler doesn't optimize the
        // calculation away.
        return qpp::internal::kernels::qubit::ptrace_psi_kq(psi, subsys, nq);
    };
}
