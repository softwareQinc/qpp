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

    auto cli =
        session.cli() | Catch::Clara::Opt(nq, "nq")["--nq"]("Number of qubits");
#ifdef QPP_OPENMP
    cli |= Catch::Clara::Opt(cli_core_count, "core_count")["--core_count"](
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

TEST_CASE("qpp::ptrace() density matrix benchmark",
          "[benchmark][ptrace-density-matrix]") {
    // Setup (NOT measured)
    REQUIRE(nq > 0);
    qpp::idx D = qpp::internal::safe_pow<qpp::idx>(2, nq);
    qpp::cmat rho = qpp::randrho(D); // D x D random density matrix
    // Test worst-case scenario
    std::vector<qpp::idx> subsys = {nq - 1};

    // Benchmarked portion (executed repeatedly)
    BENCHMARK("Partial trace (rho) nq=" + std::to_string(nq)) {
        // CRITICAL: Return the result so the compiler doesn't optimize the
        // calculation away.
        return qpp::ptrace(rho, subsys);
    };
}
