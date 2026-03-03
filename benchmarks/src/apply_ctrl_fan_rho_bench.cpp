#define CATCH_CONFIG_RUNNER

#include <string>

#include <catch2/catch_all.hpp>

#include <qpp/internal/kernels/qubit/apply_ctrl_fan.hpp>
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

TEST_CASE("qpp::applyCTRL_fan() density matrix benchmark",
          "[benchmark][applyCTRL_fan-density-matrix]") {
    // Setup (NOT measured)
    REQUIRE(nq > 1);
    qpp::idx D = qpp::internal::safe_pow<qpp::idx>(2, nq);
    qpp::cmat rho = qpp::randrho(D); // D x D random density matrix
    // Test worst-case scenario
    qpp::idx mid = nq / 2;
    std::vector<qpp::idx> ctrl(mid);
    std::vector<qpp::idx> target(mid);
    std::vector<qpp::idx> shift(nq, 0);
    for (qpp::idx i = 0; i < mid; ++i) {
        ctrl[i] = i;
        target[i] = mid + i;
    }

    // Benchmarked portion (executed repeatedly)
    BENCHMARK("CTRL-fan (rho) nq=" + std::to_string(nq)) {
        // CRITICAL: Return the result so the compiler doesn't optimize the
        // calculation away.
        return qpp::applyCTRL_fan(rho, qpp::gt.H, ctrl, target);
    };

    // Benchmarked portion (executed repeatedly)
    BENCHMARK("CTRL-fan qubit optimizations (rho) nq=" + std::to_string(nq)) {
        // CRITICAL: Return the result so the compiler doesn't optimize the
        // calculation away.
        return qpp::internal::kernels::qubit::apply_ctrl_fan_rho(
            rho, qpp::gt.H, ctrl, target, shift, nq);
    };
}
