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

TEST_CASE("qpp::apply() state vector benchmark",
          "[benchmark][apply-state-vector]") {
    // Setup (NOT measured)
    REQUIRE(nq > 3);
    qpp::idx D = qpp::internal::safe_pow<qpp::idx>(2, nq);
    qpp::ket psi = qpp::randket(D); // random D x 1 state vector
    qpp::cmat U1 = qpp::randU(2);   // random 1 qubit gate
    qpp::cmat U2 = qpp::randU(4);   // random 2 qubit gate
    qpp::cmat U3 = qpp::randU(8);   // random 3 qubit gate
    qpp::cmat U4 = qpp::randU(16);  // random 4 qubit gate

    // Benchmarked portion (executed repeatedly)
    BENCHMARK("Apply gate to 1 subsystem (psi) nq=" + std::to_string(nq)) {
        // CRITICAL: Return the result so the compiler doesn't optimize the
        // calculation away.
        return qpp::apply(psi, U1, {nq - 1});
    };
    // Benchmarked portion (executed repeatedly)
    BENCHMARK("Apply gate to 2 subsystems (psi) nq=" + std::to_string(nq)) {
        // CRITICAL: Return the result so the compiler doesn't optimize the
        // calculation away.
        return qpp::apply(psi, U2, {nq - 2, nq - 1});
    };
    // Benchmarked portion (executed repeatedly)
    BENCHMARK("Apply gate to 3 subsystems (psi) nq=" + std::to_string(nq)) {
        // CRITICAL: Return the result so the compiler doesn't optimize the
        // calculation away.
        return qpp::apply(psi, U3, {nq - 3, nq - 2, nq - 1});
    };
    // Benchmarked portion (executed repeatedly)
    BENCHMARK("Apply gate to 4 subsystems (psi) nq=" + std::to_string(nq)) {
        // CRITICAL: Return the result so the compiler doesn't optimize the
        // calculation away.
        return qpp::apply(psi, U4, {nq - 4, nq - 3, nq - 2, nq - 1});
    };

    // Benchmarked portion (executed repeatedly)
    BENCHMARK("Apply gate to 1 subsystem qubit optimizations (psi) nq=" +
              std::to_string(nq)) {
        // CRITICAL: Return the result so the compiler doesn't optimize the
        // calculation away.
        return qpp::internal::kernels::qubit::apply_psi_1q(psi, U1, nq - 1, nq);
    };
    // Benchmarked portion (executed repeatedly)
    BENCHMARK("Apply gate to 2 subsystems qubit optimizations (psi) nq=" +
              std::to_string(nq)) {
        // CRITICAL: Return the result so the compiler doesn't optimize the
        // calculation away.
        return qpp::internal::kernels::qubit::apply_psi_2q(psi, U2, nq - 2,
                                                           nq - 1, nq);
    };
    // Benchmarked portion (executed repeatedly)
    BENCHMARK("Apply gate to 3 subsystems qubit optimizations (psi) nq=" +
              std::to_string(nq)) {
        // CRITICAL: Return the result so the compiler doesn't optimize the
        // calculation away.
        return qpp::internal::kernels::qubit::apply_psi_3q(psi, U3, nq - 3,
                                                           nq - 2, nq - 1, nq);
    };
    // Benchmarked portion (executed repeatedly)
    BENCHMARK("Apply gate to 4 subsystems qubit optimizations (psi) nq=" +
              std::to_string(nq)) {
        // CRITICAL: Return the result so the compiler doesn't optimize the
        // calculation away.
        return qpp::internal::kernels::qubit::apply_psi_kq(
            psi, U4, {nq - 4, nq - 3, nq - 2, nq - 1}, nq);
    };
}
