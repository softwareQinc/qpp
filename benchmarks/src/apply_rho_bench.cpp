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

TEST_CASE("qpp::apply() density matrix benchmark",
          "[benchmark][apply-density-matrix]") {
    // Setup (NOT measured)
    REQUIRE(nq > 3);
    qpp::idx D = qpp::internal::safe_pow<qpp::idx>(2, nq);
    qpp::cmat rho = qpp::randrho(D); // random D x D density matrix
    qpp::cmat U1 = qpp::randU(2);    // random 1 qubit gate
    qpp::cmat U2 = qpp::randU(4);    // random 2 qubit gate
    qpp::cmat U3 = qpp::randU(8);    // random 3 qubit gate
    qpp::cmat U4 = qpp::randU(16);   // random 4 qubit gate

    // Benchmarked portion (executed repeatedly)
    BENCHMARK("apply/baseline/rho/nq=" + std::to_string(nq) + "/targets=1") {
        // CRITICAL: Return the result so the compiler doesn't optimize the
        // calculation away.
        return qpp::apply(rho, U1, {nq - 1});
    };
    // Benchmarked portion (executed repeatedly)
    BENCHMARK("apply/baseline/rho/nq=" + std::to_string(nq) + "/targets=2") {
        // CRITICAL: Return the result so the compiler doesn't optimize the
        // calculation away.
        return qpp::apply(rho, U2, {nq - 2, nq - 1});
    };
    // Benchmarked portion (executed repeatedly)
    BENCHMARK("apply/baseline/rho/nq=" + std::to_string(nq) + "/targets=3") {
        // CRITICAL: Return the result so the compiler doesn't optimize the
        // calculation away.
        return qpp::apply(rho, U3, {nq - 3, nq - 2, nq - 1});
    };
    // Benchmarked portion (executed repeatedly)
    BENCHMARK("apply/baseline/rho/nq=" + std::to_string(nq) + "/targets=4") {
        // CRITICAL: Return the result so the compiler doesn't optimize the
        // calculation away.
        return qpp::apply(rho, U4, {nq - 4, nq - 3, nq - 2, nq - 1});
    };

    // Benchmarked portion (executed repeatedly)
    BENCHMARK("apply/qubit-kernel/rho/nq=" + std::to_string(nq) +
              "/targets=1") {
        // CRITICAL: Return the result so the compiler doesn't optimize the
        // calculation away.
        return qpp::internal::kernels::qubit::apply_rho_1q(rho, U1, nq - 1, nq);
    };
    // Benchmarked portion (executed repeatedly)
    BENCHMARK("apply/qubit-kernel/rho/nq=" + std::to_string(nq) +
              "/targets=2") {
        // CRITICAL: Return the result so the compiler doesn't optimize the
        // calculation away.
        return qpp::internal::kernels::qubit::apply_rho_2q(rho, U2, nq - 2,
                                                           nq - 1, nq);
    };
    // Benchmarked portion (executed repeatedly)
    BENCHMARK("apply/qubit-kernel/rho/nq=" + std::to_string(nq) +
              "/targets=3") {
        // CRITICAL: Return the result so the compiler doesn't optimize the
        // calculation away.
        return qpp::internal::kernels::qubit::apply_rho_3q(rho, U3, nq - 3,
                                                           nq - 2, nq - 1, nq);
    };
    // Benchmarked portion (executed repeatedly)
    BENCHMARK("apply/qubit-kernel/rho/nq=" + std::to_string(nq) +
              "/targets=4") {
        // CRITICAL: Return the result so the compiler doesn't optimize the
        // calculation away.
        return qpp::internal::kernels::qubit::apply_rho_kq(
            rho, U4, {nq - 4, nq - 3, nq - 2, nq - 1}, nq);
    };
}
