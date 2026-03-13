#define CATCH_CONFIG_RUNNER

#include <string>

#include <catch2/catch_all.hpp>

#include <qpp/internal/kernels/qubit/apply_ctrl.hpp>
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

TEST_CASE("qpp::applyCTRL() state vector benchmark",
          "[benchmark][applyCTRL-state-vector]") {
    // Setup (NOT measured)
    REQUIRE(nq > 4);
    qpp::idx D = qpp::internal::safe_pow<qpp::idx>(2, nq);
    qpp::ket psi = qpp::randket(D); // random D x 1 state vector
    qpp::cmat U1 = qpp::randU(2);   // random 1 qubit gate
    qpp::cmat U2 = qpp::randU(4);   // random 2 qubit gate
    qpp::cmat U3 = qpp::randU(8);   // random 3 qubit gate

    // Benchmarked portion (executed repeatedly)
    BENCHMARK("applyCTRL/baseline/psi/nq=" + std::to_string(nq) +
              "/ctrls=1/targets=1") {
        // CRITICAL: Return the result so the compiler doesn't optimize the
        // calculation away.
        return qpp::applyCTRL(psi, U1, {nq - 1}, {nq - 2});
    };
    // Benchmarked portion (executed repeatedly)
    BENCHMARK("applyCTRL/baseline/psi/nq=" + std::to_string(nq) +
              "/ctrls=1/targets=2") {
        // CRITICAL: Return the result so the compiler doesn't optimize the
        // calculation away.
        return qpp::applyCTRL(psi, U2, {nq - 1}, {nq - 2, nq - 3});
    };
    // Benchmarked portion (executed repeatedly)
    BENCHMARK("applyCTRL/baseline/psi/nq=" + std::to_string(nq) +
              "/ctrls=1/targets=3") {
        // CRITICAL: Return the result so the compiler doesn't optimize the
        // calculation away.
        return qpp::applyCTRL(psi, U3, {nq - 1}, {nq - 2, nq - 3, nq - 4});
    };
    // Benchmarked portion (executed repeatedly)
    BENCHMARK("applyCTRL/baseline/psi/nq=" + std::to_string(nq) +
              "/ctrls=2/targets=1") {
        // CRITICAL: Return the result so the compiler doesn't optimize the
        // calculation away.
        return qpp::applyCTRL(psi, U1, {nq - 1, nq - 2}, {nq - 3});
    };
    // Benchmarked portion (executed repeatedly)
    BENCHMARK("applyCTRL/baseline/psi/nq=" + std::to_string(nq) +
              "/ctrls=2/targets=2") {
        // CRITICAL: Return the result so the compiler doesn't optimize the
        // calculation away.
        return qpp::applyCTRL(psi, U2, {nq - 1, nq - 2}, {nq - 3, nq - 4});
    };
    // Benchmarked portion (executed repeatedly)
    BENCHMARK("applyCTRL/baseline/psi/nq=" + std::to_string(nq) +
              "/ctrls=2/targets=3") {
        // CRITICAL: Return the result so the compiler doesn't optimize the
        // calculation away.
        return qpp::applyCTRL(psi, U3, {nq - 1, nq - 2},
                              {nq - 3, nq - 4, nq - 5});
    };

    // Benchmarked portion (executed repeatedly)
    BENCHMARK("applyCTRL/qubit-kernel/psi/nq=" + std::to_string(nq) +
              "/ctrls=1/targets=1") {
        // CRITICAL: Return the result so the compiler doesn't optimize the
        // calculation away.
        return qpp::internal::kernels::qubit::apply_ctrl_psi_1q(
            psi, U1, {nq - 1}, nq - 2, {0}, nq);
    };
    // Benchmarked portion (executed repeatedly)
    BENCHMARK("applyCTRL/qubit-kernel/psi/nq=" + std::to_string(nq) +
              "/ctrls=1/targets=2") {
        // CRITICAL: Return the result so the compiler doesn't optimize the
        // calculation away.
        return qpp::internal::kernels::qubit::apply_ctrl_psi_2q(
            psi, U1, {nq - 1}, nq - 2, nq - 3, {0}, nq);
    };
    // Benchmarked portion (executed repeatedly)
    BENCHMARK("applyCTRL/qubit-kernel/psi/nq=" + std::to_string(nq) +
              "/ctrls=1/targets=3") {
        // CRITICAL: Return the result so the compiler doesn't optimize the
        // calculation away.
        return qpp::internal::kernels::qubit::apply_ctrl_psi_kq(
            psi, U1, {nq - 1}, {nq - 2, nq - 3, nq - 4}, {0}, nq);
    };
    // Benchmarked portion (executed repeatedly)
    BENCHMARK("applyCTRL/qubit-kernel/psi/nq=" + std::to_string(nq) +
              "/ctrls=2/targets=1") {
        // CRITICAL: Return the result so the compiler doesn't optimize the
        // calculation away.
        return qpp::internal::kernels::qubit::apply_ctrl_psi_1q(
            psi, U1, {nq - 1, nq - 2}, nq - 3, {0, 0}, nq);
    };
    // Benchmarked portion (executed repeatedly)
    BENCHMARK("applyCTRL/qubit-kernel/psi/nq=" + std::to_string(nq) +
              "/ctrls=2/targets=2") {
        // CRITICAL: Return the result so the compiler doesn't optimize the
        // calculation away.
        return qpp::internal::kernels::qubit::apply_ctrl_psi_2q(
            psi, U1, {nq - 1, nq - 2}, nq - 3, nq - 4, {0, 0}, nq);
    };
    // Benchmarked portion (executed repeatedly)
    BENCHMARK("applyCTRL/qubit-kernel/psi/nq=" + std::to_string(nq) +
              "/ctrls=2/targets=3") {
        // CRITICAL: Return the result so the compiler doesn't optimize the
        // calculation away.
        return qpp::internal::kernels::qubit::apply_ctrl_psi_kq(
            psi, U1, {nq - 1, nq - 2}, {nq - 3, nq - 4, nq - 5}, {0, 0}, nq);
    };
}
