#ifndef OPENMP_UTILS_HPP_
#define OPENMP_UTILS_HPP_

#ifdef QPP_OPENMP

#include <cstdlib>
#include <iostream>

#include <omp.h>

// CLI override (default -1 -> auto-detect)
inline int cli_threads = -1;

/**
 * @brief Detects/sets the number of OpenMP threads.
 * Priority: 1. CLI (--nq) -> 2. OMP_NUM_THREADS -> 3. System Default
 *
 * @param requested Optional requested thread count; uses cli_threads by
 * default
 * @return Final number of threads used
 */
inline int set_openmp_threads(int requested = -1) {
    int threads = (requested >= 0) ? requested : cli_threads;

    // 1. CLI / Requested parameter has HIGHEST priority
    if (threads > 0) {
        std::cout << "OpenMP enabled, threads=" << threads
                  << " set by the CLI\n";
    } else {
        // 2. Environment variable is the second priority
        const char* env = std::getenv("OMP_NUM_THREADS");
        if (env && std::atoi(env) > 0) {
            threads = std::atoi(env);
            std::cout << "OpenMP enabled, threads=" << threads
                      << " set by the environment variable OMP_NUM_THREADS\n";
        }
        // 3. Fallback to OpenMP auto-detection
        else {
            threads = omp_get_max_threads();
            std::cout << "OpenMP enabled, threads=" << threads
                      << " auto-detected by OpenMP\n";
        }
    }
    omp_set_num_threads(threads);

    return threads;
}

#endif // QPP_OPENMP
#endif // OPENMP_UTILS_HPP_
