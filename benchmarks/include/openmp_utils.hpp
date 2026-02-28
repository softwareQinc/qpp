#ifndef OPENMP_UTILS_HPP_
#define OPENMP_UTILS_HPP_

#ifdef QPP_OPENMP

#include <cstdlib>
#include <iostream>

#include <omp.h>

// CLI override (default -1 -> auto-detect)
inline int cli_threads = -1;

/**
 * @brief Detects/sets the number of OpenMP threads
 * @param requested Optional requested thread count; uses \a cli_threads by
 * default
 *
 * @return Final number of threads used
 */
inline int set_openmp_threads(int requested = -1) {
    int threads = (requested >= 0) ? requested : cli_threads;

    const char* env = std::getenv("OMP_NUM_THREADS");
    const bool env_set = (env && std::atoi(env) > 0);

    // 1. Environment variable has highest priority
    if (env_set) {
        threads = std::atoi(env);
        std::cout << "OpenMP enabled, threads=" << threads
                  << " set by OMP_NUM_THREADS\n";
    }
    // 2. If NOT set by OMP_NUM_THREADS
    else if (threads <= 0) {
        threads = omp_get_max_threads();
        std::cout << "OpenMP enabled, threads=" << threads
                  << " auto-detected by OpenMP\n";
    }
    // 3. Otherwise (threads > 0) and NOT set by OMP_NUM_THREADS
    else {
        std::cout << "OpenMP enabled, threads=" << threads
                  << " set by the CLI\n";
    }

    omp_set_num_threads(threads);

    return threads;
}
#endif // QPP_OPENMP

#endif // OPENMP_UTILS_HPP_
