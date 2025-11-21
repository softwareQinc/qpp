#ifndef OPENMP_UTILS_HPP_
#define OPENMP_UTILS_HPP_

#ifdef QPP_OPENMP
#include <cstdlib>
#include <iostream>

#include <omp.h>

// CLI override (default -1 -> auto-detect)
inline int cli_core_count = -1;

/**
 * \brief Detects/sets the number of OpenMP threads
 * \param requested Optional requested thread count; uses \a cli_core_count by
 * default
 *
 * \return Final number of threads used
 */
inline int set_openmp_threads(int requested = -1) {
    int core_count = (requested >= 0) ? requested : cli_core_count;

    const char* env = std::getenv("OMP_NUM_THREADS");
    const bool env_set = (env && std::atoi(env) > 0);

    // 1. Environment variable has highest priority
    if (env_set) {
        core_count = std::atoi(env);
        std::cout << "OpenMP enabled, core_count=" << core_count
                  << " set by OMP_NUM_THREADS\n";
    }
    // 2. If NOT set by OMP_NUM_THREADS
    else if (core_count <= 0) {
        core_count = omp_get_max_threads();
        std::cout << "OpenMP enabled, core_count=" << core_count
                  << " auto-detected by OpenMP\n";
    }
    // 3. Otherwise (core_count > 0) and NOT set by OMP_NUM_THREADS
    else {
        std::cout << "OpenMP enabled, core_count=" << core_count
                  << " set by the CLI\n";
    }

    omp_set_num_threads(core_count);
    return core_count;
}
#endif // QPP_OPENMP

#endif // OPENMP_UTILS_HPP_
