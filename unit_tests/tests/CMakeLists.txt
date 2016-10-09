INCLUDE_DIRECTORIES(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})
ADD_EXECUTABLE(qpp_testing
        classes/gates.cpp
        classes/timer.cpp
        entanglement.cpp
        entropies.cpp
        functions.cpp
        input_output.cpp
        instruments.cpp
        number_theory.cpp
        operations.cpp
        random.cpp
        statistics.cpp
        testing_main.cpp
        traits.cpp)
TARGET_LINK_LIBRARIES(qpp_testing gmock pthread)
