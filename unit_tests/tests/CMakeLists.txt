# Quantum++
#
# Copyright (c) 2013 - 2019 Vlad Gheorghiu (vgheorgh@gmail.com)
#
# This file is part of Quantum++.
#
# Quantum++ is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Quantum++ is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Quantum++.  If not, see <http://www.gnu.org/licenses/>.

INCLUDE_DIRECTORIES(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})
AUX_SOURCE_DIRECTORY(. TEST_FILES)
SET(TARGET_NAME "qpp_testing")
ADD_EXECUTABLE(${TARGET_NAME} testing_main.cpp)
CMAKE_POLICY(SET CMP0076 NEW)

#### Build all tests in ${TEST_FILES}
FOREACH (FILE ${TEST_FILES})
    #### MATALB support
    IF (${BUILD_WITH_MATLAB})
        TARGET_SOURCES(${TARGET_NAME} PUBLIC ${FILE})
        TARGET_LINK_LIBRARIES(${TARGET_NAME} mx mat)
    ELSE () #### No MATLAB support
        #### Do not build tests/input_output.cpp if there's no MATLAB support
        IF (${FILE} STREQUAL "tests/input_output.cpp")
            CONTINUE()
        ELSE ()
            TARGET_SOURCES(${TARGET_NAME} PUBLIC ${FILE})
        ENDIF ()
    ENDIF ()
    #### Eigen3 was found automatically
    IF (TARGET Eigen3::Eigen)
        TARGET_LINK_LIBRARIES(${TARGET_NAME} Eigen3::Eigen)
    ENDIF ()
    #### OpenMP was enabled
    IF ($WITH_OPENMP$ AND ${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang"
            AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER "3.7")
        TARGET_LINK_LIBRARIES(${TARGET_NAME} omp)
    ENDIF ()
    TARGET_LINK_LIBRARIES(${TARGET_NAME} gmock)
ENDFOREACH ()
