#
# This file is part of Quantum++.
#
# MIT License
#
# Copyright (c) 2013 - 2019 Vlad Gheorghiu (vgheorgh@gmail.com)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#

INCLUDE_DIRECTORIES(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})
AUX_SOURCE_DIRECTORY(. TEST_FILES)
AUX_SOURCE_DIRECTORY(MATLAB TEST_FILES)
AUX_SOURCE_DIRECTORY(classes TEST_FILES)
SET(TARGET_NAME "qpp_testing")
ADD_EXECUTABLE(${TARGET_NAME} main.cpp)
IF (NOT ${CMAKE_VERSION} VERSION_LESS "3.13")
    CMAKE_POLICY(SET CMP0076 NEW)
ENDIF ()

#### Build all tests in ${TEST_FILES}
FOREACH (FILE ${TEST_FILES})
    #### MATALB support
    IF (${BUILD_WITH_MATLAB})
        TARGET_SOURCES(${TARGET_NAME} PUBLIC ${FILE})
        TARGET_LINK_LIBRARIES(${TARGET_NAME} mx mat)
    ELSE () #### No MATLAB support
        #### Do not build tests/matlab.cpp if there's no MATLAB support
        IF (${FILE} STREQUAL "MATLAB/matlab.cpp")
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
